"""Replay recorded scores through every manifest stage; never launch Boltz/GPU work."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest
import yaml
from biosim.signals import make_signal, unwrap_payload
from biosim.world import BioWorld

LAB = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("failed", [set(), {"Dasatinib"}, {"Dasatinib", "Imatinib", "Nilotinib"}])
def test_manifest_graph_keeps_outputs_and_failure_rows(monkeypatch, tmp_path, failed):
    manifest = yaml.safe_load((LAB / "lab.yaml").read_text())
    fixture = json.loads((LAB / "tests/fixtures/recorded_batch_scores.json").read_text())
    scores = {row["ligand"]: row for row in fixture["rows"]}
    calls = []

    class RecordedRunner:
        def __init__(self, **kwargs):
            self.name = kwargs["default_run_options"]["batch_ligand_name"]
            assert kwargs["default_ligand_smiles"] == scores[self.name]["smiles"]
            assert len(kwargs["default_protein_sequence"]) == 273
        def execute(self, inputs, *, context):
            calls.append(self.name)
            row = scores[self.name]
            data = {
                "run_metadata": {"status": "error" if self.name in failed else "completed", "error": "simulated failure" if self.name in failed else None},
                "affinity_summary": {"affinity_pred_value": row["affinity_like_value"], "affinity_probability_binary": row["binder_probability"]},
                "confidence_summary": {"confidence_score": row["confidence"]},
                # Only a path-routing fixture: no synthetic coordinates are supplied.
                "structure_artifacts": {"structure_file": str(tmp_path / (self.name + ".cif"))},
            }
            return {key: make_signal(source="recorded_score_replay", name=key, value=value, emitted_at=0., spec=None) for key, value in data.items()}

    monkeypatch.syspath_prepend(str(LAB / "models/core"))
    for name in list(sys.modules):
        if name == "src" or name.startswith("src."):
            monkeypatch.delitem(sys.modules, name)
    world = BioWorld(communication_step=manifest["runtime"]["communication_step"])
    instances = {}
    for entry in manifest["models"]:
        directory = LAB / entry["path"]
        model_manifest = yaml.safe_load((directory / "model.yaml").read_text())
        module_name, class_name = model_manifest["biosim"]["entrypoint"].split(":")
        spec = importlib.util.spec_from_file_location("replay_" + entry["alias"], directory / (module_name.replace(".", "/") + ".py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if entry["alias"] == "boltz_batch_ligand_ranker":
            monkeypatch.setattr(module, "Boltz2AffinityPredictor", RecordedRunner)
        kwargs = dict(model_manifest["biosim"].get("init_kwargs", {}))
        kwargs.update(entry.get("parameters", {}))
        instance = getattr(module, class_name)(**kwargs)
        instances[entry["alias"]] = instance
        world.add_biomodule(entry["alias"], instance)
    for wire in manifest["wiring"]:
        for target in wire["to"]:
            world.connect(wire["from"], target)
    world.run(manifest["runtime"]["duration"])
    assert calls == ["Imatinib", "Dasatinib", "Nilotinib"]
    batch = unwrap_payload(world.get_outputs("boltz_batch_ligand_ranker")["batch_summary"])
    evidence = unwrap_payload(world.get_outputs("ranking_interpreter")["prediction_evidence"])
    expected = [name for name in ["Dasatinib", "Imatinib", "Nilotinib"] if name not in failed]
    assert [row["ligand"] for row in batch["ranked_ligands"] if row["rank"]] == expected
    assert batch["top_ligand_name"] == (expected[0] if expected else None)
    assert evidence["top_ligand"] == batch["top_ligand_name"]
    visuals = instances["visualisation"].visualize()
    ranking_table = next(visual["data"] for visual in visuals if visual["render"] == "table" and visual["data"]["title"].startswith("Batch Ligand Ranking"))
    assert len(ranking_table["rows"]) == 3
    assert "Predicted log10(IC50 / µM)" in ranking_table["columns"]
    assert "No completed Boltz prediction" not in evidence["observed_answer"] if expected else "No completed Boltz prediction" in evidence["observed_answer"]
    if expected:
        structure = next(visual for visual in visuals if visual["render"] == "structure3d")
        assert structure["data"]["source"]["path"].endswith(expected[0] + ".cif")
    else:
        assert not any(visual["render"] == "structure3d" for visual in visuals)
