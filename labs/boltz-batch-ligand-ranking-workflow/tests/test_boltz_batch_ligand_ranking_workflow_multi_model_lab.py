from __future__ import annotations

from pathlib import Path

import yaml


def test_lab_manifest_uses_embedded_visualisation_model():
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())
    model_aliases = [entry["alias"] for entry in manifest["models"]]
    assert "boltz_batch_ligand_ranker" in model_aliases
    assert "visualisation" in model_aliases
    assert all(str(entry["path"]).startswith("models/") for entry in manifest["models"])

    assert not any(str(entry["maps_to"]).startswith("visualisation.") for entry in manifest["io"]["outputs"])

    wiring_targets = {target for entry in manifest["wiring"] for target in entry["to"]}
    assert any(target.startswith("visualisation.") for target in wiring_targets)


def test_lab_manifest_is_guided_boltz_batch_workflow():
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())

    assert manifest["title"] == "Boltz Workflow: Batch Ligand Ranking"
    assert "boltz-workflow" in manifest["tags"]
    assert "batch-ranking" in manifest["tags"]

    core = next(entry for entry in manifest["models"] if entry["alias"] == "boltz_batch_ligand_ranker")
    params = core["parameters"]
    assert params["default_protein_sequence"]
    assert params["default_ligand_csv"]
    assert params["max_ligands"] == 3
    assert params["default_run_options"]["workflow_name"] == "Batch Ligand Ranking"
