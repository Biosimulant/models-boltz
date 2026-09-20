from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
from biosim.signals import make_signal, unwrap_payload
from biosim.modules import ExecutionContext, ExecutionPolicy

from src.boltz2_batch_ligand_ranker import Boltz2BatchLigandRanker


def test_parse_ligand_csv_accepts_name_smiles_and_metadata(tmp_path):
    module = Boltz2BatchLigandRanker(work_dir=str(tmp_path))

    ligands = module._parse_ligand_csv(
        "name,smiles,source\n"
        "Imatinib,CCN,PubChem\n"
        "Dasatinib,CCC,PubChem\n"
    )

    assert [ligand["name"] for ligand in ligands] == ["Imatinib", "Dasatinib"]
    assert [ligand["smiles"] for ligand in ligands] == ["CCN", "CCC"]
    assert "PubChem" in ligands[0]["metadata"]


@pytest.mark.parametrize("mode,expected", [
    ("active_affinity", ["strong", "tie", "probable"]),
    ("binder_probability", ["probable", "strong", "tie"]),
])
def test_ranking_uses_task_and_lower_affinity(mode, expected, tmp_path):
    module = Boltz2BatchLigandRanker(default_run_options={"ranking_mode": mode})
    rows = [
        {"ligand": "probable", "binder_probability": 0.9, "affinity_like_value": 1., "status": "completed"},
        {"ligand": "tie", "binder_probability": 0.8, "affinity_like_value": -1., "status": "completed"},
        {"ligand": "strong", "binder_probability": 0.8, "affinity_like_value": -2., "status": "completed"},
        {"ligand": "failed", "binder_probability": 1., "affinity_like_value": -9., "status": "error"},
    ]
    ranked = module._rank_rows(rows)
    assert [row["ligand"] for row in ranked] == expected + ["failed"]
    assert [row["rank"] for row in ranked] == [1, 2, 3, None]


def test_missing_csv_surfaces_error_payload(tmp_path):
    module = Boltz2BatchLigandRanker(default_protein_sequence="MKT", work_dir=str(tmp_path))

    outputs = module.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0.0, run_end=0.1))

    metadata = module._output_payloads["run_metadata"]
    assert metadata["status"] == "error"
    assert "ligand_csv" in metadata["error"]


@pytest.mark.parametrize("field,value", [
    ("affinity_probability_binary", float("nan")),
    ("affinity_probability_binary", float("inf")),
    ("affinity_probability_binary", True),
    ("affinity_probability_binary", 1.1),
    ("affinity_probability_binary", -0.1),
    ("affinity_pred_value", float("inf")),
    ("affinity_pred_value", None),
    ("affinity_pred_value", "-2"),
])
def test_invalid_scores_cannot_receive_a_rank(field, value):
    module = Boltz2BatchLigandRanker()
    affinity = {"affinity_probability_binary": 0.8, "affinity_pred_value": -2.}
    affinity[field] = value
    row = module._build_row(1, {"name": "bad", "smiles": "CC"}, "completed", affinity, {}, {})
    assert row["status"] == "error"
    assert "invalid affinity" in row["error"]
    assert module._rank_rows([row])[0]["rank"] is None


@pytest.mark.parametrize("csv,error", [
    ("name,smiles\nA,CC\nB,\n", "no SMILES"),
    ("name,smiles\nA,CC\nB,CC\nC,CC\nD,CC\n", "exceeds max_ligands"),
])
def test_invalid_library_is_rejected_before_any_runner(csv, error, tmp_path):
    module = Boltz2BatchLigandRanker(default_protein_sequence="MKT", default_ligand_csv=csv, work_dir=str(tmp_path))
    outputs = module.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0., run_end=1.))
    metadata = unwrap_payload(outputs["run_metadata"])
    assert metadata["status"] == "error"
    assert error in metadata["error"]
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("all_failed", [False, True])
def test_batch_preserves_failures_selects_completed_top_and_interprets_partial(monkeypatch, tmp_path, all_failed):
    # Load this wrapper under a stable name: the repository clears src modules per test.
    path = Path(__file__).resolve().parents[1] / "src/boltz2_batch_ligand_ranker.py"
    spec = importlib.util.spec_from_file_location("curation_batch", path)
    wrapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wrapper)

    class StubRunner:
        def __init__(self, **kwargs):
            self.name = kwargs["default_run_options"]["batch_ligand_name"]
        def execute(self, inputs, *, context):
            payloads = {
                "run_metadata": {"status": "error" if all_failed or self.name == "Failed" else "completed", "error": "fixture failure"},
                "affinity_summary": {"affinity_pred_value": -9. if self.name == "Failed" else -2. if self.name == "Strong" else -1., "affinity_probability_binary": .8},
                "confidence_summary": {"confidence_score": .9},
                "structure_artifacts": {"structure_file": str(tmp_path / (self.name + ".cif"))},
            }
            return {key: make_signal(source="fixture", name=key, value=value, emitted_at=0., spec=None) for key, value in payloads.items()}

    monkeypatch.setattr(wrapper, "Boltz2AffinityPredictor", StubRunner)
    module = wrapper.Boltz2BatchLigandRanker(default_protein_sequence="MKT", default_ligand_csv="name,smiles\nFailed,CC\nWeak,CCC\nStrong,CCCC", work_dir=str(tmp_path))
    outputs = module.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0., run_end=1.))
    batch = unwrap_payload(outputs["batch_summary"])
    metadata = unwrap_payload(outputs["run_metadata"])
    assert batch["status"] == metadata["status"] == ("error" if all_failed else "partial")
    assert batch["submitted_count"] == batch["evaluated_count"] == 3
    assert batch["completed_count"] == (0 if all_failed else 2)
    assert batch["failed_count"] == (3 if all_failed else 1)
    assert len(batch["ranked_ligands"]) == 3
    assert batch["top_ligand_name"] == (None if all_failed else "Strong")
    if not all_failed:
        assert unwrap_payload(outputs["structure_artifacts"])["structure_file"].endswith("Strong.cif")
    interpreter_path = path.parents[2] / "interpreter/src/boltz_prediction_interpreter.py"
    spec = importlib.util.spec_from_file_location("curation_interpreter", interpreter_path)
    interpreter_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(interpreter_module)
    interpreter = interpreter_module.BoltzPredictionInterpreterModel("test", mode="batch", core_alias="core")
    evidence = interpreter.execute({f"core_{key}": value for key, value in outputs.items()}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0., run_end=1.))
    result = unwrap_payload(evidence["prediction_evidence"])
    assert result["top_ligand"] == (None if all_failed else "Strong")
    assert result["ranked_ligand_count"] == (0 if all_failed else 2)
    if not all_failed:
        assert "partial comparison" in result["observed_answer"]
