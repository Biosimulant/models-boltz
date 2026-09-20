"""Verify elapsed-time limits without running Boltz or contacting services."""
import importlib.util
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from biosim import ExecutionContext, ExecutionPolicy
from biosim.signals import make_signal, unwrap_payload

CORE = Path(__file__).resolve().parents[1]


def load(filename):
    spec = importlib.util.spec_from_file_location("budget_" + filename, CORE / "src" / (filename + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("first_completed", [True, False])
def test_batch_stops_launching_ligands_and_preserves_all_rows(monkeypatch, tmp_path, first_completed):
    wrapper = load("boltz2_batch_ligand_ranker")
    clock = [100.0]
    monkeypatch.setattr(wrapper, "time", SimpleNamespace(monotonic=lambda: clock[0]))
    calls = []

    class Runner:
        def __init__(self, **kwargs):
            assert kwargs["execution_deadline"] == 110.0
            calls.append(kwargs["default_run_options"]["batch_ligand_name"])

        def execute(self, inputs, *, context):
            clock[0] = 110.0
            payloads = {
                "affinity_summary": {"affinity_probability_binary": .8, "affinity_pred_value": -1.0},
                "run_metadata": {"status": "completed" if first_completed else "error", "error": None if first_completed else "timeout"},
            }
            return {k: make_signal(source="test", name=k, value=v, emitted_at=0, spec=None) for k, v in payloads.items()}

    monkeypatch.setattr(wrapper, "Boltz2AffinityPredictor", Runner)
    model = wrapper.Boltz2BatchLigandRanker(default_protein_sequence="ACD", default_ligand_csv="name,smiles\nA,C\nB,CC\nC,CCC", work_dir=str(tmp_path), execution_budget_s=10)
    outputs = model.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0, run_end=1))
    batch = unwrap_payload(outputs["batch_summary"])
    metadata = unwrap_payload(outputs["run_metadata"])
    assert calls == ["A"]
    assert batch["status"] == metadata["status"] == ("partial" if first_completed else "error")
    assert batch["submitted_count"] == 3
    assert batch["evaluated_count"] == 1
    assert batch["not_started_count"] == 2
    assert batch["completed_count"] == int(first_completed)
    assert batch["failed_count"] == 3 - int(first_completed)
    assert batch["top_ligand_name"] == ("A" if first_completed else None)
    rows = {row["ligand"]: row for row in batch["ranked_ligands"]}
    for name in ("B", "C"):
        assert rows[name]["status"] == "not_started"
        assert rows[name]["rank"] is None
        assert rows[name]["affinity_like_value"] is None


def test_setup_prediction_and_retry_share_remaining_time(monkeypatch, tmp_path):
    module = load("boltz2_affinity_predictor")
    clock = [100.0]
    monkeypatch.setattr(module, "time", SimpleNamespace(monotonic=lambda: clock[0]))
    timeouts = []

    def run(command, **kwargs):
        timeouts.append(kwargs["timeout"])
        clock[0] += 4
        return subprocess.CompletedProcess(command, 0, "3.12", "")

    monkeypatch.setattr(module, "subprocess", SimpleNamespace(run=run))
    model = module.Boltz2AffinityPredictor(execution_deadline=110.0)
    assert model._python_version_string(Path(sys.executable)) == "3.12"
    model._run_setup_command(["test-setup"], tmp_path, {"runtime_setup_commands": []})
    model._run_predict_command(["test-predict"], tmp_path)
    assert timeouts == [10, 6, 2]
    with pytest.raises(TimeoutError, match="budget exhausted"):
        model._run_predict_command(["test-retry"], tmp_path)
    assert timeouts == [10, 6, 2]


def test_deadline_terminates_a_real_subprocess(tmp_path):
    module = load("boltz2_affinity_predictor")
    model = module.Boltz2AffinityPredictor(execution_deadline=time.monotonic() + 0.15)
    with pytest.raises(subprocess.TimeoutExpired):
        model._run_predict_command([sys.executable, "-c", "import time; time.sleep(10)"], tmp_path)


@pytest.mark.parametrize("budget", [0, -1, 1501, float("nan"), float("inf"), True, "1500"])
def test_budget_must_fit_managed_limit(budget):
    wrapper = load("boltz2_batch_ligand_ranker")
    with pytest.raises(ValueError, match="execution_budget_s"):
        wrapper.Boltz2BatchLigandRanker(execution_budget_s=budget)
