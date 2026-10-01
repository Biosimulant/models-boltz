"""Evidence retention checks without GPU inference or network access."""
import hashlib
import importlib.util
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1]


def wrapper():
    spec = importlib.util.spec_from_file_location("provenance_affinity", CORE / "src" / "boltz2_affinity_predictor.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scenario(tmp_path, monkeypatch, *, server=True):
    module = wrapper()
    run = tmp_path / "run"
    run.mkdir()
    (run / "request.yaml").write_text("version: 1\n")
    root = run / "output" / "boltz_results_request"
    (root / "processed" / "msa").mkdir(parents=True)
    (root / "processed" / "msa" / "request_0.npz").write_bytes(b"exact processed alignment fixture")
    raw_bytes = b"key,sequence\n-1,ACDE\n-1,AC-E\n"
    if server:
        (root / "msa").mkdir()
        raw = root / "msa" / "request_0.csv"
    else:
        raw = tmp_path / "supplied.csv"
    raw.write_bytes(raw_bytes)
    cache = tmp_path / "weights"
    cache.mkdir()
    for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt"):
        (cache / name).write_bytes((name + " observed weight bytes").encode())
    # Exercise a real child interpreter with a minimal isolated torch fixture.
    # This proves inventory execution and retention, not native GPU inference.
    probe_modules = tmp_path / "probe_modules"
    probe_modules.mkdir()
    (probe_modules / "torch.py").write_text("__version__='fixture'\nclass version: cuda=None\nclass cuda:\n @staticmethod\n def is_available(): return False\n @staticmethod\n def device_count(): return 0\n")
    monkeypatch.setenv("PYTHONPATH", str(probe_modules))
    model = module.Boltz2AffinityPredictor(cache_dir=str(cache), use_msa_server=server, default_msa_path=None if server else str(raw))
    monkeypatch.setattr(model, "_venv_path", lambda root, name: Path(sys.executable))
    return model, run, root, raw


@pytest.mark.parametrize("server", [True, False])
def test_retains_exact_server_or_supplied_alignment_and_runtime(tmp_path, monkeypatch, server):
    model, run, root, raw = scenario(tmp_path, monkeypatch, server=server)
    evidence = model._collect_execution_provenance(run, run / "output", model._resolved_options())
    assert evidence["status"] == "complete" and not evidence["gaps"]
    assert evidence["environment"]["python_executable"] == sys.executable
    assert evidence["environment"]["packages"]
    assert evidence["environment"]["torch_version"] == "fixture"
    assert evidence["environment"]["cuda_available"] is False
    expected_files = {"raw": raw, "processed": root / "processed" / "msa" / "request_0.npz"}
    archive = Path(evidence["alignment_archive_file"])
    assert evidence["alignment_archive"]["sha256"] == hashlib.sha256(archive.read_bytes()).hexdigest()
    with zipfile.ZipFile(archive) as zip_file:
        assert zip_file.testzip() is None
        for item in evidence["alignment_files"]:
            data = zip_file.read(item["archive_member"])
            assert data == expected_files[item["role"]].read_bytes()
            assert len(data) == item["size_bytes"]
            assert hashlib.sha256(data).hexdigest() == item["sha256"]
        assert zip_file.read("request.yaml") == (run / "request.yaml").read_bytes()
        assert json.loads(zip_file.read("alignment-manifest.json")) == evidence["alignment_files"]
    for item in evidence["weights"]:
        assert item["sha256"] == hashlib.sha256((model.cache_dir / item["file_name"]).read_bytes()).hexdigest()


@pytest.mark.parametrize("missing", ["boltz2_conf.ckpt", "boltz2_aff.ckpt", "raw", "processed"])
def test_missing_evidence_is_explicit_and_never_complete(tmp_path, monkeypatch, missing):
    model, run, root, raw = scenario(tmp_path, monkeypatch)
    path = raw if missing == "raw" else root / "processed" / "msa" / "request_0.npz" if missing == "processed" else model.cache_dir / missing
    path.unlink()
    evidence = model._collect_execution_provenance(run, run / "output", model._resolved_options())
    assert evidence["status"] == "incomplete"
    assert evidence["gaps"]
    if missing in {"raw", "processed"}:
        assert "alignment_archive_file" not in evidence


def test_external_executable_is_not_mislabeled_with_host_environment(tmp_path, monkeypatch):
    model, run, root, raw = scenario(tmp_path, monkeypatch)
    resolved = {**model._resolved_options(), "runtime_mode": "external"}
    evidence = model._collect_execution_provenance(run, run / "output", resolved)
    assert evidence["status"] == "incomplete"
    assert "environment" not in evidence
    assert any("interpreter" in gap for gap in evidence["gaps"])


def test_expired_budget_is_reported_as_provenance_gap(tmp_path, monkeypatch):
    model, run, root, raw = scenario(tmp_path, monkeypatch)
    model.execution_deadline = 0
    evidence = model._collect_execution_provenance(run, run / "output", model._resolved_options())
    assert evidence["status"] == "incomplete"
    assert "alignment_archive_file" not in evidence
    assert any("budget" in gap for gap in evidence["gaps"])


def test_batch_preserves_evidence_for_candidates_other_than_the_winner(tmp_path, monkeypatch):
    from biosim import ExecutionContext, ExecutionPolicy
    from biosim.signals import make_signal, unwrap_payload
    spec = importlib.util.spec_from_file_location("provenance_batch", CORE / "src" / "boltz2_batch_ligand_ranker.py")
    batch_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(batch_module)

    class Runner:
        def __init__(self, **kwargs):
            self.index = kwargs['default_run_options']['batch_index']

        def execute(self, inputs, *, context):
            evidence = {'status': 'complete', 'alignment_archive_file': str(tmp_path / f'alignment-{self.index}.zip')}
            payloads = {
                'affinity_summary': {'affinity_pred_value': -self.index, 'affinity_probability_binary': .8},
                'confidence_summary': {'confidence_score': .9},
                'structure_artifacts': {'alignment_archive_file': evidence['alignment_archive_file']},
                'run_metadata': {'status': 'completed', 'execution_provenance': evidence},
            }
            return {key: make_signal(source='test', name=key, value=value, emitted_at=0, spec=None) for key, value in payloads.items()}

    monkeypatch.setattr(batch_module, 'Boltz2AffinityPredictor', Runner)
    model = batch_module.Boltz2BatchLigandRanker(default_protein_sequence='ACDE', default_ligand_csv='name,smiles\nfirst,C\nsecond,CC\nthird,CCC', work_dir=str(tmp_path))
    outputs = model.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0, run_end=1))
    summary = unwrap_payload(outputs['batch_summary'])
    metadata = unwrap_payload(outputs['run_metadata'])
    artifacts = unwrap_payload(outputs['structure_artifacts'])
    assert summary['top_ligand_name'] == 'third'
    assert [item['execution_provenance']['status'] for item in metadata['ligand_runs']] == ['complete'] * 3
    for index in range(1, 4):
        assert artifacts[f'ligand_{index}_alignment_archive_file'] == str(tmp_path / f'alignment-{index}.zip')
