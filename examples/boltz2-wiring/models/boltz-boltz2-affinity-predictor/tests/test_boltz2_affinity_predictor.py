from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
import yaml


def _set_required_inputs(module, BioSignal, *, msa_path: str | None = None, run_options: dict | None = None):
    signals = {
        "protein_sequence": BioSignal(source="test", name="protein_sequence", value="MKTAYIAKQRQISFVKSHFSRQ", time=0.0),
        "ligand_smiles": BioSignal(source="test", name="ligand_smiles", value="CCO", time=0.0),
    }
    if msa_path is not None:
        signals["msa_path"] = BioSignal(source="test", name="msa_path", value=msa_path, time=0.0)
    if run_options is not None:
        signals["run_options"] = BioSignal(source="test", name="run_options", value=run_options, time=0.0)
    module.set_inputs(signals)


def test_instantiation(biosim, tmp_path):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path))
    assert module.min_dt > 0
    assert module.runtime_mode == "managed"
    assert module.cache_dir is not None
    assert module.cache_dir.name == "boltz-cache"
    assert module.inputs() == {"protein_sequence", "ligand_smiles", "msa_path", "run_options"}
    assert module.outputs() == {"affinity_summary", "confidence_summary", "structure_artifacts", "run_metadata"}


def test_request_document_with_explicit_msa(biosim, tmp_path):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    msa = tmp_path / "example.a3m"
    msa.write_text(">query\nMKTAYIAKQRQISFVKSHFSRQ\n", encoding="utf-8")

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path))
    _set_required_inputs(module, BioSignal, msa_path=str(msa))

    document = module._build_request_document(module._resolved_options())
    request = yaml.safe_load(document)
    assert request["sequences"][0]["protein"]["msa"] == str(msa.resolve())
    assert request["properties"][0]["affinity"]["binder"] == "B"


def test_request_document_with_server_mode(biosim, tmp_path):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), use_msa_server=True)
    _set_required_inputs(module, BioSignal)

    document = module._build_request_document(module._resolved_options())
    request = yaml.safe_load(document)
    assert "msa" not in request["sequences"][0]["protein"]


def test_request_document_supports_empty_msa(biosim, tmp_path):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path))
    _set_required_inputs(module, BioSignal, msa_path="empty")

    document = module._build_request_document(module._resolved_options())
    request = yaml.safe_load(document)
    assert request["sequences"][0]["protein"]["msa"] == "empty"


def test_missing_inputs_surface_error_metadata(biosim, tmp_path):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), use_msa_server=True)
    module.set_inputs({
        "ligand_smiles": BioSignal(source="test", name="ligand_smiles", value="CCO", time=0.0),
    })
    module.advance_to(0.1)

    outputs = module.get_outputs()
    assert outputs["run_metadata"].value["status"] == "error"
    assert "protein_sequence" in outputs["run_metadata"].value["error"]
    assert module.visualize() is None


def test_managed_runtime_bootstraps_and_parses_outputs(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    commands: list[list[str]] = []

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        commands.append(command)
        run_root = Path(cwd)

        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")

        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        prediction_dir = run_root / "output" / "predictions" / "request"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "request_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_request_model_0.json").write_text(
            json.dumps({"confidence_score": 0.91, "ptm": 0.82, "iptm": 0.78, "complex_plddt": 0.88}),
            encoding="utf-8",
        )
        (prediction_dir / "affinity_request.json").write_text(
            json.dumps({"affinity_pred_value": -1.2, "affinity_probability_binary": 0.97}),
            encoding="utf-8",
        )
        (prediction_dir / "pae_request_model_0.npz").write_text("placeholder", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    runtime_dir = tmp_path / "managed-runtime"
    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(runtime_dir), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    outputs = module.get_outputs()
    metadata = outputs["run_metadata"].value
    assert metadata["status"] == "completed"
    assert metadata["runtime_bootstrapped"] is True
    assert metadata["cache_dir"].endswith("boltz-cache")
    assert outputs["affinity_summary"].value["affinity_probability_binary"] == 0.97
    assert outputs["confidence_summary"].value["confidence_score"] == 0.91
    assert Path(outputs["structure_artifacts"].value["structure_file"]).is_absolute()
    assert Path(outputs["structure_artifacts"].value["affinity_file"]).is_absolute()
    assert any("-m" in command and "venv" in command for command in commands)
    assert any("-m" in command and "pip" in command and any(item.startswith("boltz") for item in command) for command in commands)
    assert any("--cache" in command for command in commands if command and command[0].endswith("boltz"))
    assert metadata["resolved_boltz_executable"].endswith("/bin/boltz")

    visuals = module.visualize()
    assert visuals is not None
    assert visuals[0]["render"] == "structure3d"
    assert visuals[0]["data"]["format"] == "mmcif"
    assert visuals[0]["data"]["source"]["kind"] == "artifact"
    assert Path(visuals[0]["data"]["source"]["path"]).name == "request_model_0.cif"
    assert visuals[1]["render"] == "table"
    assert visuals[1]["data"]["columns"] == ["Metric", "Value"]


def test_runtime_bootstrap_failure_surfaces_metadata(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            return subprocess.CompletedProcess(command, 2, stdout="", stderr="install failed")
        raise AssertionError("predict command should not run when bootstrap fails")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(tmp_path / "managed-runtime"), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    metadata = module.get_outputs()["run_metadata"].value
    assert metadata["status"] == "error"
    assert "prepare Boltz runtime" in metadata["error"]


def test_subprocess_failure_surfaces_metadata(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")
        return subprocess.CompletedProcess(command, 3, stdout="", stderr="boom")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(tmp_path / "managed-runtime"), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    metadata = module.get_outputs()["run_metadata"].value
    assert metadata["status"] == "error"
    assert metadata["returncode"] == 3
    assert "non-zero" in metadata["error"]


def test_corrupted_cache_is_purged_and_retried_once(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    predict_calls = {"count": 0}

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        predict_calls["count"] += 1
        cache_dir = Path(command[command.index("--cache") + 1])
        if predict_calls["count"] == 1:
            (cache_dir / "mols").mkdir(parents=True, exist_ok=True)
            (cache_dir / "mols.tar").write_text("corrupt", encoding="utf-8")
            (cache_dir / "boltz2_conf.ckpt").write_text("partial", encoding="utf-8")
            return subprocess.CompletedProcess(command, 1, stdout="Extracting the CCD data", stderr="tarfile.ReadError: unexpected end of data")

        run_root = Path(cwd)
        prediction_dir = run_root / "output" / "predictions" / "request"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "request_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_request_model_0.json").write_text(json.dumps({"confidence_score": 0.72}), encoding="utf-8")
        (prediction_dir / "affinity_request.json").write_text(json.dumps({"affinity_probability_binary": 0.61}), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    cache_dir = tmp_path / "local-cache"
    module = Boltz2AffinityPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "managed-runtime"),
        cache_dir=str(cache_dir),
        use_msa_server=True,
    )
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    metadata = module.get_outputs()["run_metadata"].value
    assert metadata["status"] == "completed"
    assert metadata["cache_repaired"] is True
    assert metadata["retry_count"] == 1
    assert predict_calls["count"] == 2
    assert not (cache_dir / "mols.tar").exists()
    assert not (cache_dir / "boltz2_conf.ckpt").exists()


def test_missing_expected_files_becomes_error(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")
        run_root = Path(cwd)
        (run_root / "output" / "predictions" / "request").mkdir(parents=True, exist_ok=True)
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(tmp_path / "managed-runtime"), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    metadata = module.get_outputs()["run_metadata"].value
    assert metadata["status"] == "error"
    assert "expected Boltz outputs" in metadata["error"]


def test_legacy_success_layout_without_affinity_json_is_accepted(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        run_root = Path(cwd)
        output_dir = run_root / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "affinity_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (output_dir / "confidence_affinity_model_0.json").write_text(
            json.dumps({"confidence_score": 0.84, "ptm": 0.8}),
            encoding="utf-8",
        )
        (output_dir / "plddt_affinity_model_0.npz").write_text("placeholder", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "managed-runtime"),
        use_msa_server=True,
    )
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    outputs = module.get_outputs()
    metadata = outputs["run_metadata"].value
    assert metadata["status"] == "completed"
    assert outputs["affinity_summary"].value == {}
    assert outputs["confidence_summary"].value["confidence_score"] == 0.84
    assert Path(outputs["structure_artifacts"].value["structure_file"]).is_absolute()
    assert Path(outputs["structure_artifacts"].value["confidence_file"]).is_absolute()
    assert "affinity_file" not in outputs["structure_artifacts"].value
    assert Path(outputs["structure_artifacts"].value["plddt_file"]).is_absolute()


def test_recursive_output_layout_is_accepted(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        run_root = Path(cwd)
        prediction_dir = run_root / "output" / "request" / "predictions"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "affinity_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_affinity_model_0.json").write_text(
            json.dumps({"confidence_score": 0.92}),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "managed-runtime"),
        use_msa_server=True,
    )
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.5)

    outputs = module.get_outputs()
    assert outputs["run_metadata"].value["status"] == "completed"
    assert outputs["confidence_summary"].value["confidence_score"] == 0.92
    assert Path(outputs["structure_artifacts"].value["structure_file"]).name == "affinity_model_0.cif"


def test_repeat_advance_does_not_rerun_until_reset(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    calls = {"predict": 0}

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        calls["predict"] += 1
        run_root = Path(cwd)
        prediction_dir = run_root / "output" / "predictions" / "request"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "request_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_request_model_0.json").write_text(
            json.dumps({"confidence_score": 0.8}),
            encoding="utf-8",
        )
        (prediction_dir / "affinity_request.json").write_text(
            json.dumps({"affinity_pred_value": 0.1}),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(tmp_path / "managed-runtime"), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.2)
    module.advance_to(0.3)
    assert calls["predict"] == 1

    module.reset()
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.4)
    assert calls["predict"] == 2


def test_constructor_defaults_allow_space_style_usage(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")

        run_root = Path(cwd)
        prediction_dir = run_root / "output" / "predictions" / "request"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "request_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_request_model_0.json").write_text(json.dumps({"confidence_score": 0.77}), encoding="utf-8")
        (prediction_dir / "affinity_request.json").write_text(json.dumps({"affinity_probability_binary": 0.66}), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(
        work_dir=str(tmp_path),
        runtime_dir=str(tmp_path / "managed-runtime"),
        use_msa_server=True,
        default_protein_sequence="MKTAYIAKQRQISFVKSHFSRQ",
        default_ligand_smiles="CCO",
    )
    module.advance_to(0.1)
    outputs = module.get_outputs()
    assert outputs["run_metadata"].value["status"] == "completed"
    assert outputs["affinity_summary"].value["affinity_probability_binary"] == 0.66


def test_managed_runtime_selects_supported_python_when_host_python_is_unsupported(biosim, tmp_path, monkeypatch):
    from src.boltz2_affinity_predictor import Boltz2AffinityPredictor
    from biosim.signals import BioSignal

    real_which = shutil.which

    def fake_which(name):
        if name == "python3.12":
            return "/opt/homebrew/bin/python3.12"
        return real_which(name)

    def fake_run(command, cwd, capture_output, text, timeout, check):  # noqa: ARG001
        command = [str(item) for item in command]
        if "-c" in command:
            if command[0] == sys.executable:
                return subprocess.CompletedProcess(command, 0, stdout="3.14\n", stderr="")
            if command[0].endswith("python3.12") or command[0].endswith("/bin/python"):
                return subprocess.CompletedProcess(command, 0, stdout="3.12\n", stderr="")
        if "-m" in command and "venv" in command:
            runtime_root = Path(command[-1])
            (runtime_root / "bin").mkdir(parents=True, exist_ok=True)
            (runtime_root / "bin" / "python").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="venv", stderr="")
        if "-m" in command and "pip" in command and "install" in command:
            if any(item.startswith("boltz") for item in command):
                runtime_python = Path(command[0])
                runtime_root = runtime_python.parent.parent
                (runtime_root / "bin" / "boltz").write_text("", encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="pip", stderr="")
        run_root = Path(cwd)
        prediction_dir = run_root / "output" / "predictions" / "request"
        prediction_dir.mkdir(parents=True, exist_ok=True)
        (prediction_dir / "request_model_0.cif").write_text("data_mock\n", encoding="utf-8")
        (prediction_dir / "confidence_request_model_0.json").write_text(json.dumps({"confidence_score": 0.77}), encoding="utf-8")
        (prediction_dir / "affinity_request.json").write_text(json.dumps({"affinity_probability_binary": 0.66}), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, stdout="ok", stderr="")

    monkeypatch.setattr(shutil, "which", fake_which)
    monkeypatch.setattr(subprocess, "run", fake_run)

    module = Boltz2AffinityPredictor(work_dir=str(tmp_path), runtime_dir=str(tmp_path / "managed-runtime"), use_msa_server=True)
    _set_required_inputs(module, BioSignal)
    module.advance_to(0.1)

    metadata = module.get_outputs()["run_metadata"].value
    assert metadata["status"] == "completed"
    assert metadata["runtime_base_python"].endswith("python3.12")


def test_example_files_parse_and_reference_real_interface(biosim):
    repo_root = Path(__file__).resolve().parents[3]
    minimal = yaml.safe_load((repo_root / "examples" / "boltz2-minimal" / "config.yaml").read_text(encoding="utf-8"))
    explicit = yaml.safe_load((repo_root / "examples" / "boltz2-explicit-msa" / "config.yaml").read_text(encoding="utf-8"))
    short_no_msa = yaml.safe_load((repo_root / "examples" / "boltz2-short-no-msa" / "config.yaml").read_text(encoding="utf-8"))
    wiring = yaml.safe_load((repo_root / "examples" / "boltz2-wiring" / "space.yaml").read_text(encoding="utf-8"))

    assert minimal["model"]["parameters"]["runtime_mode"] == "managed"
    assert minimal["model"]["parameters"]["use_msa_server"] is True
    assert explicit["model"]["inputs"]["msa_path"] == "./assets/seq1.a3m"
    assert short_no_msa["model"]["inputs"]["msa_path"] == "empty"
    assert short_no_msa["model"]["parameters"]["sampling_steps"] == 1
    assert minimal["model"]["path"] == "../../models/boltz-boltz2-affinity-predictor"
    assert wiring["models"][0]["path"] == "../../models/boltz-boltz2-affinity-predictor"
    assert wiring["models"][0]["parameters"]["runtime_mode"] == "managed"
    assert "default_protein_sequence" in wiring["models"][0]["parameters"]


@pytest.mark.skipif(
    os.getenv("BIOSIM_BOLTZ_RUN_REAL_SMOKE") != "1",
    reason="Set BIOSIM_BOLTZ_RUN_REAL_SMOKE=1 to run the real Boltz smoke test.",
)
def test_real_smoke_example_runs(tmp_path):
    repo_root = Path(__file__).resolve().parents[3]
    output_json = tmp_path / "real-smoke-output.json"
    completed = subprocess.run(
        [
            sys.executable,
            str(repo_root / "examples" / "run_example.py"),
            "boltz2-short-no-msa",
            "--work-dir",
            str(tmp_path / "runs"),
            "--runtime-dir",
            str(tmp_path / "runtime"),
            "--output-json",
            str(output_json),
        ],
        cwd=repo_root,
        capture_output=True,
        text=True,
        timeout=7200,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    outputs = payload["outputs"]
    assert outputs["run_metadata"]["value"]["status"] == "completed"
    structure_file = Path(outputs["structure_artifacts"]["value"]["structure_file"])
    confidence_file = Path(outputs["structure_artifacts"]["value"]["confidence_file"])
    assert structure_file.exists()
    assert confidence_file.exists()
    affinity_file = outputs["structure_artifacts"]["value"].get("affinity_file")
    if affinity_file is not None:
        assert Path(affinity_file).exists()
