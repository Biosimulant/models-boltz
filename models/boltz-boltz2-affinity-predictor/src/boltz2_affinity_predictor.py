# SPDX-FileCopyrightText: 2025-present Demi <bjaiye1@gmail.com>
#
# SPDX-License-Identifier: MIT
"""Boltz-2 affinity-focused BioModule wrapper with managed local runtime."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Dict, Optional, Set

import yaml

from biosim import BioModule
from biosim.signals import BioSignal, SignalMetadata


def _coerce_string(value: Any, preferred_key: str) -> Optional[str]:
    if isinstance(value, str):
        text = value.strip()
        return text or None
    if isinstance(value, Mapping):
        candidate = value.get(preferred_key)
        if isinstance(candidate, str):
            text = candidate.strip()
            return text or None
    return None


def _coerce_run_options(value: Any) -> Dict[str, Any]:
    if not isinstance(value, Mapping):
        return {}
    out: Dict[str, Any] = {}
    for key, item in value.items():
        if isinstance(key, str):
            out[key] = item
    return out


class Boltz2AffinityPredictor(BioModule):
    """Run `boltz predict` once and surface compact structured outputs."""

    def __init__(
        self,
        boltz_executable: str = "boltz",
        runtime_mode: str = "managed",
        runtime_dir: Optional[str] = None,
        runtime_python: Optional[str] = None,
        boltz_package_spec: Optional[str] = None,
        upgrade_runtime: bool = False,
        work_dir: Optional[str] = None,
        cache_dir: Optional[str] = None,
        use_msa_server: bool = False,
        accelerator: str = "gpu",
        devices: int = 1,
        output_format: str = "mmcif",
        recycling_steps: int = 3,
        sampling_steps: int = 200,
        diffusion_samples: int = 1,
        override: bool = True,
        command_timeout_s: float = 3600.0,
        runtime_setup_timeout_s: float = 1800.0,
        min_dt: float = 0.01,
    ) -> None:
        self.min_dt = min_dt
        self.boltz_executable = boltz_executable
        self.runtime_mode = runtime_mode
        self.runtime_python = runtime_python
        self.boltz_package_spec = boltz_package_spec
        self.upgrade_runtime = upgrade_runtime
        self.work_dir = Path(work_dir).resolve() if work_dir else None
        self.cache_dir = Path(cache_dir).resolve() if cache_dir else None
        self.use_msa_server = use_msa_server
        self.accelerator = accelerator
        self.devices = devices
        self.output_format = output_format
        self.recycling_steps = recycling_steps
        self.sampling_steps = sampling_steps
        self.diffusion_samples = diffusion_samples
        self.override = override
        self.command_timeout_s = command_timeout_s
        self.runtime_setup_timeout_s = runtime_setup_timeout_s

        repo_root = Path(__file__).resolve().parents[3]
        self.runtime_dir = (
            Path(runtime_dir).expanduser().resolve()
            if runtime_dir
            else (repo_root / ".runtime" / "boltz2").resolve()
        )

        self._protein_sequence: Optional[str] = None
        self._ligand_smiles: Optional[str] = None
        self._msa_path: Optional[str] = None
        self._run_options: Dict[str, Any] = {}
        self._outputs: Dict[str, BioSignal] = {}
        self._cached_payloads: Dict[str, Any] = {}
        self._last_signature: Optional[str] = None

    def inputs(self) -> Set[str]:
        return {"protein_sequence", "ligand_smiles", "msa_path", "run_options"}

    def outputs(self) -> Set[str]:
        return {"affinity_summary", "confidence_summary", "structure_artifacts", "run_metadata"}

    def reset(self) -> None:
        self._outputs = {}
        self._cached_payloads = {}
        self._last_signature = None

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        changed = False

        protein_signal = signals.get("protein_sequence")
        if protein_signal is not None:
            protein_sequence = _coerce_string(protein_signal.value, "sequence")
            if protein_sequence != self._protein_sequence:
                self._protein_sequence = protein_sequence
                changed = True

        ligand_signal = signals.get("ligand_smiles")
        if ligand_signal is not None:
            ligand_smiles = _coerce_string(ligand_signal.value, "smiles")
            if ligand_smiles != self._ligand_smiles:
                self._ligand_smiles = ligand_smiles
                changed = True

        msa_signal = signals.get("msa_path")
        if msa_signal is not None:
            msa_path = _coerce_string(msa_signal.value, "path")
            if msa_path != self._msa_path:
                self._msa_path = msa_path
                changed = True

        run_signal = signals.get("run_options")
        if run_signal is not None:
            run_options = _coerce_run_options(run_signal.value)
            if run_options != self._run_options:
                self._run_options = run_options
                changed = True

        if changed:
            self._last_signature = None

    def advance_to(self, t: float) -> None:
        resolved = self._resolved_options()
        signature = json.dumps(
            {
                "protein_sequence": self._protein_sequence,
                "ligand_smiles": self._ligand_smiles,
                "msa_path": self._msa_path,
                "run_options": resolved,
            },
            sort_keys=True,
            default=str,
        )

        if signature == self._last_signature and self._cached_payloads:
            self._emit_outputs(t)
            return

        if not self._protein_sequence:
            self._set_error_payload("protein_sequence input is required")
            self._last_signature = signature
            self._emit_outputs(t)
            return
        if not self._ligand_smiles:
            self._set_error_payload("ligand_smiles input is required")
            self._last_signature = signature
            self._emit_outputs(t)
            return
        if not resolved["use_msa_server"] and not self._msa_path:
            self._set_error_payload("msa_path is required unless use_msa_server is enabled")
            self._last_signature = signature
            self._emit_outputs(t)
            return

        run_root = self._create_run_root()
        request_path = run_root / "request.yaml"
        request_path.write_text(self._build_request_document(resolved), encoding="utf-8")
        output_dir = run_root / "output"

        metadata: Dict[str, Any] = {
            "status": "running",
            "command": [],
            "cwd": str(run_root),
            "input_yaml_path": str(request_path),
            "output_dir": str(output_dir),
            "use_msa_server": bool(resolved["use_msa_server"]),
            "stdout": "",
            "stderr": "",
            "runtime_mode": resolved["runtime_mode"],
            "runtime_dir": str(self.runtime_dir),
            "runtime_bootstrapped": False,
            "runtime_setup_commands": [],
        }

        try:
            resolved_executable = self._resolve_boltz_executable(run_root, resolved, metadata)
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = f"failed to prepare Boltz runtime: {exc}"
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        command = self._build_command(resolved_executable, request_path, output_dir, resolved)
        metadata["command"] = command

        try:
            completed = subprocess.run(
                command,
                cwd=run_root,
                capture_output=True,
                text=True,
                timeout=self.command_timeout_s,
                check=False,
            )
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = f"failed to execute boltz: {exc}"
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        metadata["stdout"] = completed.stdout
        metadata["stderr"] = completed.stderr
        metadata["returncode"] = completed.returncode

        if completed.returncode != 0:
            metadata["status"] = "error"
            metadata["error"] = "boltz predict returned a non-zero exit code"
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        try:
            prediction_dir = self._find_prediction_dir(output_dir)
            confidence_path = self._require_single(prediction_dir.glob("confidence_*_model_0.json"))
            affinity_path = self._require_single(prediction_dir.glob("affinity_*.json"))
            structure_path = self._find_structure_file(prediction_dir)
        except Exception as exc:  # noqa: BLE001
            metadata["status"] = "error"
            metadata["error"] = f"expected Boltz outputs were not found: {exc}"
            self._set_error_payload(metadata["error"], metadata=metadata)
            self._last_signature = signature
            self._emit_outputs(t)
            return

        confidence_summary = self._load_json(confidence_path)
        affinity_summary = self._load_json(affinity_path)
        artifacts = {
            "prediction_dir": str(prediction_dir.resolve()),
            "structure_file": str(structure_path.resolve()),
            "confidence_file": str(confidence_path.resolve()),
            "affinity_file": str(affinity_path.resolve()),
        }

        optional_files = {
            "pae_file": next(prediction_dir.glob("pae_*_model_0.npz"), None),
            "pde_file": next(prediction_dir.glob("pde_*_model_0.npz"), None),
            "plddt_file": next(prediction_dir.glob("plddt_*_model_0.npz"), None),
        }
        for key, maybe_path in optional_files.items():
            if maybe_path is not None:
                artifacts[key] = str(Path(maybe_path).resolve())

        metadata["status"] = "completed"
        metadata["prediction_dir"] = str(prediction_dir.resolve())
        self._cached_payloads = {
            "affinity_summary": affinity_summary,
            "confidence_summary": confidence_summary,
            "structure_artifacts": artifacts,
            "run_metadata": metadata,
        }
        self._last_signature = signature
        self._emit_outputs(t)

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def _resolved_options(self) -> Dict[str, Any]:
        resolved: Dict[str, Any] = {
            "use_msa_server": self.use_msa_server,
            "accelerator": self.accelerator,
            "devices": self.devices,
            "output_format": self.output_format,
            "recycling_steps": self.recycling_steps,
            "sampling_steps": self.sampling_steps,
            "diffusion_samples": self.diffusion_samples,
            "override": self.override,
            "template_path": None,
            "msa_server_url": None,
            "max_parallel_samples": None,
            "affinity_mw_correction": False,
            "runtime_mode": self.runtime_mode,
            "upgrade_runtime": self.upgrade_runtime,
            "boltz_package_spec": self.boltz_package_spec,
        }
        for key, value in self._run_options.items():
            if key in resolved:
                resolved[key] = value

        if not isinstance(resolved["boltz_package_spec"], str) or not resolved["boltz_package_spec"].strip():
            resolved["boltz_package_spec"] = "boltz[cuda]" if resolved["accelerator"] == "gpu" else "boltz"
        return resolved

    def _create_run_root(self) -> Path:
        base_dir = self.work_dir
        if base_dir is not None:
            base_dir.mkdir(parents=True, exist_ok=True)
        root = Path(tempfile.mkdtemp(prefix="boltz2-run-", dir=str(base_dir) if base_dir else None))
        return root.resolve()

    def _build_request_document(self, resolved: Mapping[str, Any]) -> str:
        request: Dict[str, Any] = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": self._protein_sequence,
                    }
                },
                {
                    "ligand": {
                        "id": "B",
                        "smiles": self._ligand_smiles,
                    }
                },
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "B",
                    }
                }
            ],
        }
        if self._msa_path and not resolved["use_msa_server"]:
            request["sequences"][0]["protein"]["msa"] = str(Path(self._msa_path).expanduser().resolve())

        template_path = resolved.get("template_path")
        if isinstance(template_path, str) and template_path.strip():
            template = Path(template_path).expanduser().resolve()
            key = "pdb" if template.suffix.lower() == ".pdb" else "cif"
            request["templates"] = [{key: str(template)}]

        return yaml.safe_dump(request, sort_keys=False)

    def _build_command(
        self,
        boltz_executable: str,
        request_path: Path,
        output_dir: Path,
        resolved: Mapping[str, Any],
    ) -> list[str]:
        command = [
            boltz_executable,
            "predict",
            str(request_path),
            "--out_dir",
            str(output_dir),
            "--accelerator",
            str(resolved["accelerator"]),
            "--devices",
            str(resolved["devices"]),
            "--output_format",
            str(resolved["output_format"]),
            "--recycling_steps",
            str(resolved["recycling_steps"]),
            "--sampling_steps",
            str(resolved["sampling_steps"]),
            "--diffusion_samples",
            str(resolved["diffusion_samples"]),
        ]
        if self.cache_dir is not None:
            command.extend(["--cache", str(self.cache_dir)])
        if resolved["use_msa_server"]:
            command.append("--use_msa_server")
        if isinstance(resolved.get("msa_server_url"), str) and resolved["msa_server_url"]:
            command.extend(["--msa_server_url", str(resolved["msa_server_url"])])
        if resolved["override"]:
            command.append("--override")
        if resolved.get("max_parallel_samples") is not None:
            command.extend(["--max_parallel_samples", str(resolved["max_parallel_samples"])])
        if bool(resolved.get("affinity_mw_correction")):
            command.append("--affinity_mw_correction")
        return command

    def _resolve_boltz_executable(
        self,
        run_root: Path,
        resolved: Mapping[str, Any],
        metadata: Dict[str, Any],
    ) -> str:
        runtime_mode = str(resolved["runtime_mode"]).strip().lower()
        if runtime_mode == "external":
            resolved_exec = shutil.which(self.boltz_executable) if not os.path.isabs(self.boltz_executable) else self.boltz_executable
            if not resolved_exec:
                raise FileNotFoundError(f"could not find boltz executable: {self.boltz_executable}")
            metadata["resolved_boltz_executable"] = str(resolved_exec)
            return str(resolved_exec)
        if runtime_mode != "managed":
            raise ValueError(f"unsupported runtime_mode: {resolved['runtime_mode']}")
        return self._ensure_managed_runtime(run_root, resolved, metadata)

    def _ensure_managed_runtime(
        self,
        run_root: Path,
        resolved: Mapping[str, Any],
        metadata: Dict[str, Any],
    ) -> str:
        runtime_root = self.runtime_dir
        runtime_root.parent.mkdir(parents=True, exist_ok=True)

        venv_python = self._venv_path(runtime_root, "python")
        boltz_bin = self._venv_path(runtime_root, "boltz")
        base_python = self.runtime_python or sys.executable or shutil.which("python3")
        if not base_python:
            raise RuntimeError("no Python executable available to create managed Boltz runtime")

        if not venv_python.exists():
            self._run_setup_command(
                [base_python, "-m", "venv", str(runtime_root)],
                run_root,
                metadata,
            )
            metadata["runtime_bootstrapped"] = True

        if bool(resolved["upgrade_runtime"]) or not boltz_bin.exists():
            self._run_setup_command(
                [str(venv_python), "-m", "pip", "install", "--upgrade", "pip"],
                run_root,
                metadata,
            )
            self._run_setup_command(
                [str(venv_python), "-m", "pip", "install", str(resolved["boltz_package_spec"])],
                run_root,
                metadata,
            )
            metadata["runtime_bootstrapped"] = True

        if not boltz_bin.exists():
            raise FileNotFoundError(f"managed Boltz executable was not created at {boltz_bin}")

        metadata["runtime_python_executable"] = str(venv_python)
        metadata["resolved_boltz_executable"] = str(boltz_bin)
        return str(boltz_bin)

    def _run_setup_command(
        self,
        command: list[str],
        run_root: Path,
        metadata: Dict[str, Any],
    ) -> None:
        metadata["runtime_setup_commands"].append(command)
        completed = subprocess.run(
            command,
            cwd=run_root,
            capture_output=True,
            text=True,
            timeout=self.runtime_setup_timeout_s,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                "runtime setup command failed: "
                f"{command} :: stdout={completed.stdout!r} stderr={completed.stderr!r}"
            )

    def _venv_path(self, runtime_root: Path, executable: str) -> Path:
        bin_dir = "Scripts" if os.name == "nt" else "bin"
        suffix = ".exe" if os.name == "nt" and not executable.endswith(".exe") else ""
        return runtime_root / bin_dir / f"{executable}{suffix}"

    def _find_prediction_dir(self, output_dir: Path) -> Path:
        predictions_root = output_dir / "predictions"
        candidates = sorted(p for p in predictions_root.iterdir() if p.is_dir())
        if not candidates:
            raise FileNotFoundError(f"no prediction folders under {predictions_root}")
        return candidates[0]

    def _find_structure_file(self, prediction_dir: Path) -> Path:
        for pattern in ("*_model_0.cif", "*_model_0.pdb"):
            candidate = next(prediction_dir.glob(pattern), None)
            if candidate is not None:
                return Path(candidate)
        raise FileNotFoundError("no top-ranked structure file was produced")

    def _require_single(self, paths: Any) -> Path:
        items = sorted(Path(p) for p in paths)
        if not items:
            raise FileNotFoundError("required file is missing")
        return items[0]

    def _load_json(self, path: Path) -> Dict[str, Any]:
        loaded = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(loaded, dict):
            return loaded
        raise ValueError(f"expected JSON object in {path}")

    def _set_error_payload(self, error: str, metadata: Optional[Dict[str, Any]] = None) -> None:
        payload = metadata or {
            "status": "error",
            "command": None,
            "cwd": None,
            "input_yaml_path": None,
            "output_dir": None,
            "stdout": "",
            "stderr": "",
            "runtime_mode": self.runtime_mode,
            "runtime_dir": str(self.runtime_dir),
            "runtime_bootstrapped": False,
            "runtime_setup_commands": [],
        }
        payload["status"] = "error"
        payload["error"] = error
        self._cached_payloads = {
            "affinity_summary": {},
            "confidence_summary": {},
            "structure_artifacts": {},
            "run_metadata": payload,
        }

    def _emit_outputs(self, t: float) -> None:
        source = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "affinity_summary": BioSignal(
                source=source,
                name="affinity_summary",
                value=self._cached_payloads.get("affinity_summary", {}),
                time=t,
                metadata=SignalMetadata(description="Parsed Boltz affinity summary for the latest run", kind="metric"),
            ),
            "confidence_summary": BioSignal(
                source=source,
                name="confidence_summary",
                value=self._cached_payloads.get("confidence_summary", {}),
                time=t,
                metadata=SignalMetadata(description="Parsed Boltz confidence summary for the top-ranked prediction", kind="metric"),
            ),
            "structure_artifacts": BioSignal(
                source=source,
                name="structure_artifacts",
                value=self._cached_payloads.get("structure_artifacts", {}),
                time=t,
                metadata=SignalMetadata(description="Absolute paths to the latest Boltz output artifacts", kind="state"),
            ),
            "run_metadata": BioSignal(
                source=source,
                name="run_metadata",
                value=self._cached_payloads.get("run_metadata", {}),
                time=t,
                metadata=SignalMetadata(description="Execution metadata and captured logs for the latest Boltz invocation", kind="metric"),
            ),
        }
