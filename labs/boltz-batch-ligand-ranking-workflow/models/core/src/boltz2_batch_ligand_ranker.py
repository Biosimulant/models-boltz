"""Batch ligand ranking wrapper around the single-ligand Boltz-2 runner."""

from __future__ import annotations

import csv
import io
import json
import math
import tempfile
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Dict, Optional

from biosim import BioModule, ExecutionContext, ExecutionPolicy
from biosim.signals import AcceptedSignalProfile, BioSignal, SignalSpec
from biosim.signals import make_signal as _make_signal
from biosim.signals import unwrap_payload as _signal_value

from src.boltz2_affinity_predictor import Boltz2AffinityPredictor


def _generic_input_spec(description: str | None = None) -> SignalSpec:
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


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
    return {key: item for key, item in value.items() if isinstance(key, str)}


class Boltz2BatchLigandRanker(BioModule):
    """Run a small ligand CSV against one protein and rank Boltz-2 outputs."""

    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN

    def __init__(
        self,
        default_protein_sequence: Optional[str] = None,
        default_ligand_csv: Optional[str] = None,
        default_msa_path: Optional[str] = None,
        default_run_options: Optional[Mapping[str, Any]] = None,
        max_ligands: int = 3,
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
        progress_heartbeat_s: float = 30.0,
        integration_step: float = 0.01,
        execution_budget_s: float = 1500.0,
    ) -> None:
        self.integration_step = float(integration_step)
        self._protein_sequence = default_protein_sequence.strip() if isinstance(default_protein_sequence, str) and default_protein_sequence.strip() else None
        self._ligand_csv = default_ligand_csv.strip() if isinstance(default_ligand_csv, str) and default_ligand_csv.strip() else None
        self._msa_path = default_msa_path.strip() if isinstance(default_msa_path, str) and default_msa_path.strip() else None
        self._run_options = _coerce_run_options(default_run_options)
        self.max_ligands = max(1, int(max_ligands))
        if isinstance(execution_budget_s, bool) or not isinstance(execution_budget_s, (int, float)) or not math.isfinite(execution_budget_s) or not 0 < execution_budget_s <= 1500:
            raise ValueError("execution_budget_s must be finite and in (0, 1500]")
        self.execution_budget_s = float(execution_budget_s)

        self.runner_kwargs = {
            "boltz_executable": boltz_executable,
            "runtime_mode": runtime_mode,
            "runtime_dir": runtime_dir,
            "runtime_python": runtime_python,
            "boltz_package_spec": boltz_package_spec,
            "upgrade_runtime": upgrade_runtime,
            "cache_dir": cache_dir,
            "use_msa_server": use_msa_server,
            "accelerator": accelerator,
            "devices": devices,
            "output_format": output_format,
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "override": override,
            "command_timeout_s": command_timeout_s,
            "runtime_setup_timeout_s": runtime_setup_timeout_s,
            "progress_heartbeat_s": progress_heartbeat_s,
            "integration_step": integration_step,
        }
        self.work_dir = Path(work_dir).resolve() if work_dir else None
        self._outputs: Dict[str, BioSignal] = {}
        self._output_payloads: Dict[str, Any] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        return {
            "protein_sequence": _generic_input_spec("Amino-acid sequence for the shared target protein."),
            "ligand_csv": _generic_input_spec("CSV text with name, smiles, and optional metadata columns."),
            "msa_path": _generic_input_spec("Optional path to a precomputed MSA."),
            "run_options": _generic_input_spec("Optional workflow/runtime options."),
        }

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            "batch_summary": SignalSpec.record(schema={"payload": "json"}, description="Ranked ligand rows and batch-level interpretation flags."),
            "affinity_summary": SignalSpec.record(schema={"payload": "json"}, description="Affinity summary for the top-ranked completed ligand."),
            "confidence_summary": SignalSpec.record(schema={"payload": "json"}, description="Confidence summary for the top-ranked completed ligand."),
            "structure_artifacts": SignalSpec.record(schema={"payload": "json"}, description="Structure artifacts for the top-ranked completed ligand."),
            "run_metadata": SignalSpec.record(schema={"payload": "json"}, description="Batch execution metadata and per-ligand run statuses."),
        }

    def reset(self) -> None:
        super().reset()
        self._outputs = {}
        self._output_payloads = {}

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        protein_signal = signals.get("protein_sequence")
        if protein_signal is not None:
            protein_sequence = _coerce_string(_signal_value(protein_signal), "sequence")
            self._protein_sequence = protein_sequence

        csv_signal = signals.get("ligand_csv")
        if csv_signal is not None:
            ligand_csv = _coerce_string(_signal_value(csv_signal), "csv")
            self._ligand_csv = ligand_csv

        msa_signal = signals.get("msa_path")
        if msa_signal is not None:
            msa_path = _coerce_string(_signal_value(msa_signal), "path")
            self._msa_path = msa_path

        run_signal = signals.get("run_options")
        if run_signal is not None:
            run_options = _coerce_run_options(_signal_value(run_signal))
            self._run_options = run_options

    def execute(self, inputs: Mapping[str, BioSignal], *, context: ExecutionContext) -> Mapping[str, BioSignal]:
        self.set_inputs(dict(inputs))
        result = self._execute_at_time(0.0, 0.0)
        return dict(result if result is not None else getattr(self, "_outputs", {}))

    def _execute_at_time(self, start: float, end: float) -> None:
        t = float(end)
        if not self._protein_sequence:
            self._set_error_payload("protein_sequence input is required")
            self._emit_outputs(t)
            return
        if not self._ligand_csv:
            self._set_error_payload("ligand_csv input is required")
            self._emit_outputs(t)
            return

        try:
            ligands = self._parse_ligand_csv(self._ligand_csv)
        except Exception as exc:  # noqa: BLE001
            self._set_error_payload(f"failed to parse ligand_csv: {exc}")
            self._emit_outputs(t)
            return

        if not ligands:
            self._set_error_payload("ligand_csv must contain at least one ligand row")
            self._emit_outputs(t)
            return
        if len(ligands) > self.max_ligands:
            self._set_error_payload(
                f"submitted {len(ligands)} ligands exceeds max_ligands={self.max_ligands}; split the CSV into smaller batches",
                metadata={"submitted_count": len(ligands), "evaluated_count": 0},
            )
            self._emit_outputs(t)
            return
        if self._run_options.get("ranking_mode", "active_affinity") not in {"active_affinity", "binder_probability"}:
            self._set_error_payload("ranking_mode must be active_affinity or binder_probability")
            self._emit_outputs(t)
            return

        run_root = self._create_run_root()
        rows: list[dict[str, Any]] = []
        per_ligand_metadata: list[dict[str, Any]] = []
        top_payload: dict[str, Any] | None = None
        deadline = time.monotonic() + self.execution_budget_s
        evaluated_count = 0

        for index, ligand in enumerate(ligands, start=1):
            name = ligand["name"]
            ligand_run_dir = run_root / f"ligand-{index:02d}"
            ligand_options = dict(self._run_options)
            ligand_options.update(
                {
                    "batch_index": index,
                    "batch_ligand_name": name,
                    "workflow_context": "Guided Boltz-2 batch ligand ranking workflow",
                }
            )
            started = time.monotonic() < deadline
            if not started:
                outputs = {}
            else:
                evaluated_count += 1
                runner = Boltz2AffinityPredictor(
                    default_protein_sequence=self._protein_sequence,
                    default_ligand_smiles=ligand["smiles"],
                    default_msa_path=self._msa_path,
                    default_run_options=ligand_options,
                    work_dir=str(ligand_run_dir),
                    execution_deadline=deadline,
                    **self.runner_kwargs,
                )
                outputs = runner.execute(
                    {},
                    context=ExecutionContext(
                        policy=ExecutionPolicy.ONCE_BEFORE_RUN,
                        run_start=0.0,
                        run_end=1.0,
                    ),
                )
            affinity = self._payload(outputs.get("affinity_summary"))
            confidence = self._payload(outputs.get("confidence_summary"))
            artifacts = self._payload(outputs.get("structure_artifacts"))
            metadata = self._payload(outputs.get("run_metadata"))
            if not started:
                metadata = {"status": "not_started", "error": "batch execution budget exhausted before this ligand could start"}
            status = metadata.get("status") if isinstance(metadata, Mapping) else "error"

            row = self._build_row(index, ligand, status, affinity, confidence, metadata)
            status = row["status"]
            rows.append(row)
            per_ligand_metadata.append(
                {
                    "rank_input_order": index,
                    "ligand": name,
                    "status": status,
                    "error": row["error"],
                }
            )
            if status == "completed":
                candidate = {
                    "ligand": ligand,
                    "row": row,
                    "affinity_summary": dict(affinity) if isinstance(affinity, Mapping) else {},
                    "confidence_summary": dict(confidence) if isinstance(confidence, Mapping) else {},
                    "structure_artifacts": dict(artifacts) if isinstance(artifacts, Mapping) else {},
                    "run_metadata": dict(metadata) if isinstance(metadata, Mapping) else {},
                }
                if top_payload is None or self._ranking_key(row) > self._ranking_key(top_payload["row"]):
                    top_payload = candidate

        ranked_rows = self._rank_rows(rows)
        completed_count = sum(row["status"] == "completed" for row in rows)
        batch_status = "completed" if completed_count == len(ligands) else "partial" if completed_count else "error"
        batch_summary = {
            "status": batch_status,
            "ranking_mode": self._run_options.get("ranking_mode", "active_affinity"),
            "ranking_basis": (
                "affinity_pred_value ascending (lower is stronger); binder probability descending breaks ties"
                if self._run_options.get("ranking_mode", "active_affinity") == "active_affinity"
                else "affinity_probability_binary descending; affinity_pred_value ascending breaks ties"
            ),
            "affinity_unit": "log10(IC50 / micromolar)",
            "top_ligand_name": top_payload["ligand"]["name"] if top_payload else None,
            "ligand_count": len(ligands),
            "submitted_count": len(ligands),
            "evaluated_count": evaluated_count,
            "not_started_count": len(ligands) - evaluated_count,
            "execution_budget_s": self.execution_budget_s,
            "completed_count": completed_count,
            "failed_count": len(rows) - completed_count,
            "ranked_ligands": ranked_rows,
            "interpretation": "Predicted scores only; no experimental validation or uncertainty interval. Active-affinity mode assumes known active ligands for the submitted target.",
        }
        if top_payload is None:
            self._set_error_payload(
                "all ligand runs failed",
                metadata={**batch_summary, "ligand_runs": per_ligand_metadata},
            )
            self._output_payloads["batch_summary"] = {**batch_summary, "error": "all ligand runs failed"}
            self._emit_outputs(t)
            return

        top_name = top_payload["ligand"]["name"]
        top_affinity = dict(top_payload["affinity_summary"])
        top_affinity["top_ligand_name"] = top_name
        top_confidence = dict(top_payload["confidence_summary"])
        top_confidence["top_ligand_name"] = top_name
        top_artifacts = dict(top_payload["structure_artifacts"])
        top_metadata = dict(top_payload["run_metadata"])
        top_metadata.update(
            {
                "status": batch_status,
                "workflow_name": self._run_options.get("workflow_name", "Batch Ligand Ranking"),
                "batch_status": batch_status,
                "ligand_count": len(ligands),
                "submitted_count": len(ligands),
                "evaluated_count": evaluated_count,
                "not_started_count": len(ligands) - evaluated_count,
                "execution_budget_s": self.execution_budget_s,
                "completed_count": sum(1 for row in rows if row["status"] == "completed"),
                "failed_count": sum(1 for row in rows if row["status"] != "completed"),
                "top_ligand_name": top_name,
                "ligand_runs": per_ligand_metadata,
            }
        )
        self._output_payloads = {
            "batch_summary": batch_summary,
            "affinity_summary": top_affinity,
            "confidence_summary": top_confidence,
            "structure_artifacts": top_artifacts,
            "run_metadata": top_metadata,
        }
        self._emit_outputs(t)

    def visualize(self) -> None:
        return None

    def _create_run_root(self) -> Path:
        base_dir = self.work_dir
        if base_dir is not None:
            base_dir.mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix="boltz2-batch-", dir=str(base_dir) if base_dir else None)).resolve()

    def _parse_ligand_csv(self, text: str) -> list[dict[str, str]]:
        reader = csv.DictReader(io.StringIO(text.strip()))
        fieldnames = {name.lower(): name for name in (reader.fieldnames or []) if name}
        name_key = fieldnames.get("name") or fieldnames.get("ligand") or fieldnames.get("ligand_name")
        smiles_key = fieldnames.get("smiles")
        if not name_key or not smiles_key:
            raise ValueError("CSV must include name and smiles columns")
        ligands: list[dict[str, str]] = []
        for index, row in enumerate(reader, start=1):
            name = str(row.get(name_key) or f"Ligand {index}").strip()
            smiles = str(row.get(smiles_key) or "").strip()
            if not smiles:
                raise ValueError(f"row {index} ({name}) has no SMILES; no ligands were evaluated")
            metadata = {
                key: str(value).strip()
                for key, value in row.items()
                if key not in {name_key, smiles_key} and value is not None and str(value).strip()
            }
            ligands.append({"name": name, "smiles": smiles, "metadata": json.dumps(metadata, sort_keys=True)})
        return ligands

    def _build_row(
        self,
        index: int,
        ligand: Mapping[str, str],
        status: Any,
        affinity: Any,
        confidence: Any,
        metadata: Any,
    ) -> dict[str, Any]:
        affinity = affinity if isinstance(affinity, Mapping) else {}
        confidence = confidence if isinstance(confidence, Mapping) else {}
        metadata = metadata if isinstance(metadata, Mapping) else {}
        binder_probability = affinity.get("affinity_probability_binary")
        affinity_value = affinity.get("affinity_pred_value")
        error = metadata.get("error")
        if status == "completed" and (
            not self._finite_number(binder_probability)
            or not 0 <= binder_probability <= 1
            or not self._finite_number(affinity_value)
        ):
            status = "error"
            error = "invalid affinity output: require finite affinity and binder probability in [0, 1]"
        if not self._finite_number(binder_probability) or not 0 <= binder_probability <= 1:
            binder_probability = None
        if not self._finite_number(affinity_value):
            affinity_value = None
        confidence_score = confidence.get("confidence_score")
        if not self._finite_number(confidence_score):
            confidence_score = None
        flags = self._flags(status, binder_probability, confidence_score, metadata)
        return {
            "input_order": index,
            "rank": None,
            "ligand": ligand["name"],
            "smiles": ligand["smiles"],
            "binder_probability": binder_probability,
            "affinity_like_value": affinity_value,
            "confidence": confidence_score,
            "status": status,
            "error": error,
            "flags": flags,
        }

    def _rank_rows(self, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        ranked = sorted(rows, key=self._ranking_key, reverse=True)
        rank = 0
        for row in ranked:
            if row.get("status") == "completed":
                rank += 1
                row["rank"] = rank
            else:
                row["rank"] = None
        return ranked

    @staticmethod
    def _finite_number(value: Any) -> bool:
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)

    def _ranking_key(self, row: Mapping[str, Any]) -> tuple[bool, float, float]:
        probability = row.get("binder_probability")
        affinity = row.get("affinity_like_value")
        if row.get("status") != "completed" or not self._finite_number(probability) or not self._finite_number(affinity):
            return False, float("-inf"), float("-inf")
        if self._run_options.get("ranking_mode", "active_affinity") == "binder_probability":
            return True, float(probability), -float(affinity)
        return True, -float(affinity), float(probability)

    def _flags(self, status: Any, binder_probability: Any, confidence_score: Any, metadata: Mapping[str, Any]) -> list[str]:
        flags: list[str] = []
        if status != "completed":
            flags.append("run failed")
            if metadata.get("error"):
                flags.append("check run metadata")
            return flags
        flags.append("review pose before follow-up; score thresholds are not calibrated")
        return flags

    def _payload(self, signal: BioSignal | None) -> Any:
        if signal is None:
            return {}
        value = _signal_value(signal)
        return value if value is not None else {}

    def _set_error_payload(self, error: str, metadata: Optional[dict[str, Any]] = None) -> None:
        run_metadata = metadata or {}
        run_metadata.update({"status": "error", "error": error})
        self._output_payloads = {
            "batch_summary": {"status": "error", "error": error, "ranked_ligands": []},
            "affinity_summary": {},
            "confidence_summary": {},
            "structure_artifacts": {},
            "run_metadata": run_metadata,
        }

    def _emit_outputs(self, t: float) -> None:
        source = getattr(self, "_world_name", self.__class__.__name__)
        specs = self.outputs()
        self._outputs = {
            name: _make_signal(source=source, name=name, value=self._output_payloads.get(name, {}), emitted_at=t, spec=specs.get(name))
            for name in specs
        }
