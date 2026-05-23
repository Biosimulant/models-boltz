# SPDX-FileCopyrightText: 2026-present Biosimulant Team
# SPDX-License-Identifier: Apache-2.0
"""Prediction interpretation stage for source-backed Boltz workflow labs."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from biosim import BioModule
from biosim.signals import AcceptedSignalProfile, BioSignal, SignalSpec, make_signal
from biosim.signals import unwrap_payload as _signal_value


def _record_input_spec(description: str) -> SignalSpec:
    return SignalSpec.record(
        schema={"payload": "json"},
        accepted_profiles=(
            AcceptedSignalProfile(signal_type="record", schema={"payload": "json"}),
            AcceptedSignalProfile(signal_type="scalar"),
        ),
        description=description,
    )


def _coerce_mapping(signal: BioSignal | None) -> dict[str, Any]:
    value = _signal_value(signal)
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}
    return {}


class BoltzPredictionInterpreterModel(BioModule):
    """Convert raw Boltz records into conservative workflow evidence."""

    def __init__(
        self,
        workflow_name: str,
        mode: str = "single",
        core_alias: str = "boltz_boltz2_affinity_predictor",
        caveat: str | None = None,
        integration_step: float = 0.01,
    ) -> None:
        self.integration_step = float(integration_step)
        self.workflow_name = workflow_name
        self.mode = mode
        self.core_alias = core_alias
        self.caveat = caveat or (
            "Boltz-2 outputs are model predictions for hypothesis generation and require independent experimental or computational follow-up."
        )
        self._inputs: dict[str, BioSignal] = {}
        self._outputs: dict[str, BioSignal] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        specs = {
            "scenario_context": _record_input_spec("Curated workflow context."),
            "assembled_boltz_request": _record_input_spec("Resolved Boltz request summary from the input assembler."),
            f"{self.core_alias}_affinity_summary": _record_input_spec("Raw Boltz affinity-style output."),
            f"{self.core_alias}_confidence_summary": _record_input_spec("Raw Boltz confidence output."),
            f"{self.core_alias}_structure_artifacts": _record_input_spec("Raw Boltz structure artifact output."),
            f"{self.core_alias}_run_metadata": _record_input_spec("Raw Boltz run metadata."),
        }
        if self.mode == "batch":
            specs[f"{self.core_alias}_batch_summary"] = _record_input_spec("Raw Boltz batch ranking output.")
        return specs

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            "prediction_evidence": SignalSpec.record(
                schema={"payload": "json"},
                description="Conservative interpretation of Boltz outputs, request provenance, and caveats.",
            )
        }

    def reset(self) -> None:
        self._inputs = {}
        self._outputs = {}

    def set_inputs(self, signals: dict[str, BioSignal]) -> None:
        self._inputs.update(signals or {})

    def advance_window(
        self,
        start: float | None = None,
        end: float | None = None,
        inputs: dict[str, BioSignal] | None = None,
    ) -> dict[str, BioSignal]:
        if inputs:
            self.set_inputs(inputs)
        emitted_at = float(end if end is not None else self.integration_step)
        context = _coerce_mapping(self._inputs.get("scenario_context"))
        assembled = _coerce_mapping(self._inputs.get("assembled_boltz_request"))
        affinity = _coerce_mapping(self._inputs.get(f"{self.core_alias}_affinity_summary"))
        confidence = _coerce_mapping(self._inputs.get(f"{self.core_alias}_confidence_summary"))
        artifacts = _coerce_mapping(self._inputs.get(f"{self.core_alias}_structure_artifacts"))
        run_metadata = _coerce_mapping(self._inputs.get(f"{self.core_alias}_run_metadata"))
        batch = _coerce_mapping(self._inputs.get(f"{self.core_alias}_batch_summary")) if self.mode == "batch" else {}

        metrics = self._metrics(affinity, confidence, batch)
        status = str(run_metadata.get("status") or "not_available")
        evidence = {
            "workflow_name": self.workflow_name,
            "mode": self.mode,
            "target_name": context.get("target_name") or assembled.get("target_name"),
            "ligand_name": context.get("ligand_name") or assembled.get("ligand_name"),
            "run_status": status,
            "scientific_question": context.get("workflow_question"),
            "observed_answer": self._observed_answer(status, metrics),
            "evidence": metrics,
            "artifact_summary": {
                "prediction_dir_name": artifacts.get("prediction_dir_name"),
                "has_structure_file": bool(artifacts.get("structure_file")),
                "has_confidence_file": bool(artifacts.get("confidence_file")),
                "has_affinity_file": bool(artifacts.get("affinity_file")),
            },
            "request_summary": assembled,
            "dominant_module": self.core_alias,
            "caveat": context.get("caveat") or self.caveat,
        }
        if self.mode == "batch":
            evidence["ranked_ligand_count"] = len(batch.get("ranked_ligands") or [])
            evidence["top_ligand"] = batch.get("top_ligand") or metrics.get("top_ligand")

        source = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "prediction_evidence": make_signal(
                source=source,
                name="prediction_evidence",
                value=evidence,
                emitted_at=emitted_at,
                spec=self.outputs()["prediction_evidence"],
            )
        }
        return dict(self._outputs)

    def get_outputs(self) -> dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> list[dict[str, Any]] | None:
        return None

    @staticmethod
    def _metrics(affinity: Mapping[str, Any], confidence: Mapping[str, Any], batch: Mapping[str, Any]) -> dict[str, Any]:
        metrics = {
            "binder_probability": affinity.get("affinity_probability_binary"),
            "affinity_like_value": affinity.get("affinity_pred_value", affinity.get("predicted_affinity_value")),
            "affinity_unit": affinity.get("affinity_unit", affinity.get("predicted_affinity_unit")),
            "confidence_score": confidence.get("confidence_score"),
            "complex_plddt": confidence.get("complex_plddt", confidence.get("mean_plddt")),
            "iptm": confidence.get("iptm", confidence.get("mean_iptm")),
            "ligand_iptm": confidence.get("ligand_iptm", confidence.get("mean_ligand_iptm")),
        }
        ranked = batch.get("ranked_ligands")
        if isinstance(ranked, list) and ranked:
            first = ranked[0] if isinstance(ranked[0], Mapping) else {}
            metrics["top_ligand"] = first.get("ligand") or first.get("name")
            metrics["top_rank_binder_probability"] = first.get("binder_probability")
            metrics["top_rank_affinity_like_value"] = first.get("affinity_like_value")
        return {key: value for key, value in metrics.items() if value not in (None, "")}

    @staticmethod
    def _observed_answer(status: str, metrics: Mapping[str, Any]) -> str:
        if status != "completed":
            return "No completed Boltz prediction is available for this run."
        if "top_ligand" in metrics:
            return f"{metrics['top_ligand']} is the top-ranked completed ligand in this configured batch run."
        if "binder_probability" in metrics:
            return f"Boltz-2 emitted a binder-probability style score of {metrics['binder_probability']} for this configured pair."
        if "confidence_score" in metrics:
            return f"Boltz-2 emitted a structure-confidence score of {metrics['confidence_score']} for this configured pair."
        return "Boltz-2 completed and emitted structure artifacts, but no scalar affinity/confidence metric was available in the parsed summaries."
