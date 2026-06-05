# SPDX-FileCopyrightText: 2026-present Biosimulant Team
# SPDX-License-Identifier: Apache-2.0
"""Input assembly stage for source-backed Boltz workflow labs."""

from __future__ import annotations

import csv
import hashlib
import io
from collections.abc import Mapping
from pathlib import Path
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


def _coerce_string(value: Any, preferred_key: str) -> str | None:
    value = _signal_value(value)
    if isinstance(value, str):
        text = value.strip()
        return text or None
    if isinstance(value, Mapping):
        candidate = value.get(preferred_key) or value.get("payload")
        if isinstance(candidate, str):
            text = candidate.strip()
            return text or None
    return None


def _coerce_ligand_value(value: Any, preferred_key: str) -> str | None:
    text = _coerce_string(value, preferred_key)
    if preferred_key == "csv" and text:
        candidate_path = Path(text).expanduser()
        if candidate_path.is_file():
            return candidate_path.read_text(encoding="utf-8")
    return text


def _coerce_mapping(value: Any) -> dict[str, Any]:
    value = _signal_value(value)
    if isinstance(value, Mapping):
        return {str(key): item for key, item in value.items()}
    return {}


class BoltzInputAssemblerModel(BioModule):
    """Prepare exact Boltz inputs from curated defaults and optional public overrides."""

    def __init__(
        self,
        workflow_name: str,
        workflow_kind: str = "single",
        default_protein_sequence: str | None = None,
        default_ligand_smiles: str | None = None,
        default_ligand_csv: str | None = None,
        default_msa_path: str | None = None,
        default_run_options: Mapping[str, Any] | None = None,
        integration_step: float = 0.01,
    ) -> None:
        self.integration_step = float(integration_step)
        self.workflow_name = workflow_name
        self.workflow_kind = workflow_kind
        self.default_protein_sequence = (default_protein_sequence or "").strip() or None
        self.default_ligand_smiles = (default_ligand_smiles or "").strip() or None
        self.default_ligand_csv = (default_ligand_csv or "").strip() or None
        self.default_msa_path = (default_msa_path or "").strip() or None
        self.default_run_options = dict(default_run_options or {})
        self._inputs: dict[str, BioSignal] = {}
        self._outputs: dict[str, BioSignal] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        ligand_name = "ligand_csv" if self.workflow_kind == "batch" else "ligand_smiles"
        ligand_description = (
            "Candidate ligand CSV text passed to the Boltz batch ranker."
            if self.workflow_kind == "batch"
            else "Ligand SMILES string passed to the Boltz-2 predictor."
        )
        return {
            "scenario_context": _record_input_spec("Curated workflow context emitted by the context stage."),
            "protein_sequence": _record_input_spec("Optional target protein sequence override."),
            ligand_name: _record_input_spec(ligand_description),
            "msa_path": _record_input_spec("Optional precomputed MSA path override."),
            "run_options": _record_input_spec("Optional structured Boltz run-option overrides."),
        }

    def outputs(self) -> dict[str, SignalSpec]:
        ligand_name = "ligand_csv" if self.workflow_kind == "batch" else "ligand_smiles"
        ligand_description = (
            "Resolved candidate ligand CSV for the Boltz batch ranker."
            if self.workflow_kind == "batch"
            else "Resolved ligand SMILES for the Boltz-2 predictor."
        )
        return {
            "protein_sequence": SignalSpec.scalar(dtype="str", description="Resolved amino-acid sequence passed to Boltz."),
            ligand_name: SignalSpec.scalar(dtype="str", description=ligand_description),
            "msa_path": SignalSpec.scalar(dtype="str", description="Resolved precomputed MSA path, or an empty string when the MSA server is used."),
            "run_options": SignalSpec.record(schema={"payload": "json"}, description="Merged Boltz run options passed to the prediction stage."),
            "assembled_boltz_request": SignalSpec.record(
                schema={"payload": "json"},
                description="Traceable non-secret summary of the exact Boltz request assembled for this workflow.",
            ),
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
        protein_sequence = _coerce_string(self._inputs.get("protein_sequence"), "sequence") or self.default_protein_sequence or ""
        ligand_name = "ligand_csv" if self.workflow_kind == "batch" else "ligand_smiles"
        ligand_key = "csv" if self.workflow_kind == "batch" else "smiles"
        ligand_value = (
            _coerce_ligand_value(self._inputs.get(ligand_name), ligand_key)
            or (self.default_ligand_csv if self.workflow_kind == "batch" else self.default_ligand_smiles)
            or ""
        )
        msa_path = _coerce_string(self._inputs.get("msa_path"), "path") or self.default_msa_path or ""
        run_options = dict(self.default_run_options)
        run_options.update(_coerce_mapping(self._inputs.get("run_options")))
        if context:
            run_options.setdefault("workflow_context", context.get("workflow_context"))
            run_options.setdefault("interpretation_scope", context.get("interpretation_scope"))
            run_options.setdefault("target_name", context.get("target_name"))
            run_options.setdefault("ligand_name", context.get("ligand_name"))

        assembled = {
            "workflow_name": self.workflow_name,
            "workflow_kind": self.workflow_kind,
            "target_name": context.get("target_name") or run_options.get("target_name"),
            "ligand_name": context.get("ligand_name") or run_options.get("ligand_name"),
            "protein_sequence_length": len(protein_sequence),
            "protein_sequence_sha256": hashlib.sha256(protein_sequence.encode("utf-8")).hexdigest() if protein_sequence else "",
            "msa_path_supplied": bool(msa_path),
            "run_option_keys": sorted(str(key) for key in run_options.keys()),
            "scientific_caveat": context.get("caveat"),
        }
        if self.workflow_kind == "batch":
            assembled["ligand_count"] = self._ligand_count(ligand_value)
            assembled["ligand_library_sha256"] = hashlib.sha256(ligand_value.encode("utf-8")).hexdigest() if ligand_value else ""
        else:
            assembled["ligand_smiles"] = ligand_value
            assembled["ligand_smiles_sha256"] = hashlib.sha256(ligand_value.encode("utf-8")).hexdigest() if ligand_value else ""

        source = getattr(self, "_world_name", self.__class__.__name__)
        specs = self.outputs()
        self._outputs = {
            "protein_sequence": make_signal(source=source, name="protein_sequence", value=protein_sequence, emitted_at=emitted_at, spec=specs["protein_sequence"]),
            ligand_name: make_signal(source=source, name=ligand_name, value=ligand_value, emitted_at=emitted_at, spec=specs[ligand_name]),
            "msa_path": make_signal(source=source, name="msa_path", value=msa_path, emitted_at=emitted_at, spec=specs["msa_path"]),
            "run_options": make_signal(source=source, name="run_options", value=run_options, emitted_at=emitted_at, spec=specs["run_options"]),
            "assembled_boltz_request": make_signal(
                source=source,
                name="assembled_boltz_request",
                value=assembled,
                emitted_at=emitted_at,
                spec=specs["assembled_boltz_request"],
            ),
        }
        return dict(self._outputs)

    def get_outputs(self) -> dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> list[dict[str, Any]] | None:
        return None

    @staticmethod
    def _ligand_count(csv_text: str) -> int:
        if not csv_text.strip():
            return 0
        try:
            rows = list(csv.DictReader(io.StringIO(csv_text)))
        except csv.Error:
            return 0
        return len([row for row in rows if any((value or "").strip() for value in row.values())])
