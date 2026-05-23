# SPDX-FileCopyrightText: 2026-present Biosimulant Team
# SPDX-License-Identifier: Apache-2.0
"""Dedicated visualisation model for docking and Boltz workflow labs."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

from biosim import BioModule
from biosim.signals import AcceptedSignalProfile, BioSignal, SignalSpec
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


class DockingVisualisationModel(BioModule):
    def __init__(
        self,
        integration_step: float = 0.01,
        source_alias: str = "core",
        mode: str = "vina",
        lab_title: str = "Docking Lab",
        context_alias: str | None = None,
        assembler_alias: str | None = None,
        interpreter_alias: str | None = None,
    ) -> None:
        self.integration_step = float(integration_step)
        self.source_alias = source_alias
        self.mode = mode
        self.lab_title = lab_title
        self.context_alias = context_alias
        self.assembler_alias = assembler_alias
        self.interpreter_alias = interpreter_alias
        self._inputs: Dict[str, BioSignal] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        names_by_mode = {
            "vina": ["pose_summary", "docking_summary", "structure_artifacts", "run_metadata"],
            "boltz": ["affinity_summary", "confidence_summary", "structure_artifacts", "run_metadata"],
            "boltz_batch": ["batch_summary", "affinity_summary", "confidence_summary", "structure_artifacts", "run_metadata"],
            "diffdock": ["pose_summary", "confidence_summary", "structure_artifacts", "run_metadata"],
        }
        specs = {
            f"{self.source_alias}_{name}": _record_input_spec(f"Internal {name} input from the sibling core model.")
            for name in names_by_mode[self.mode]
        }
        if self.context_alias:
            specs[f"{self.context_alias}_scenario_context"] = _record_input_spec("Workflow target and caveat context.")
        if self.assembler_alias:
            specs[f"{self.assembler_alias}_assembled_boltz_request"] = _record_input_spec("Resolved Boltz request summary.")
        if self.interpreter_alias:
            specs[f"{self.interpreter_alias}_prediction_evidence"] = _record_input_spec("Conservative interpreted prediction evidence.")
        return specs

    def outputs(self) -> dict[str, SignalSpec]:
        return {}

    def reset(self) -> None:
        self._inputs = {}

    def set_inputs(self, signals: dict[str, BioSignal]) -> None:
        self._inputs.update(signals or {})

    def advance_window(self, start: float | None = None, end: float | None = None, inputs: dict[str, BioSignal] | None = None) -> dict[str, BioSignal]:
        if inputs:
            self.set_inputs(inputs)
        return {}

    def get_outputs(self) -> dict[str, BioSignal]:
        return {}

    def visualize(self) -> Optional[list[dict[str, Any]]]:
        if self.mode == "vina":
            primary = self._visualize_vina()
        elif self.mode == "boltz":
            primary = self._visualize_boltz()
        elif self.mode == "boltz_batch":
            primary = self._visualize_boltz_batch()
        else:
            primary = self._visualize_diffdock()
        workflow = self._workflow_visuals()
        if primary and workflow:
            return [*primary, *workflow]
        if primary:
            return primary
        return workflow or None

    def _input_value(self, name: str) -> Any:
        return _signal_value(self._inputs.get(f"{self.source_alias}_{name}"))

    def _stage_value(self, alias: str | None, name: str) -> Any:
        if not alias:
            return None
        return _signal_value(self._inputs.get(f"{alias}_{name}"))

    def _workflow_visuals(self) -> list[dict[str, Any]]:
        visuals: list[dict[str, Any]] = []
        context = self._stage_value(self.context_alias, "scenario_context")
        assembled = self._stage_value(self.assembler_alias, "assembled_boltz_request")
        evidence = self._stage_value(self.interpreter_alias, "prediction_evidence")
        if isinstance(evidence, Mapping):
            rows = [
                ["Scientific question", str(evidence.get("scientific_question") or "")],
                ["Observed answer", str(evidence.get("observed_answer") or "")],
                ["Dominant module", str(evidence.get("dominant_module") or self.source_alias)],
                ["Run status", str(evidence.get("run_status") or "")],
                ["Caveat", str(evidence.get("caveat") or "")],
            ]
            visuals.append(
                {
                    "render": "table",
                    "description": "Conservative interpretation of the current Boltz workflow run.",
                    "data": {"title": f"{self.lab_title} - workflow answer", "columns": ["Prompt", "Answer"], "rows": rows},
                }
            )
        if isinstance(context, Mapping):
            fields = [
                ("Target", context.get("target_name")),
                ("Target family", context.get("target_family")),
                ("Disease or use case", context.get("disease_area")),
                ("Ligand", context.get("ligand_name")),
                ("Ligand role", context.get("ligand_role")),
                ("Source PDB", context.get("source_pdb")),
                ("PubChem CID", context.get("source_pubchem_cid")),
                ("Scope", context.get("interpretation_scope")),
            ]
            rows = [[label, str(value)] for label, value in fields if value not in (None, "")]
            if rows:
                visuals.append(
                    {
                        "render": "table",
                        "description": "Source-backed target, ligand, and use-case context for this workflow.",
                        "data": {"title": "Workflow target and ligand context", "columns": ["Field", "Value"], "rows": rows},
                    }
                )
        if isinstance(assembled, Mapping):
            fields = [
                ("Workflow kind", assembled.get("workflow_kind")),
                ("Protein sequence length", assembled.get("protein_sequence_length")),
                ("Protein sequence SHA-256", assembled.get("protein_sequence_sha256")),
                ("Ligand count", assembled.get("ligand_count")),
                ("Ligand SMILES SHA-256", assembled.get("ligand_smiles_sha256")),
                ("Ligand library SHA-256", assembled.get("ligand_library_sha256")),
                ("MSA path supplied", assembled.get("msa_path_supplied")),
            ]
            rows = [[label, str(value)] for label, value in fields if value not in (None, "")]
            if rows:
                visuals.append(
                    {
                        "render": "table",
                        "description": "Traceability summary for the assembled Boltz request.",
                        "data": {"title": "Assembled Boltz request evidence", "columns": ["Field", "Value"], "rows": rows},
                    }
                )
        return visuals

    def _visualize_vina(self) -> Optional[list[dict[str, Any]]]:
        run_metadata = self._input_value("run_metadata")
        artifacts = self._input_value("structure_artifacts")
        docking_summary = self._input_value("docking_summary")
        poses = self._input_value("pose_summary")
        if not isinstance(run_metadata, Mapping) or run_metadata.get("status") != "completed":
            return None
        if not isinstance(artifacts, Mapping) or not isinstance(docking_summary, Mapping) or not isinstance(poses, list):
            return None
        top_complex_path = self._resolved_path(artifacts.get("top_complex_file"))
        if top_complex_path is None:
            return None
        rows = []
        for row in poses:
            if not isinstance(row, Mapping):
                continue
            rows.append([
                str(row.get("rank", "")),
                "" if row.get("affinity_kcal_mol") is None else str(row.get("affinity_kcal_mol")),
                "" if row.get("rmsd_lb") is None else str(row.get("rmsd_lb")),
                "" if row.get("rmsd_ub") is None else str(row.get("rmsd_ub")),
                Path(str(row.get("pose_file", ""))).name,
            ])
        return [
            {
                "render": "structure3d",
                "description": "Top-ranked AutoDock Vina complex for the latest docking run.",
                "data": {
                    "title": "Top-Ranked Docked Complex",
                    "source": {"kind": "artifact", "artifact_id": self._artifact_id(top_complex_path), "path": str(top_complex_path)},
                    "format": "pdb",
                    "annotations": [
                        {"label": "Top Pose Affinity (kcal/mol)", "value": docking_summary.get("top_pose_affinity_kcal_mol")},
                        {"label": "Scoring", "value": docking_summary.get("scoring")},
                        {"label": "Pose Count", "value": docking_summary.get("pose_count")},
                    ],
                    "initial_view": {"reset_camera": True},
                },
            },
            {
                "render": "table",
                "description": "Ranked pose summary from the latest AutoDock Vina run.",
                "data": {"title": "AutoDock Vina Pose Summary", "columns": ["Rank", "Affinity", "RMSD l.b.", "RMSD u.b.", "Pose File"], "rows": rows},
            },
        ]

    def _visualize_boltz(self) -> Optional[list[dict[str, Any]]]:
        run_metadata = self._input_value("run_metadata")
        artifacts = self._input_value("structure_artifacts")
        confidence = self._input_value("confidence_summary")
        affinity = self._input_value("affinity_summary")
        if not isinstance(run_metadata, Mapping) or run_metadata.get("status") != "completed":
            return None
        if not isinstance(artifacts, Mapping):
            return None
        structure_path = self._resolved_path(artifacts.get("structure_file"))
        if structure_path is None:
            return None
        structure_format = self._structure_format(structure_path)
        if structure_format is None:
            return None
        annotations = self._build_boltz_annotations(confidence, affinity)
        return [
            {
                "render": "structure3d",
                "description": "Top-ranked Boltz structure prediction for the latest protein-ligand run.",
                "data": {
                    "title": "Predicted Complex Structure",
                    "source": {"kind": "artifact", "artifact_id": self._artifact_id(structure_path), "path": str(structure_path)},
                    "format": structure_format,
                    "annotations": [{"label": label, "value": value} for label, value in annotations],
                    "initial_view": {"reset_camera": True},
                },
            },
            {
                "render": "table",
                "description": "Key affinity and confidence metrics extracted from the latest Boltz outputs.",
                "data": {"title": "Boltz Summary", "columns": ["Metric", "Value"], "rows": [[label, str(value)] for label, value in annotations]},
            },
        ]

    def _visualize_boltz_batch(self) -> Optional[list[dict[str, Any]]]:
        run_metadata = self._input_value("run_metadata")
        artifacts = self._input_value("structure_artifacts")
        confidence = self._input_value("confidence_summary")
        affinity = self._input_value("affinity_summary")
        batch = self._input_value("batch_summary")
        if not isinstance(run_metadata, Mapping) or run_metadata.get("status") != "completed":
            return None
        if not isinstance(artifacts, Mapping) or not isinstance(batch, Mapping):
            return None
        structure_path = self._resolved_path(artifacts.get("structure_file"))
        if structure_path is None:
            return None
        structure_format = self._structure_format(structure_path)
        if structure_format is None:
            return None
        annotations = self._build_boltz_annotations(confidence, affinity)
        rows = []
        ranked = batch.get("ranked_ligands")
        if isinstance(ranked, list):
            for row in ranked:
                if not isinstance(row, Mapping):
                    continue
                rows.append(
                    [
                        str(row.get("rank") or ""),
                        str(row.get("ligand") or ""),
                        "" if row.get("binder_probability") is None else str(row.get("binder_probability")),
                        "" if row.get("affinity_like_value") is None else str(row.get("affinity_like_value")),
                        "" if row.get("confidence") is None else str(row.get("confidence")),
                        "; ".join(str(item) for item in row.get("flags", []) if item),
                    ]
                )
        return [
            {
                "render": "structure3d",
                "description": "Top-ranked Boltz-2 complex from the latest batch ligand ranking run.",
                "data": {
                    "title": "Top-Ranked Batch Complex",
                    "source": {"kind": "artifact", "artifact_id": self._artifact_id(structure_path), "path": str(structure_path)},
                    "format": structure_format,
                    "annotations": [{"label": label, "value": value} for label, value in annotations],
                    "initial_view": {"reset_camera": True},
                },
            },
            {
                "render": "table",
                "description": "Ranked ligand table from the latest Boltz-2 batch run.",
                "data": {
                    "title": "Batch Ligand Ranking",
                    "columns": ["Rank", "Ligand", "Binder Probability", "Affinity-Like Value", "Confidence", "Flags"],
                    "rows": rows,
                },
            },
        ]

    def _visualize_diffdock(self) -> Optional[list[dict[str, Any]]]:
        run_metadata = self._input_value("run_metadata")
        artifacts = self._input_value("structure_artifacts")
        confidence = self._input_value("confidence_summary")
        poses = self._input_value("pose_summary")
        if not isinstance(run_metadata, Mapping) or run_metadata.get("status") != "completed":
            return None
        if not isinstance(artifacts, Mapping) or not isinstance(poses, list):
            return None
        top_complex_path = self._resolved_path(artifacts.get("top_complex_file"))
        if top_complex_path is None:
            return None
        rows = []
        for row in poses:
            if not isinstance(row, Mapping):
                continue
            rows.append([
                str(row.get("rank", "")),
                "" if row.get("confidence") is None else str(row.get("confidence")),
                str(row.get("confidence_band") or ""),
                Path(str(row.get("file_path", ""))).name,
            ])
        return [
            {
                "render": "structure3d",
                "description": "Top-ranked DiffDock-L complex for the latest docking run.",
                "data": {
                    "title": "Top-Ranked Docked Complex",
                    "source": {"kind": "artifact", "artifact_id": self._artifact_id(top_complex_path), "path": str(top_complex_path)},
                    "format": "pdb",
                    "annotations": [
                        {"label": "Top Pose Confidence", "value": confidence.get("top_pose_confidence") if isinstance(confidence, Mapping) else None},
                        {"label": "Confidence Band", "value": confidence.get("confidence_band") if isinstance(confidence, Mapping) else None},
                        {"label": "Pose Count", "value": confidence.get("pose_count") if isinstance(confidence, Mapping) else None},
                    ],
                    "initial_view": {"reset_camera": True},
                },
            },
            {
                "render": "table",
                "description": "Ranked pose summary from the latest DiffDock-L run.",
                "data": {"title": "DiffDock Pose Summary", "columns": ["Rank", "Confidence", "Band", "Pose File"], "rows": rows},
            },
        ]

    @staticmethod
    def _artifact_id(path: Path) -> str:
        return hashlib.sha1(str(path).encode("utf-8")).hexdigest()

    @staticmethod
    def _resolved_path(value: Any) -> Optional[Path]:
        if not isinstance(value, str) or not value:
            return None
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        return path.resolve()

    @staticmethod
    def _structure_format(path: Path) -> Optional[str]:
        suffix = path.suffix.lower()
        if suffix in {".cif", ".mmcif"}:
            return "mmcif"
        if suffix in {".pdb", ".ent"}:
            return "pdb"
        return None

    @staticmethod
    def _build_boltz_annotations(confidence: Any, affinity: Any) -> list[tuple[str, Any]]:
        confidence = confidence if isinstance(confidence, Mapping) else {}
        affinity = affinity if isinstance(affinity, Mapping) else {}
        candidates = [
            ("Affinity-like value", affinity.get("affinity_pred_value", affinity.get("predicted_affinity_value"))),
            ("Affinity unit", affinity.get("affinity_unit", affinity.get("predicted_affinity_unit"))),
            ("Binder probability", affinity.get("affinity_probability_binary")),
            ("Confidence score", confidence.get("confidence_score")),
            ("pTM", confidence.get("ptm", confidence.get("mean_ptm"))),
            ("ipTM", confidence.get("iptm", confidence.get("mean_iptm"))),
            ("Complex pLDDT", confidence.get("complex_plddt", confidence.get("mean_plddt"))),
            ("Ligand ipTM", confidence.get("ligand_iptm", confidence.get("mean_ligand_iptm"))),
        ]
        return [(label, value) for label, value in candidates if value not in (None, "")]
