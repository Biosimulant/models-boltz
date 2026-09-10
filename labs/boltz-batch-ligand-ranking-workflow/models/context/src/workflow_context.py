# SPDX-FileCopyrightText: 2026-present Biosimulant Team
# SPDX-License-Identifier: Apache-2.0
"""Source-backed scenario context stage for Boltz workflow labs."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from biosim import BioModule, ExecutionContext, ExecutionPolicy
from biosim.signals import BioSignal, SignalSpec, make_signal


class WorkflowContextModel(BioModule):
    """Emit curated target/ligand provenance and caveats for a Boltz workflow."""

    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN

    def __init__(self, scenario: Mapping[str, Any] | None = None, integration_step: float = 0.01) -> None:
        self.integration_step = float(integration_step)
        self.scenario = dict(scenario or {})
        self._outputs: dict[str, BioSignal] = {}

    def inputs(self) -> dict[str, SignalSpec]:
        return {}

    def outputs(self) -> dict[str, SignalSpec]:
        return {
            "scenario_context": SignalSpec.record(
                schema={"payload": "json"},
                description="Curated source-backed target, ligand, workflow scope, and caveat context.",
            )
        }

    def reset(self) -> None:
        super().reset()
        self._outputs = {}

    def set_inputs(self, signals: dict[str, BioSignal]) -> None:
        return None

    def execute(self, inputs: Mapping[str, BioSignal], *, context: ExecutionContext) -> Mapping[str, BioSignal]:
        self.set_inputs(dict(inputs))
        result = self._execute_at_time(0.0, 0.0)
        return dict(result if result is not None else getattr(self, "_outputs", {}))

    def _execute_at_time(
        self,
        start: float | None = None,
        end: float | None = None,
        inputs: dict[str, BioSignal] | None = None,
    ) -> dict[str, BioSignal]:
        emitted_at = float(end if end is not None else self.integration_step)
        payload = dict(self.scenario)
        payload.setdefault("scientific_status", "source-backed workflow context")
        payload.setdefault(
            "caveat",
            "Boltz-2 outputs are computational predictions for hypothesis generation; they are not experimental binding, potency, selectivity, clinical, or efficacy evidence.",
        )
        source = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "scenario_context": make_signal(
                source=source,
                name="scenario_context",
                value=payload,
                emitted_at=emitted_at,
                spec=self.outputs()["scenario_context"],
            )
        }
        return dict(self._outputs)

    def visualize(self) -> list[dict[str, Any]] | None:
        return None
