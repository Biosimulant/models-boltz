"""Custom BioModule implementations.

These modules can be referenced in YAML configs:
    modules:
      counter:
        class: my_pack.Counter
        args:
          name: "step_counter"
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

from biosim import BioModule
from biosim.signals import BioSignal, SignalMetadata

import logging

logger = logging.getLogger(__name__)


class Counter(BioModule):
    """Counts simulation steps and emits the count.

    Outputs:
        count: {"count": int, "t": float}

    Parameters:
        name: Display name for this counter
    """

    def __init__(self, name: str = "counter", min_dt: float = 0.1) -> None:
        self.min_dt = min_dt
        self.name = name
        self._count = 0
        self._history: List[List[float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return set()

    def outputs(self) -> Set[str]:
        return {"count"}

    def reset(self) -> None:
        self._count = 0
        self._history = []

    def advance_to(self, t: float) -> None:
        self._count += 1
        self._history.append([t, self._count])
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "count": BioSignal(
                source=source_name,
                name="count",
                value={"count": self._count, "t": t},
                time=t,
                metadata=SignalMetadata(units=None, description="Step count", kind="state"),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[Dict[str, Any]]:
        if not self._history:
            return None
        return {
            "render": "timeseries",
            "data": {
                "series": [{"name": self.name, "points": self._history}],
                "title": f"Counter: {self.name}",
            },
        }


class Accumulator(BioModule):
    """Receives values and accumulates them over time.

    Inputs:
        value: float - adds amount to total

    Outputs:
        total: {"total": float, "t": float}

    Parameters:
        initial: Starting value for accumulator
    """

    def __init__(self, initial: float = 0.0, min_dt: float = 0.1) -> None:
        self.min_dt = min_dt
        self._initial = initial
        self._total = initial
        self._history: List[List[float]] = []
        self._outputs: Dict[str, BioSignal] = {}

    def inputs(self) -> Set[str]:
        return {"value"}

    def outputs(self) -> Set[str]:
        return {"total"}

    def reset(self) -> None:
        self._total = self._initial
        self._history = []

    def set_inputs(self, signals: Dict[str, BioSignal]) -> None:
        signal = signals.get("value")
        if signal is None:
            return
        try:
            self._total += float(signal.value)
        except (KeyError, ValueError, TypeError):  # narrowed from bare Exception
            pass

    def advance_to(self, t: float) -> None:
        self._history.append([t, self._total])
        source_name = getattr(self, "_world_name", self.__class__.__name__)
        self._outputs = {
            "total": BioSignal(
                source=source_name,
                name="total",
                value={"total": self._total, "t": t},
                time=t,
                metadata=SignalMetadata(units=None, description="Accumulated total", kind="state"),
            )
        }

    def get_outputs(self) -> Dict[str, BioSignal]:
        return dict(self._outputs)

    def visualize(self) -> Optional[Dict[str, Any]]:
        if not self._history:
            return None
        return {
            "render": "timeseries",
            "data": {
                "series": [{"name": "total", "points": self._history}],
                "title": "Accumulator",
            },
        }
