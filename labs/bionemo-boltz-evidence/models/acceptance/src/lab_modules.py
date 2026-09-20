"""Finite BioSimulant Lab for captured BioNeMo results and enforced ranking checks."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from biosimulant import BioModule, ExecutionPolicy, SignalSpec, unwrap_payload
from evidence import evaluate, digest


class ResultSource(BioModule):
    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN
    def __init__(self, task, records=None):
        self.task = task
        self.records = records or []
    def inputs(self):
        return {}
    def outputs(self):
        return {"provider_results": SignalSpec.record(schema={"payload": "json"})}
    def execute(self, inputs, *, context):
        return {"provider_results": {"task": self.task, "records": self.records}}


class Acceptance(BioModule):
    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN
    def inputs(self):
        return {"provider_results": SignalSpec.record(schema={"payload": "json"})}
    def outputs(self):
        return {"acceptance": SignalSpec.record(schema={"payload": "json"})}
    def execute(self, inputs, *, context):
        data = unwrap_payload(inputs["provider_results"])
        report = evaluate(data["task"], data["records"])
        return {"acceptance": report}


class AcceptedRanking(BioModule):
    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN
    def inputs(self):
        return {"acceptance": SignalSpec.record(schema={"payload": "json"})}
    def outputs(self):
        return {"report": SignalSpec.record(schema={"payload": "json"}),
                "accepted": SignalSpec.scalar(dtype="bool")}
    def execute(self, inputs, *, context):
        report = unwrap_payload(inputs["acceptance"])
        accepted = report["status"] == "passed"
        return {"accepted": accepted, "report": {**report,
                "ranking": report["ranking"] if accepted else [], "report_sha256": digest(report),
                "execution_boundary": "NVIDIA inference runs in the authorized agent environment. This Lab validates captured responses and gates downstream ranking."}}
