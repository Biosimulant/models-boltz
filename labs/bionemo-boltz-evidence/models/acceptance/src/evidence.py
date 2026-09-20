"""Provider-neutral inspection of actual BioNeMo Boltz request/response records.

Shared with the native baseline. Hashes establish byte identity, not provider
authentication. Acceptance covers declared technical checks, not drug efficacy.
"""
import hashlib
import json
import math
import tempfile
from pathlib import Path

from structure_checks import inspect_structure

SCHEMA = "biosimulant.bionemo-boltz-evidence/v1"
LIMITATIONS = [
    "Research predictions, not experimental binding or clinical evidence.",
    "Ligand heavy-element inventory does not establish connectivity, stereochemistry or pose accuracy.",
    "NIM seed, checkpoint hash and hosting hardware are not exposed.",
    "Receipt hashes identify captured bytes; NVIDIA has not cryptographically signed these responses.",
    "Query-only MSA supplies no evolutionary homologs.",
    "Native affinity score and provider affinity_pic50 remain separate fields.",
]


def digest(value):
    raw = value if isinstance(value, bytes) else json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def request_for(bundle, ligand):
    return {
        "polymers": [{"id": "A", "molecule_type": "protein", "sequence": bundle["sequence"],
                      "msa": {"benchmark_query": {"a3m": {"alignment": bundle["msa"], "format": "a3m", "rank": 0}}}}],
        "ligands": [{"id": "B", "smiles": ligand["smiles"], "predict_affinity": True}],
        **bundle["settings"], "step_scale": 1.5, "output_format": "mmcif",
        "without_potentials": True, "affinity_mw_correction": False,
    }


def task_for(bundle, names):
    candidates = [next(x for x in bundle["ligands"] if x["name"] == n) for n in names]
    if len(names) != len(set(names)):
        raise ValueError("Candidate names must be unique")
    msa = "".join(line.strip() for line in bundle["msa"].splitlines() if not line.startswith(">"))
    if msa != bundle["sequence"]:
        raise ValueError("MSA query does not match requested sequence")
    return {"schema": SCHEMA, "target": "UniProt:P62942", "sequence": bundle["sequence"],
            "candidates": [{**x, "request": request_for(bundle, x)} for x in candidates],
            "acceptance": "Exact candidate coverage, request and response hashes, HTTP 200, one structure, finite scores/confidence, bounded probability, exact protein sequence, ligand heavy-element inventory.",
            "limitations": LIMITATIONS}


def scalar(value, label):
    if isinstance(value, list):
        if len(value) != 1:
            raise ValueError(label + " requires exactly one value")
        value = value[0]
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(label + " must be a finite number")
    return value


def validate_record(candidate, record):
    checks = []
    def require(condition, name):
        if not condition:
            raise ValueError(name)
        checks.append(name)
    require(record["candidate"] == candidate["name"], "candidate_identity")
    require(record["request"] == candidate["request"], "request_matches_frozen_task")
    require(record["request_sha256"] == digest(record["request"]), "request_hash")
    raw = record["response_text"].encode()
    require(record["response_sha256"] == digest(raw), "response_hash")
    require(record["http_status"] == 200, "provider_http_200")
    response = json.loads(raw)
    require(len(response.get("structures", [])) == 1, "one_structure")
    structure = response["structures"][0]["structure"]
    require(isinstance(structure, str) and bool(structure.strip()), "nonempty_structure")
    values = response["affinities"]["B"]
    score = scalar(values["affinity_pred_value"], "raw_affinity")
    probability = scalar(values["affinity_probability_binary"], "binder_probability")
    require(0 <= probability <= 1, "bounded_binder_probability")
    confidence = scalar(response["confidence_scores"], "aggregate_confidence")
    require(0 <= confidence <= 1, "bounded_aggregate_confidence")
    checks.append("finite_affinity")
    with tempfile.TemporaryDirectory(prefix="boltz-inspect-") as directory:
        path = Path(directory) / "prediction.cif"
        path.write_text(structure)
        structural = inspect_structure(str(path), candidate["request"]["polymers"][0]["sequence"], candidate["smiles"])
    checks.extend(["finite_coordinates", "exact_target_sequence", "ligand_heavy_element_inventory", "one_model_unique_atom_ids"])
    return {"candidate": candidate["name"], "cid": candidate["cid"], "raw_affinity_score": score,
            "binder_probability": probability, "aggregate_confidence": confidence,
            "provider_affinity_pic50": scalar(values["affinity_pic50"], "provider_pic50") if "affinity_pic50" in values else None,
            "structure_sha256": digest(structure.encode()), "response_sha256": record["response_sha256"],
            "origin_execution_id": record["execution_id"], "provenance_mode": record["provenance_mode"],
            "structure_checks": structural, "passed_checks": checks}


def evaluate(task, records):
    """Recompute acceptance from source records; never trust stored pass flags."""
    expected = [x["name"] for x in task["candidates"]]
    task_errors = []
    if task.get("schema") != SCHEMA or not expected or len(expected) != len(set(expected)):
        task_errors.append({"check": "task_contract", "message": "Recognized schema and nonempty unique candidates required"})
    for candidate in task["candidates"]:
        request = candidate.get("request", {})
        polymers, ligands = request.get("polymers", []), request.get("ligands", [])
        if (len(polymers) != 1 or len(ligands) != 1 or polymers[0].get("id") != "A"
            or ligands[0].get("id") != "B" or polymers[0].get("sequence") != task.get("sequence")
            or ligands[0].get("smiles") != candidate.get("smiles")):
            task_errors.append({"check": "task_request_consistency", "message": "Declared task must match the exact provider request"})
    observed = [x.get("candidate") for x in records]
    failures, rows = task_errors, []
    if len(observed) != len(set(observed)):
        failures.append({"check": "unique_candidates", "message": "Duplicate candidate records"})
    if set(observed) != set(expected) or len(observed) != len(expected):
        failures.append({"check": "candidate_coverage", "message": f"Expected {len(expected)} predictions, received {len(observed)}", "missing": sorted(set(expected) - set(observed))})
    for candidate in task["candidates"]:
        matching = [x for x in records if x.get("candidate") == candidate["name"]]
        if len(matching) != 1:
            continue
        try:
            rows.append(validate_record(candidate, matching[0]))
        except (ValueError, KeyError, TypeError, IndexError, OSError) as exc:
            failures.append({"candidate": candidate["name"], "check": "record_acceptance", "message": str(exc)})
    passed = not failures and len(rows) == len(expected)
    return {"schema": SCHEMA, "task_sha256": digest(task), "status": "passed" if passed else "failed",
            "scope": "declared_technical_requirements", "expected_count": len(expected), "received_count": len(records),
            "accepted_count": len(rows), "failures": failures, "rows": rows,
            "ranking": [r["candidate"] for r in sorted(rows, key=lambda x: x["raw_affinity_score"])] if passed else [],
            "ranking_status": "available" if passed else "blocked_by_acceptance", "limitations": LIMITATIONS}


def record_from_files(candidate, request_path, response_path, execution_id, mode, http_status=200):
    request = json.loads(Path(request_path).read_text())
    raw = Path(response_path).read_bytes()
    return {"candidate": candidate, "execution_id": execution_id, "provenance_mode": mode,
            "request": request, "request_sha256": digest(request), "response_text": raw.decode(),
            "response_sha256": digest(raw), "http_status": http_status}
