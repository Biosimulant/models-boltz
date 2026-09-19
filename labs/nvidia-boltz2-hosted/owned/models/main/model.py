"""BIO-001: fresh NVIDIA Boltz-2 inference, never captured-response replay."""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import re
import time
import uuid
from datetime import datetime, timezone

import httpx
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from biosim import BioModule, ExecutionPolicy, SignalSpec

ENDPOINT = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
POLL_ENDPOINT = "https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/"
SYSTEM_SCHEMA = {"polymers": "json", "ligands": "json", "constraints": "json"}
OPTION_SCHEMA = {"parameters": "json"}
STRUCTURE_SCHEMA = {"items": "json", "coordinate_unit": "str", "request_sha256": "str"}
AFFINITY_SCHEMA = {"by_ligand": "json", "field_definitions": "json", "request_sha256": "str"}
QUALITY_SCHEMA = {"values": "json", "field_units": "json", "request_sha256": "str"}
PROVENANCE_SCHEMA = {
    "request": "json", "request_sha256": "str", "response_sha256": "str",
    "response_bytes": "int", "provider_request_id": "str", "endpoint": "str",
    "captured_at": "str", "wall_seconds": "float", "resumed": "bool",
    "provider_headers": "json", "backend_identity": "json", "limitations": "json",
}
INT_OPTIONS = {"recycling_steps": (1, 10), "sampling_steps": (10, 1000),
               "diffusion_samples": (1, 25), "sampling_steps_affinity": (10, 1000),
               "diffusion_samples_affinity": (1, 10)}
BOOL_OPTIONS = {"without_potentials", "concatenate_msas", "affinity_mw_correction"}
QUALITY_FIELDS = (
    "confidence_scores", "ptm_scores", "iptm_scores", "ligand_iptm_scores",
    "protein_iptm_scores", "complex_plddt_scores", "complex_iplddt_scores",
    "complex_pde_scores", "complex_ipde_scores", "chains_ptm_scores",
    "pair_chains_iptm_scores", "pae", "pde",
)


class ProviderObservationError(RuntimeError):
    """Safe diagnostic; a failure to observe does not cancel provider compute."""


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _keys(value, required, allowed, label):
    require(isinstance(value, dict), label + " must be an object")
    require(set(required) <= set(value) <= set(allowed), label + " has missing or unknown fields")


def _positive_index(value, size):
    return type(value) is int and 1 <= value <= size


def _a3m_query(text):
    require(isinstance(text, str) and text.startswith(">"), "A3M must start with a FASTA header")
    sequences = []
    for line in text.splitlines():
        if line.startswith(">"):
            sequences.append("")
        elif line.strip():
            require(bool(sequences), "A3M sequence has no header")
            require(re.fullmatch(r"[A-Za-z.\-]+", line) is not None, "Invalid A3M sequence characters")
            sequences[-1] += line
    require(bool(sequences) and all(sequences), "Empty A3M sequence")
    columns = [re.sub(r"[a-z.]", "", s) for s in sequences]
    require(len({len(s) for s in columns}) == 1, "A3M alignment columns disagree")
    return columns[0].replace("-", "")


def prepare_request(system, options):
    """Validate the molecular request without changing scientific payloads."""
    _keys(system, SYSTEM_SCHEMA, SYSTEM_SCHEMA, "molecular_system")
    _keys(options, OPTION_SCHEMA, OPTION_SCHEMA, "run_options")
    require(all(isinstance(system[k], list) for k in SYSTEM_SCHEMA), "System fields must be lists")
    polymers, ligands, constraints = (system[k] for k in SYSTEM_SCHEMA)
    require(1 <= len(polymers) <= 12, "Require 1 to 12 polymers")
    require(len(ligands) <= 20, "At most 20 ligands")
    identities = set()
    polymer_by_id = {}
    for polymer in polymers:
        _keys(polymer, {"id", "molecule_type", "sequence"},
              {"id", "molecule_type", "sequence", "cyclic", "msa", "modifications", "structural_templates"}, "polymer")
        molecule_type = polymer["molecule_type"]
        alphabet = {"protein": "ACDEFGHIKLMNPQRSTVWYX", "dna": "ACGT", "rna": "ACGU"}
        require(molecule_type in alphabet, "Unknown molecule_type")
        sequence = polymer["sequence"]
        require(isinstance(sequence, str) and 1 <= len(sequence) <= 4096
                and set(sequence) <= set(alphabet[molecule_type]), "Invalid polymer sequence; uppercase required")
        if "cyclic" in polymer:
            require(type(polymer["cyclic"]) is bool, "cyclic must be boolean")
        modifications = polymer.get("modifications") or []
        require(isinstance(modifications, list), "modifications must be a list")
        positions = set()
        for modification in modifications:
            _keys(modification, {"ccd", "position"}, {"ccd", "position"}, "modification")
            require(isinstance(modification["ccd"], str) and re.fullmatch(r"[A-Z0-9]{1,5}", modification["ccd"]) is not None, "Invalid modification CCD")
            require(_positive_index(modification["position"], len(sequence)), "Modification position outside sequence")
            require(modification["position"] not in positions, "Duplicate modification position")
            positions.add(modification["position"])
        msa = polymer.get("msa")
        if msa is not None:
            require(molecule_type == "protein" and isinstance(msa, dict) and len(msa) <= 3, "Invalid MSA mapping")
            for formats in msa.values():
                require(isinstance(formats, dict) and bool(formats), "Empty MSA format mapping")
                for fmt, alignment in formats.items():
                    require(fmt in {"a3m", "csv", "fasta", "sto"}, "Unknown MSA format")
                    _keys(alignment, {"format", "alignment"}, {"format", "alignment", "rank"}, "alignment")
                    require(alignment["format"] == fmt and isinstance(alignment["alignment"], str), "Alignment format/text mismatch")
                    if "rank" in alignment:
                        require(type(alignment["rank"]) is int, "MSA rank must be an integer")
                    if fmt == "a3m":
                        require(_a3m_query(alignment["alignment"]) == sequence, "MSA query differs from polymer sequence")
        templates = polymer.get("structural_templates") or []
        require(isinstance(templates, list) and len(templates) <= 4, "Invalid structural templates")
        require(not templates or molecule_type == "protein", "Templates require a protein")
        for template in templates:
            _keys(template, {"structure"}, {"structure", "format", "name", "chain_id"}, "template")
            require(isinstance(template["structure"], str) and bool(template["structure"]), "Empty template")
            require(template.get("format", "cif") in {"cif", "pdb"}, "Invalid template format")
        polymer_by_id[polymer["id"]] = polymer
    active_affinity = []
    for ligand in ligands:
        _keys(ligand, {"id"}, {"id", "ccd", "smiles", "predict_affinity", "output_affinity_embedding"}, "ligand")
        require(sum(ligand.get(k) is not None for k in ("ccd", "smiles")) == 1, "Ligand requires exactly one of CCD or SMILES")
        chemistry = ligand.get("ccd") or ligand.get("smiles")
        require(isinstance(chemistry, str) and bool(chemistry.strip()), "Empty ligand representation")
        if ligand.get("ccd") is not None:
            require(re.fullmatch(r"[A-Z0-9]{1,5}", chemistry) is not None, "Invalid ligand CCD")
        for key in ("predict_affinity", "output_affinity_embedding"):
            require(type(ligand.get(key, False)) is bool, key + " must be boolean")
        require(not ligand.get("output_affinity_embedding") or ligand.get("predict_affinity"), "Embedding requires affinity on the same ligand")
        if ligand.get("predict_affinity"):
            active_affinity.append(ligand["id"])
    require(len(active_affinity) <= 1, "At most one affinity ligand")
    if active_affinity:
        chemicals = [x.get("ccd") or x.get("smiles") for x in ligands]
        require(len(set(chemicals)) == len(chemicals), "Affinity requests cannot duplicate ligand representations")
    for entity in polymers + ligands:
        entity_id = entity["id"]
        require(isinstance(entity_id, str) and re.fullmatch(r"[A-Za-z0-9]{1,4}", entity_id) is not None, "Invalid entity ID")
        require(entity_id not in identities, "Duplicate entity ID")
        identities.add(entity_id)
    pocket_binders = set()
    for constraint in constraints:
        require(isinstance(constraint, dict), "Constraint must be an object")
        if constraint.get("constraint_type") == "bond":
            _keys(constraint, {"constraint_type", "atoms"}, {"constraint_type", "atoms"}, "bond")
            require(isinstance(constraint["atoms"], list) and len(constraint["atoms"]) == 2, "Bond requires two atoms")
            contacts = constraint["atoms"]
            fields = {"id", "residue_index", "atom_name"}
        elif constraint.get("constraint_type") == "pocket":
            _keys(constraint, {"constraint_type", "binder", "contacts"}, {"constraint_type", "binder", "contacts"}, "pocket")
            require(constraint["binder"] in {x["id"] for x in ligands}, "Unknown pocket ligand")
            pocket_binders.add(constraint["binder"])
            contacts, fields = constraint["contacts"], {"id", "residue_index"}
            require(isinstance(contacts, list) and bool(contacts), "Pocket requires contacts")
        else:
            raise ValueError("Unknown constraint type")
        for contact in contacts:
            _keys(contact, fields, fields, "contact")
            require(contact["id"] in identities, "Unknown constraint entity")
            size = len(polymer_by_id[contact["id"]]["sequence"]) if contact["id"] in polymer_by_id else 1
            require(_positive_index(contact["residue_index"], size), "Constraint residue outside entity")
            if "atom_name" in contact:
                require(isinstance(contact["atom_name"], str) and bool(contact["atom_name"]), "Invalid atom name")
            else:
                require(contact["id"] in polymer_by_id, "Pocket contact must be a polymer")
    require(len(pocket_binders) <= 1, "Only one pocket binder per request")
    parameters = options["parameters"]
    require(isinstance(parameters, dict), "parameters must be an object")
    allowed = set(INT_OPTIONS) | BOOL_OPTIONS | {"step_scale", "output_format"}
    require(set(parameters) <= allowed, "Unknown or unsupported hosted inference control")
    for key, value in parameters.items():
        if key in INT_OPTIONS:
            lo, hi = INT_OPTIONS[key]
            require(type(value) is int and lo <= value <= hi, "Invalid " + key)
        elif key in BOOL_OPTIONS:
            require(type(value) is bool, "Invalid " + key)
        elif key == "step_scale":
            require(type(value) in (int, float) and math.isfinite(value) and 0.5 <= value <= 5.0, "Invalid step_scale")
        elif key == "output_format":
            require(value == "mmcif", "Only mmCIF output is supported")
    # Serialization also rejects NaN/Infinity and detaches the caller's objects.
    return json.loads(canonical_json({**system, **parameters}))


def inspect_structure(text, request):
    require(isinstance(text, str) and bool(text), "Empty structure")
    data = MMCIF2Dict(io.StringIO(text))
    atom_ids = data.get("_atom_site.id", [])
    require(bool(atom_ids) and len(set(atom_ids)) == len(atom_ids), "Missing or duplicate atom IDs")
    count = len(atom_ids)
    for axis in "xyz":
        values = data.get("_atom_site.Cartn_" + axis, [])
        require(len(values) == count and all(math.isfinite(float(x)) for x in values), "Invalid coordinates")
    chains = data.get("_atom_site.label_asym_id", [])
    require(len(chains) == count, "Missing chain mapping")
    expected = {x["id"] for x in request["polymers"] + request["ligands"]}
    require(set(chains) == expected, "Output chain set differs from input entity IDs")
    positions = data.get("_atom_site.label_seq_id", [])
    residues = data.get("_atom_site.label_comp_id", [])
    require(len(positions) == len(residues) == count, "Missing residue mapping")
    for polymer in request["polymers"]:
        observed = {}
        for chain, position, residue in zip(chains, positions, residues):
            if chain != polymer["id"]:
                continue
            index = int(position)
            require(index not in observed or observed[index] == residue, "Conflicting residue identity")
            observed[index] = residue
        require(set(observed) == set(range(1, len(polymer["sequence"]) + 1)), "Missing or extra polymer residues")
        modified = {m["position"]: m["ccd"] for m in (polymer.get("modifications") or [])}
        for index, expected_letter in enumerate(polymer["sequence"], 1):
            residue = observed[index]
            if index in modified:
                require(residue == modified[index], "Modified residue mismatch")
            elif polymer["molecule_type"] == "protein":
                require(expected_letter == "X" or seq1(residue) == expected_letter, "Protein residue identity mismatch")
            else:
                native_letter = residue[1:] if polymer["molecule_type"] == "dna" and residue.startswith("D") else residue
                require(native_letter == expected_letter, "Nucleic-acid residue identity mismatch")
    return sorted(set(chains))


def parse_response(raw, request, observation):
    text = raw.decode("utf-8")
    response = json.loads(text)
    canonical_json(response)  # Reject non-finite provider JSON without rewriting bytes.
    require(isinstance(response, dict), "Response must be an object")
    structures = response.get("structures")
    require(isinstance(structures, list) and bool(structures), "No returned structures")
    require(len(structures) == request.get("diffusion_samples", 1), "Returned structure count differs from requested samples")
    request_hash = digest(canonical_json(request).encode())
    items = []
    for index, structure in enumerate(structures):
        require(isinstance(structure, dict) and structure.get("format") == "mmcif", "Unexpected structure format")
        mmcif = structure.get("structure")
        chain_ids = inspect_structure(mmcif, request)
        encoded = mmcif.encode("utf-8")
        items.append({"sample_index": index, "content": mmcif, "format": "mmcif",
                      "sha256": digest(encoded), "bytes": len(encoded), "chain_ids": chain_ids,
                      "name": structure.get("name"), "source": structure.get("source")})
    affinities = response.get("affinities")
    require(isinstance(affinities, dict), "Missing affinity object")
    targets = {x["id"] for x in request["ligands"] if x.get("predict_affinity")}
    require(set(affinities) == targets, "Affinity ligand identity mismatch")
    for ligand_id, scores in affinities.items():
        require(isinstance(scores, dict), "Invalid affinity entry")
        embedding_requested = next(x for x in request["ligands"] if x["id"] == ligand_id).get("output_affinity_embedding", False)
        if embedding_requested:
            require(all(scores.get(k) is not None for k in ("affinity_embedding", "model_1_affinity_embedding", "model_2_affinity_embedding")), "Requested affinity embeddings are missing")
        for key in ("affinity_pred_value", "affinity_probability_binary"):
            require(isinstance(scores.get(key), list) and bool(scores[key]), "Missing native affinity samples")
        for key, values in scores.items():
            if "probability" in key:
                require(isinstance(values, list) and all(type(v) in (int, float) and 0 <= v <= 1 for v in values), "Invalid native probability")
            elif "embedding" in key:
                if values is None and not embedding_requested:
                    continue
                require(isinstance(values, list) and bool(values) and all(isinstance(row, list) and len(row) == 384 and all(type(v) in (int, float) for v in row) for row in values), "Invalid affinity embedding")
            else:
                require(isinstance(values, list) and all(type(v) in (int, float) for v in values), "Invalid affinity values")
    quality = {key: response[key] for key in QUALITY_FIELDS if key in response}
    for key in QUALITY_FIELDS[:9]:
        require(key in quality and isinstance(quality[key], list) and len(quality[key]) == len(items)
                and all(type(v) in (int, float) for v in quality[key]), "Invalid per-structure quality field: " + key)
    units = {key: ("angstrom" if key in {"complex_pde_scores", "complex_ipde_scores", "pae", "pde"}
                   else "1; native scale preserved") for key in quality}
    return {
        "predicted_complexes": {"items": items, "coordinate_unit": "angstrom", "request_sha256": request_hash},
        "affinity_predictions": {"by_ligand": affinities, "field_definitions": {
            "affinity_pred_value": "Native model affinity-like log-micromolar score; not measured affinity",
            "affinity_probability_binary": "Native model probability in [0,1]; not empirical accuracy",
            "affinity_pic50": "Separate provider field; conversion/calibration unverified",
            "embeddings": "Optional native affinity vectors, width 384; unit 1; no cross-model compatibility claimed",
        }, "request_sha256": request_hash},
        "quality_metrics": {"values": quality, "field_units": units, "request_sha256": request_hash},
        "run_provenance": {"request": request, "request_sha256": request_hash,
            "response_sha256": digest(raw), "response_bytes": len(raw), **observation,
            "backend_identity": {"checkpoint_sha256": None, "container_digest": None, "seed": None},
            "limitations": ["Scientific qualification incomplete", "Exact hosted backend identity unavailable",
                "Confidence is not measured accuracy", "Ligand connectivity/stereochemistry not checked",
                "Non-A3M alignment query checks delegated to provider", "Full PAE/PDE artifact retrieval not implemented",
                "Native numeric chain-score indices are preserved; mapping to input IDs requires provider verification"]},
        "native_response": text,
    }


class NvidiaBoltz2(BioModule):
    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN

    def __init__(self, observation_timeout_s=1200.0, max_response_bytes=20_000_000,
                 resume_request_id=None, resume_request_sha256=None):
        require(type(observation_timeout_s) in (int, float) and math.isfinite(observation_timeout_s)
                and 1 <= observation_timeout_s <= 3600, "Invalid observation timeout")
        require(type(max_response_bytes) is int and 1024 <= max_response_bytes <= 100_000_000, "Invalid response limit")
        self.observation_timeout_s = float(observation_timeout_s)
        self.max_response_bytes = max_response_bytes
        self.resume_request_id = str(uuid.UUID(resume_request_id)) if resume_request_id else None
        require(bool(resume_request_id) == bool(resume_request_sha256), "Resume requires the original request ID and SHA-256 together")
        if resume_request_sha256:
            require(isinstance(resume_request_sha256, str) and re.fullmatch(r"[0-9a-f]{64}", resume_request_sha256) is not None, "Invalid resume request digest")
        self.resume_request_sha256 = resume_request_sha256
        self._observation = {"status": "not_submitted", "provider_request_id": None, "request_sha256": None}
        self._last_raw_response = None

    def snapshot(self):
        return dict(self._observation)

    def inputs(self):
        return {
            "molecular_system": SignalSpec.record(schema=SYSTEM_SCHEMA, accepted_units=("1",), required=True,
                description="Polymers, ligands and constraints with explicit entity IDs."),
            "run_options": SignalSpec.record(schema=OPTION_SCHEMA, accepted_units=("1",), required=False,
                default={"parameters": {}}, description="Declared native inference controls; no credentials."),
        }

    def outputs(self):
        return {
            "predicted_complexes": SignalSpec.record(schema=STRUCTURE_SCHEMA, emitted_unit="1"),
            "affinity_predictions": SignalSpec.record(schema=AFFINITY_SCHEMA, emitted_unit="1"),
            "quality_metrics": SignalSpec.record(schema=QUALITY_SCHEMA, emitted_unit="1"),
            "run_provenance": SignalSpec.record(schema=PROVENANCE_SCHEMA, emitted_unit="1"),
            "native_response": SignalSpec.scalar(dtype="str", format="json", emitted_unit="1"),
        }

    def execute(self, inputs, *, context):
        require("molecular_system" in inputs, "molecular_system is required")
        options = inputs["run_options"].value if "run_options" in inputs else {"parameters": {}}
        request = prepare_request(inputs["molecular_system"].value, options)
        request_hash = digest(canonical_json(request).encode())
        require(not self.resume_request_id or request_hash == self.resume_request_sha256, "Resume inputs differ from original request digest")
        key = os.environ.get("NVIDIA_API_KEY", "").strip()
        if not key:
            raise ProviderObservationError("NVIDIA_API_KEY is unavailable in this execution environment")
        started = time.monotonic()
        request_id = self.resume_request_id
        self._observation = {"status": "observing" if request_id else "submitting",
                             "provider_request_id": request_id, "request_sha256": request_hash}
        headers = {"Authorization": "Bearer " + key, "Accept": "application/json", "Content-Type": "application/json"}
        method = "GET" if request_id else "POST"
        url = POLL_ENDPOINT + request_id if request_id else ENDPOINT
        body = None if request_id else canonical_json(request).encode()
        safe_headers = {}
        with httpx.Client(follow_redirects=False, trust_env=False) as client:
            while True:
                remaining = self.observation_timeout_s - (time.monotonic() - started)
                if remaining <= 0:
                    raise ProviderObservationError("Observation timed out; provider may still run. request_id=" + (request_id or "unknown") + "; request_sha256=" + request_hash)
                call_headers = {**headers, **({"NVCF-POLL-SECONDS": str(min(30, max(1, int(remaining))))} if method == "GET" else {})}
                try:
                    with client.stream(method, url, content=body, headers=call_headers, timeout=min(60.0, remaining)) as response:
                        returned_id = response.headers.get("nvcf-reqid")
                        if returned_id:
                            parsed_id = str(uuid.UUID(returned_id))
                            require(not request_id or request_id == parsed_id, "Provider request ID changed")
                            request_id = parsed_id
                            self._observation["provider_request_id"] = request_id
                        safe_headers.update({k: v for k, v in response.headers.items() if k in {"nvcf-reqid", "nvcf-status", "nvcf-function-id", "nvcf-function-version-id", "date", "content-type"}})
                        chunks, size = [], 0
                        for chunk in response.iter_bytes():
                            size += len(chunk)
                            if size > self.max_response_bytes:
                                raise ProviderObservationError("Response size limit exceeded; request_id=" + (request_id or "unknown"))
                            if time.monotonic() - started > self.observation_timeout_s:
                                raise ProviderObservationError("Observation timed out during response; request_id=" + (request_id or "unknown"))
                            chunks.append(chunk)
                        raw = b"".join(chunks)
                        status = response.status_code
                        self._observation["http_status"] = status
                        self._observation["status"] = "provider_pending" if status == 202 else "response_received"
                except httpx.HTTPError:
                    raise ProviderObservationError("Transport observation failed; no automatic resubmission. request_id=" + (request_id or "unknown") + "; request_sha256=" + request_hash) from None
                if status == 200:
                    break
                if status != 202:
                    raise ProviderObservationError("Provider HTTP " + str(status) + "; request_id=" + (request_id or "unknown"))
                if not request_id:
                    raise ProviderObservationError("Accepted request has no observable ID; do not resubmit automatically")
                method, url, body = "GET", POLL_ENDPOINT + request_id, None
                time.sleep(min(2.0, max(0.0, self.observation_timeout_s - (time.monotonic() - started))))
        observation = {"provider_request_id": request_id or "", "endpoint": ENDPOINT,
            "captured_at": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.monotonic() - started, "resumed": bool(self.resume_request_id),
            "provider_headers": safe_headers}
        if key.encode() in raw:
            raise ProviderObservationError("Provider response contains credential material; capture rejected")
        self._last_raw_response = raw
        result = parse_response(raw, request, observation)
        self._observation["status"] = "outputs_validated"
        return result
