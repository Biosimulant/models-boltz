"""Technical tests only. Synthetic transport payloads are not inference evidence."""
import copy
import importlib.util
import json
from pathlib import Path

import httpx
import pytest
import yaml
from biosim import BioModule, BioWorld, ExecutionContext, ExecutionPolicy, SignalSpec, make_signal
from biosim.pack import validate_lab_source

LAB = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("nvidia_boltz2_model", LAB / "owned/models/main/model.py")
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)
CONTEXT = ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0, run_end=1)
REQUEST_ID = "12345678-1234-5678-1234-567812345678"


@pytest.fixture
def molecular_system():
    return {"polymers": [{"id": "A", "molecule_type": "protein", "sequence": "A"}], "ligands": [], "constraints": []}


@pytest.fixture
def native_response():
    # Deliberately minimal synthetic coordinate fragment for parser failure tests.
    cif = """data_technical_fixture
loop_
_atom_site.id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 A 1 ALA 0.0 1.0 2.0
#
"""
    return {"structures": [{"structure": cif, "format": "mmcif", "name": "synthetic", "source": "technical-test"}],
            "affinities": {}, **{key: [0.5] for key in model.QUALITY_FIELDS[:9]},
            "chains_ptm_scores": [0.5], "pair_chains_iptm_scores": [{"0": {"0": 0.5}}], "pae": None, "pde": None}


def signals(system, options=None):
    module = model.NvidiaBoltz2()
    values = {"molecular_system": system, "run_options": {"parameters": options or {}}}
    return {name: make_signal(spec=module.inputs()[name], value=value, source="fixture", name=name, emitted_at=0.0) for name, value in values.items()}


def mock_client(monkeypatch, handler):
    client_class = httpx.Client
    monkeypatch.setattr(model.httpx, "Client", lambda **kwargs: client_class(transport=httpx.MockTransport(handler), **kwargs))
    monkeypatch.setenv("NVIDIA_API_KEY", "test-credential-not-valid")
    monkeypatch.setattr(model.time, "sleep", lambda _: None)


def observation():
    return {"provider_request_id": REQUEST_ID, "endpoint": model.ENDPOINT, "captured_at": "2026-09-20T00:00:00+00:00",
            "wall_seconds": 0.0, "resumed": False, "provider_headers": {}}


def test_manifest_runtime_parity():
    module = model.NvidiaBoltz2()
    manifest = yaml.safe_load((LAB / "owned/models/main/model.yaml").read_text())
    for direction, actual in [("inputs", module.inputs()), ("outputs", module.outputs())]:
        declared = {port["name"]: port for port in manifest["io"][direction]}
        assert set(actual) == set(declared)
        for name, runtime_spec in actual.items():
            port = declared[name]
            assert port["signal_type"] == runtime_spec.signal_type
            assert port.get("schema") == runtime_spec.schema
            assert port.get("dtype") == runtime_spec.dtype
            assert port.get("format") == runtime_spec.format
            assert port.get("emitted_unit") == runtime_spec.emitted_unit
            assert tuple(port.get("accepted_units", ())) == (runtime_spec.accepted_units or ())
            assert port.get("required") == runtime_spec.required
            assert port.get("default") == runtime_spec.default
    lab = yaml.safe_load((LAB / "lab.yaml").read_text())
    assert len(lab["models"]) == 1
    for direction, ports in [("inputs", module.inputs()), ("outputs", module.outputs())]:
        assert {p["maps_to"] for p in lab["io"][direction]} == {"main." + p for p in ports}
    result = validate_lab_source(LAB)
    assert result.valid, result.errors
    assert module.execution_policy == ExecutionPolicy.ONCE_BEFORE_RUN


@pytest.mark.parametrize("case", ["duplicate_id", "bad_sequence", "wrong_msa", "bad_msa_columns", "missing_id", "bad_modification", "nan_option", "seed_option", "bool_integer", "two_affinities", "wrong_bond", "empty_ligand"])
def test_rejects_invalid_request(molecular_system, case):
    options = {"parameters": {}}
    polymer = molecular_system["polymers"][0]
    if case == "duplicate_id": molecular_system["polymers"].append(copy.deepcopy(polymer))
    elif case == "bad_sequence": polymer["sequence"] = "A*"
    elif case in {"wrong_msa", "bad_msa_columns"}:
        polymer["msa"] = {"supplied": {"a3m": {"format": "a3m", "alignment": ">q\nC\n" if case == "wrong_msa" else ">q\nA\n>x\nAA\n"}}}
    elif case == "missing_id": del polymer["id"]
    elif case == "bad_modification": polymer["modifications"] = [{"ccd": "SEP", "position": 2}]
    elif case == "nan_option": options["parameters"]["step_scale"] = float("nan")
    elif case == "seed_option": options["parameters"]["seed"] = 1
    elif case == "bool_integer": options["parameters"]["sampling_steps"] = True
    elif case == "two_affinities": molecular_system["ligands"] = [{"id": "B", "smiles": "C", "predict_affinity": True}, {"id": "C", "smiles": "CC", "predict_affinity": True}]
    elif case == "wrong_bond": molecular_system["constraints"] = [{"constraint_type": "bond", "atoms": [{"id": "Z", "residue_index": 1, "atom_name": "CA"}, {"id": "A", "residue_index": 1, "atom_name": "CA"}]}]
    elif case == "empty_ligand": molecular_system["ligands"] = [{"id": "B", "smiles": ""}]
    with pytest.raises(ValueError): model.prepare_request(molecular_system, options)


def test_request_preserves_dna_rna_modifications_and_options():
    system = {"polymers": [{"id": "D", "molecule_type": "dna", "sequence": "ACGT"},
                           {"id": "R", "molecule_type": "rna", "sequence": "ACGU"},
                           {"id": "P", "molecule_type": "protein", "sequence": "AS", "modifications": [{"ccd": "SEP", "position": 2}]}],
              "ligands": [], "constraints": []}
    assert model.prepare_request(system, {"parameters": {"diffusion_samples": 2}}) == {**system, "diffusion_samples": 2}


def test_missing_secret_fails_before_network(monkeypatch, molecular_system):
    monkeypatch.delenv("NVIDIA_API_KEY", raising=False)
    monkeypatch.setattr(model.httpx, "Client", lambda **_: pytest.fail("Network attempted without key"))
    with pytest.raises(model.ProviderObservationError, match="unavailable"):
        model.NvidiaBoltz2().execute(signals(molecular_system), context=CONTEXT)


def test_single_submission_then_poll_and_portable_outputs(monkeypatch, molecular_system, native_response):
    calls = []
    raw = json.dumps(native_response, indent=2).encode()
    def handler(request):
        calls.append(request)
        assert request.headers["authorization"] == "Bearer test-credential-not-valid"
        if len(calls) == 1: return httpx.Response(202, headers={"nvcf-reqid": REQUEST_ID})
        return httpx.Response(200, content=raw, headers={"nvcf-reqid": REQUEST_ID})
    mock_client(monkeypatch, handler)
    result = model.NvidiaBoltz2().execute(signals(molecular_system), context=CONTEXT)
    assert [c.method for c in calls] == ["POST", "GET"]
    assert str(calls[1].url) == model.POLL_ENDPOINT + REQUEST_ID
    assert json.loads(calls[0].content) == molecular_system
    assert result["native_response"].encode() == raw
    assert result["run_provenance"]["response_sha256"] == model.digest(raw)
    assert result["predicted_complexes"]["items"][0]["sha256"] == model.digest(native_response["structures"][0]["structure"].encode())
    for name, port in model.NvidiaBoltz2().outputs().items():
        assert make_signal(spec=port, value=result[name], source="main", name=name, emitted_at=0.0).value == result[name]
    assert "test-credential-not-valid" not in json.dumps(result)


@pytest.mark.parametrize("status", [302, 401, 422, 429, 500])
def test_http_error_no_retry_or_credential_redirect(monkeypatch, molecular_system, status):
    calls = []
    def handler(request):
        calls.append(request)
        return httpx.Response(status, headers={"location": "https://untrusted.invalid/steal"}, content=b"test-credential-not-valid")
    mock_client(monkeypatch, handler)
    with pytest.raises(model.ProviderObservationError) as error:
        model.NvidiaBoltz2().execute(signals(molecular_system), context=CONTEXT)
    assert len(calls) == 1
    assert "test-credential-not-valid" not in str(error.value)
    assert "untrusted" not in str(error.value)


def test_transport_failure_does_not_resubmit(monkeypatch, molecular_system):
    calls = []
    def handler(request):
        calls.append(request)
        raise httpx.ReadTimeout("secret: test-credential-not-valid")
    mock_client(monkeypatch, handler)
    with pytest.raises(model.ProviderObservationError) as error:
        model.NvidiaBoltz2().execute(signals(molecular_system), context=CONTEXT)
    assert len(calls) == 1
    assert "test-credential-not-valid" not in str(error.value)
    assert "request_sha256=" in str(error.value)


def test_resume_polls_only_same_request(monkeypatch, molecular_system, native_response):
    calls = []
    def handler(request):
        calls.append(request)
        return httpx.Response(200, json=native_response, headers={"nvcf-reqid": REQUEST_ID})
    mock_client(monkeypatch, handler)
    sha = model.digest(model.canonical_json(molecular_system).encode())
    module = model.NvidiaBoltz2(resume_request_id=REQUEST_ID, resume_request_sha256=sha)
    assert module.execute(signals(molecular_system), context=CONTEXT)["run_provenance"]["resumed"]
    assert [x.method for x in calls] == ["GET"]
    changed = copy.deepcopy(molecular_system)
    changed["polymers"][0]["sequence"] = "C"
    with pytest.raises(ValueError, match="Resume inputs differ"):
        module.execute(signals(changed), context=CONTEXT)
    assert len(calls) == 1


@pytest.mark.parametrize("case", ["nan", "lost_chain", "wrong_residue", "duplicate_atom", "no_structures", "sample_count", "probability", "score_shape"])
def test_output_corruption_rejected(molecular_system, native_response, case):
    if case == "nan": native_response["confidence_scores"][0] = float("nan")
    elif case == "lost_chain": native_response["structures"][0]["structure"] = native_response["structures"][0]["structure"].replace("1 A 1", "1 B 1")
    elif case == "wrong_residue": native_response["structures"][0]["structure"] = native_response["structures"][0]["structure"].replace("ALA", "CYS")
    elif case == "duplicate_atom": native_response["structures"][0]["structure"] += "1 A 1 ALA 3.0 4.0 5.0\n"
    elif case == "no_structures": native_response["structures"] = []
    elif case == "sample_count": molecular_system["diffusion_samples"] = 2
    elif case == "probability":
        molecular_system["ligands"] = [{"id": "B", "smiles": "C", "predict_affinity": True}]
        native_response["structures"][0]["structure"] += "2 B . LIG 3.0 4.0 5.0\n"
        native_response["affinities"] = {"B": {"affinity_pred_value": [1.0], "affinity_probability_binary": [1.5]}}
    elif case == "score_shape": native_response["ptm_scores"] = []
    with pytest.raises(ValueError): model.parse_response(json.dumps(native_response).encode(), molecular_system, observation())


def test_size_limit(monkeypatch, molecular_system):
    mock_client(monkeypatch, lambda _: httpx.Response(200, content=b"x" * 1025))
    with pytest.raises(model.ProviderObservationError, match="size limit"):
        model.NvidiaBoltz2(max_response_bytes=1024).execute(signals(molecular_system), context=CONTEXT)


def test_null_embedding_preserved_but_requested_embedding_required(molecular_system, native_response):
    molecular_system["ligands"] = [{"id": "B", "smiles": "C", "predict_affinity": True}]
    native_response["structures"][0]["structure"] += "2 B . LIG 3.0 4.0 5.0\n"
    native_response["affinities"] = {"B": {"affinity_pred_value": [2.59375], "affinity_probability_binary": [0.4453125],
        "affinity_pic50": [4.646125], "affinity_embedding": None, "model_1_affinity_embedding": None, "model_2_affinity_embedding": None}}
    result = model.parse_response(json.dumps(native_response).encode(), molecular_system, observation())
    assert result["affinity_predictions"]["by_ligand"] == native_response["affinities"]
    molecular_system["ligands"][0]["output_affinity_embedding"] = True
    with pytest.raises(ValueError, match="embeddings are missing"):
        model.parse_response(json.dumps(native_response).encode(), molecular_system, observation())


def test_bioworld_executes_once_and_commits_all_ports(monkeypatch, molecular_system, native_response):
    calls = []
    def handler(request):
        calls.append(request)
        return httpx.Response(200, json=native_response, headers={"nvcf-reqid": REQUEST_ID})
    mock_client(monkeypatch, handler)
    class InputFixture(BioModule):
        execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN
        def outputs(self):
            return {"system": SignalSpec.record(schema=model.SYSTEM_SCHEMA, emitted_unit="1")}
        def execute(self, inputs, *, context):
            return {"system": molecular_system}
    world = BioWorld(communication_step=1.0)
    module = model.NvidiaBoltz2()
    world.add_biomodule("fixture", InputFixture())
    world.add_biomodule("main", module)
    world.connect("fixture.system", "main.molecular_system")
    world.run(duration=3.0)
    assert len(calls) == 1
    assert set(module.get_outputs()) == set(module.outputs())
    assert all(s.emitted_at == 0.0 for s in module.get_outputs().values())
    assert module.get_outputs()["predicted_complexes"].value["items"][0]["chain_ids"] == ["A"]
