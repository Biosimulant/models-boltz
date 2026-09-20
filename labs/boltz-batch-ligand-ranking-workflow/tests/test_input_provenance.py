"""Input-routing and attribution checks; no biological prediction is performed."""
import hashlib
import importlib.util
from pathlib import Path

import pytest
import yaml
from biosim import ExecutionContext, ExecutionPolicy
from biosim.signals import make_signal, unwrap_payload

LAB = Path(__file__).resolve().parents[1]


def load_stage(directory, filename, classname, **kwargs):
    path = LAB / "models" / directory / "src" / (filename + ".py")
    spec = importlib.util.spec_from_file_location("provenance_" + filename, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, classname)(**kwargs)


def signal(name, value):
    return make_signal(source="test", name=name, value=value, emitted_at=0, spec=None)


def execute(stage, inputs):
    return {k: unwrap_payload(v) for k, v in stage.execute(
        {k: signal(k, v) for k, v in inputs.items()},
        context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0, run_end=1),
    ).items()}


@pytest.mark.parametrize("protein_changed,ligand_changed", [(False, False), (True, False), (False, True), (True, True)])
def test_resolved_context_reaches_interpreter_and_visualisation(protein_changed, ligand_changed):
    manifest = yaml.safe_load((LAB / "lab.yaml").read_text())
    entries = {x["alias"]: x for x in manifest["models"]}
    params = entries["ligand_library_loader"]["parameters"]
    context = entries["batch_target_context"]["parameters"]["scenario"]
    assembler = load_stage("input_assembler", "boltz_input_assembler", "BoltzInputAssemblerModel", **params)
    inputs = {"scenario_context": context}
    if protein_changed:
        # Deliberately truncated routing fixture; no inference or validity claim.
        inputs["protein_sequence"] = params["default_protein_sequence"][:-1]
    if ligand_changed:
        inputs["ligand_csv"] = "\n".join(params["default_ligand_csv"].splitlines()[:2])
    result = execute(assembler, inputs)
    request = result["assembled_boltz_request"]
    resolved = request["effective_context"]
    assert request["input_provenance"]["protein_matches_packaged_example"] is not protein_changed
    assert request["input_provenance"]["ligands_match_packaged_example"] is not ligand_changed
    assert request["protein_sequence_sha256"] == hashlib.sha256(result["protein_sequence"].encode()).hexdigest()
    if protein_changed:
        assert resolved["target_name"] == "User-supplied protein"
        assert "source_pdb" not in resolved
        assert "source_pdb" not in result["run_options"]
        assert resolved["protein_sequence_length"] == 272
    else:
        assert resolved["source_pdb"] == "2HYY"
    if ligand_changed:
        assert resolved["ligand_name"] == "User-supplied ligand library"
        assert "ligand_examples" not in result["run_options"]
    interpreter = load_stage("interpreter", "boltz_prediction_interpreter", "BoltzPredictionInterpreterModel", **entries["ranking_interpreter"]["parameters"])
    evidence = execute(interpreter, {"scenario_context": context, "assembled_boltz_request": request})["prediction_evidence"]
    assert evidence["target_name"] == request["target_name"]
    assert evidence["ligand_name"] == request["ligand_name"]
    assert evidence["scientific_question"] == resolved["workflow_question"]
    visual = load_stage("visualisation", "docking_visualisation", "DockingVisualisationModel", **entries["visualisation"]["parameters"])
    execute(visual, {
        "batch_target_context_scenario_context": context,
        "ligand_library_loader_assembled_boltz_request": request,
        "ranking_interpreter_prediction_evidence": evidence,
    })
    table = next(x["data"] for x in visual._workflow_visuals() if x["data"]["title"] == "Workflow target and ligand context")
    rows = dict(table["rows"])
    assert rows["Target"] == request["target_name"]
    assert ("Source PDB" not in rows) is protein_changed


def test_blank_inputs_do_not_silently_select_packaged_example():
    stage = load_stage("input_assembler", "boltz_input_assembler", "BoltzInputAssemblerModel", workflow_name="test", workflow_kind="batch", default_protein_sequence="ACD", default_ligand_csv="name,smiles\nfixture,C", default_msa_path="example.a3m")
    result = execute(stage, {"protein_sequence": " ", "ligand_csv": " ", "msa_path": " "})
    assert result["protein_sequence"] == result["ligand_csv"] == result["msa_path"] == ""
    assert result["assembled_boltz_request"]["input_provenance"]["protein_matches_packaged_example"] is False


def test_long_csv_text_is_not_treated_as_a_filename(tmp_path):
    stage = load_stage("input_assembler", "boltz_input_assembler", "BoltzInputAssemblerModel", workflow_name="test", workflow_kind="batch")
    csv_text = "name,smiles\n" + "\n".join(f"fixture{i},C" for i in range(100))
    assert execute(stage, {"ligand_csv": csv_text})["ligand_csv"] == csv_text
    csv_path = tmp_path / "ligands.csv"
    csv_path.write_text(csv_text)
    assert execute(stage, {"ligand_csv": str(csv_path)})["ligand_csv"] == csv_text
