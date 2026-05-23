from __future__ import annotations

from pathlib import Path

import yaml


EXPECTED_ALIASES = ["malaria_target_context", "falcipain_ligand_setup", "boltz_boltz2_affinity_predictor", "malaria_prediction_interpreter", "visualisation"]


def _manifest():
    lab_dir = Path(__file__).resolve().parents[1]
    return yaml.safe_load((lab_dir / "lab.yaml").read_text())


def test_lab_manifest_uses_multi_stage_workflow_graph():
    manifest = _manifest()
    model_aliases = [entry["alias"] for entry in manifest["models"]]
    assert model_aliases == EXPECTED_ALIASES
    assert all(str(entry["path"]).startswith("models/") for entry in manifest["models"])
    assert all(not str(entry["maps_to"]).startswith("visualisation.") for entry in manifest["io"]["outputs"])

    wiring = manifest["wiring"]
    wiring_sources = {entry["from"] for entry in wiring}
    wiring_targets = {target for entry in wiring for target in entry["to"]}
    assert "malaria_target_context.scenario_context" in wiring_sources
    assert "falcipain_ligand_setup.assembled_boltz_request" in wiring_sources
    assert "malaria_prediction_interpreter.prediction_evidence" in wiring_sources
    assert any(target.startswith("visualisation.") for target in wiring_targets)


def test_lab_manifest_keeps_boltz_prediction_as_only_scientific_runner():
    manifest = _manifest()
    assert manifest["title"].startswith("Boltz Workflow:")
    assert "boltz-workflow" in manifest["tags"]
    assert "guided-workflow" in manifest["tags"]

    context_model = next(entry for entry in manifest["models"] if entry["alias"] == "malaria_target_context")
    assembler_model = next(entry for entry in manifest["models"] if entry["alias"] == "falcipain_ligand_setup")
    core_model = next(entry for entry in manifest["models"] if entry["alias"] == "boltz_boltz2_affinity_predictor")
    interpreter_model = next(entry for entry in manifest["models"] if entry["alias"] == "malaria_prediction_interpreter")

    assert context_model["parameters"]["scenario"]["caveat"]
    assert assembler_model["parameters"]["default_protein_sequence"]
    assert assembler_model["parameters"]["default_ligand_smiles"]
    assert assembler_model["parameters"]["default_run_options"]["workflow_name"]
    assert "default_protein_sequence" not in core_model["parameters"]
    assert "default_ligand_smiles" not in core_model["parameters"]
    assert "default_ligand_csv" not in core_model["parameters"]
    assert core_model["parameters"]["accelerator"] == "gpu"
    assert interpreter_model["parameters"]["core_alias"] == "boltz_boltz2_affinity_predictor"


def test_public_ports_route_through_workflow_stages():
    manifest = _manifest()
    input_maps = [entry["maps_to"] for entry in manifest["io"]["inputs"]]
    assert all(item.startswith("falcipain_ligand_setup.") for item in input_maps)
    output_maps = [entry["maps_to"] for entry in manifest["io"]["outputs"]]
    assert "malaria_target_context.scenario_context" in output_maps
    assert "falcipain_ligand_setup.assembled_boltz_request" in output_maps
    assert "malaria_prediction_interpreter.prediction_evidence" in output_maps
