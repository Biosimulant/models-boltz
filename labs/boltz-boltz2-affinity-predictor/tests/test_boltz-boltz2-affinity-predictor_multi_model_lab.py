from __future__ import annotations

from pathlib import Path

import yaml


def test_lab_manifest_uses_embedded_visualisation_model():
    lab_dir = Path(__file__).resolve().parents[1]
    manifest = yaml.safe_load((lab_dir / "lab.yaml").read_text())
    model_aliases = [entry["alias"] for entry in manifest["models"]]
    assert "visualisation" in model_aliases
    assert all(str(entry["path"]).startswith("models/") for entry in manifest["models"])

    assert not any(str(entry["maps_to"]).startswith("visualisation.") for entry in manifest["io"]["outputs"])

    wiring_targets = {target for entry in manifest["wiring"] for target in entry["to"]}
    assert any(target.startswith("visualisation.") for target in wiring_targets)
    assert manifest["version"] == "1.1.1"
    core = next(entry for entry in manifest["models"] if entry["alias"] == "boltz_boltz2_affinity_predictor")
    assert core["parameters"]["use_msa_server"] is False
    assert core["parameters"]["default_msa_path"] == "assets/seq1.a3m"
    assert core["parameters"]["output_format"] == "mmcif"
    inputs = {entry["name"]: entry for entry in manifest["io"]["inputs"]}
    assert inputs["msa_path"]["file"]["accept"] == [".a3m"]
    outputs = {entry["name"]: entry["maps_to"] for entry in manifest["io"]["outputs"]}
    assert outputs["binding_probability"].endswith(".binding_probability")
    assert outputs["affinity_log10_ic50_micromolar"].endswith(
        ".affinity_log10_ic50_micromolar"
    )
    assert outputs["predicted_structure"].endswith(".predicted_structure")
