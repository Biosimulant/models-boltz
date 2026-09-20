"""Real package/CLI plumbing against an explicit synthetic Boltz executable.

This fixture is never scientific inference evidence and allocates no GPU.
"""

import json
import shutil
import sys
from pathlib import Path

import pytest
import yaml
from biosim.__main__ import main
from biosim.pack import validate_lab_source

LAB = Path(__file__).parents[1]


def synthetic_cif(sequence, smiles):
    from Bio.SeqUtils import seq3
    from rdkit import Chem

    columns = [
        "label_asym_id",
        "label_seq_id",
        "label_comp_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "id",
        "pdbx_PDB_model_num",
    ]
    lines = [
        "data_explicit_synthetic_fixture",
        "loop_",
        *["_atom_site." + c for c in columns],
    ]
    for index, letter in enumerate(sequence, 1):
        lines.append(f"A {index} {seq3(letter).upper()} C {index}.0 0.0 0.0 {index} 1")
    for index, atom in enumerate(
        Chem.MolFromSmiles(smiles).GetAtoms(), len(sequence) + 1
    ):
        lines.append(f"B . LIG {atom.GetSymbol()} 0.0 {index}.0 0.0 {index} 1")
    return "\n".join(lines) + "\n"


def prepare(tmp_path, *, omit_affinity=False):
    lab = tmp_path / "lab"
    shutil.copytree(LAB, lab, ignore=shutil.ignore_patterns("__pycache__", ".runtime"))
    executable = tmp_path / "fixture_boltz"
    affinity = (
        {}
        if omit_affinity
        else {"affinity_probability_binary": 0.8, "affinity_pred_value": -1.2}
    )
    manifest = yaml.safe_load((lab / "lab.yaml").read_text())
    parameters = manifest["models"][0]["parameters"]
    cif = synthetic_cif(
        parameters["default_protein_sequence"], parameters["default_ligand_smiles"]
    )
    executable.write_text(f"""#!{sys.executable}
import json, pathlib, sys
out=pathlib.Path(sys.argv[sys.argv.index('--out_dir')+1])/'predictions'/'request'
out.mkdir(parents=True,exist_ok=True)
(out/'request_model_0.cif').write_text({cif!r})
(out/'confidence_request_model_0.json').write_text(json.dumps({{'confidence_score':.9}}))
(out/'affinity_request.json').write_text(json.dumps({affinity!r}))
""")
    executable.chmod(0o755)
    manifest = yaml.safe_load((lab / "lab.yaml").read_text())
    manifest["runtime"].pop("python_version", None)
    parameters = manifest["models"][0]["parameters"]
    parameters.update(
        runtime_mode="external",
        boltz_executable=str(executable),
        work_dir=str(tmp_path / "work"),
        cache_dir=str(tmp_path / "cache"),
        progress_heartbeat_s=0,
    )
    (lab / "lab.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    return lab


def test_profiled_boltz_pipeline_package_and_cli(tmp_path, capsys):
    lab = prepare(tmp_path)
    assert validate_lab_source(lab).valid
    result = tmp_path / "result.json"
    main(
        [
            "labs",
            "run",
            str(lab),
            "--no-install-deps",
            "--require-acceptance",
            "--results-file",
            str(result),
            "--json",
        ]
    )
    capsys.readouterr()
    report = json.loads(result.read_text())
    assert report["compatibility"]["summary"]["verified"] == 5
    assert report["compatibility"]["summary"]["partial"] == 0
    assert report["acceptance"]["status"] == "passed"
    assert report["acceptance"]["summary"]["passed"] == 17
    assert report["outputs"]["analysis"]["completed"]["value"] is True


def test_wrong_scientific_consumer_rejected_before_execution(tmp_path):
    lab = prepare(tmp_path)
    path = lab / "models/analysis/model.yaml"
    manifest = yaml.safe_load(path.read_text())
    manifest["io"]["inputs"][0]["contract"]["profile"] = (
        "boltz.log10-ic50-micromolar/v1"
    )
    path.write_text(yaml.safe_dump(manifest))
    assert not validate_lab_source(lab).valid
    assert not (tmp_path / "work").exists()


def test_missing_prediction_is_not_accepted(tmp_path, capsys):
    lab = prepare(tmp_path, omit_affinity=True)
    result = tmp_path / "result.json"
    # Depending on graph scheduling, absent required inputs can fail execution
    # before acceptance is evaluated. Neither outcome is an acceptance pass.
    try:
        main(
            [
                "labs",
                "run",
                str(lab),
                "--no-install-deps",
                "--require-acceptance",
                "--results-file",
                str(result),
                "--json",
            ]
        )
    except (SystemExit, ValueError, RuntimeError, KeyError):
        pass
    assert (
        not result.exists()
        or json.loads(result.read_text())["acceptance"]["status"] != "passed"
    )


def load_analyzer():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "boltz_acceptance_analysis", LAB / "models/analysis/src/analysis.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_structure_content_checks_reject_empty_wrong_identity_and_nonfinite(tmp_path):
    analyzer = load_analyzer()
    path = tmp_path / "structure.cif"
    good = synthetic_cif("MA", "CCO")
    path.write_text(good)
    assert analyzer.inspect_structure(str(path), "MA", "CCO")["protein_residues"] == 2
    for text, sequence, smiles in [
        ("data_empty\n", "MA", "CCO"),
        (good, "MG", "CCO"),
        (good, "MA", "CCC"),
        (good.replace("1.0 0.0 0.0", "nan 0.0 0.0"), "MA", "CCO"),
        (good.replace("B . LIG", "C . LIG"), "MA", "CCO"),
    ]:
        path.write_text(text)
        with pytest.raises(ValueError):
            analyzer.inspect_structure(str(path), sequence, smiles)


def test_declared_chain_mapping_and_multimer(tmp_path):
    import importlib.util

    module_spec = importlib.util.spec_from_file_location(
        "chain_analysis", LAB / "models/analysis/src/analysis.py"
    )
    analysis = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(analysis)
    original = synthetic_cif("MA", "CO")
    mapped = original.replace("\nA ", "\nX ").replace("\nB ", "\nL ")
    path = tmp_path / "mapped.cif"
    path.write_text(mapped)
    result = analysis.inspect_structure(
        str(path), protein_chains={"X": "MA"}, ligand_chains={"L": "CO"}
    )
    assert result["protein_chains"] == ["X"]
    with pytest.raises(ValueError):
        analysis.inspect_structure(str(path), "MA", "CO")
    # Add a second protein and ligand with unique atom IDs.
    lines = [
        line.split() for line in mapped.splitlines() if line.startswith(("X ", "L "))
    ]
    extras = []
    for row in lines:
        row[0] = {"X": "Y", "L": "J"}[row[0]]
        row[-2] = str(int(row[-2]) + 100)
        extras.append(" ".join(row))
    path.write_text(mapped + "\n".join(extras) + "\n")
    result = analysis.inspect_structure(
        str(path),
        protein_chains={"X": "MA", "Y": "MA"},
        ligand_chains={"L": "CO", "J": "CO"},
    )
    assert result["protein_residues"] == 4
    with pytest.raises(ValueError):
        analysis.inspect_structure(
            str(path),
            protein_chains={"X": "MA", "Y": "MA"},
            ligand_chains={"L": "CO", "J": "CC"},
        )


def test_boltz_request_uses_configured_chain_ids(tmp_path):
    import importlib.util

    module_path = next((LAB / "models/core/src").glob("*.py"))
    module_spec = importlib.util.spec_from_file_location("chain_core", module_path)
    core = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(core)
    predictor = core.Boltz2AffinityPredictor(protein_chain_id="X", ligand_chain_id="L")
    predictor._protein_sequence = "MA"
    predictor._ligand_smiles = "CO"
    request = yaml.safe_load(
        predictor._build_request_document({"use_msa_server": False})
    )
    assert request["sequences"][0]["protein"]["id"] == "X"
    assert request["sequences"][1]["ligand"]["id"] == "L"
    assert request["properties"][0]["affinity"]["binder"] == "L"
    with pytest.raises(ValueError):
        core.Boltz2AffinityPredictor(protein_chain_id="X", ligand_chain_id="X")
