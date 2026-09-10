from __future__ import annotations

from biosim.modules import ExecutionContext, ExecutionPolicy

from src.boltz2_batch_ligand_ranker import Boltz2BatchLigandRanker


def test_parse_ligand_csv_accepts_name_smiles_and_metadata(tmp_path):
    module = Boltz2BatchLigandRanker(work_dir=str(tmp_path))

    ligands = module._parse_ligand_csv(
        "name,smiles,source\n"
        "Imatinib,CCN,PubChem\n"
        "Dasatinib,CCC,PubChem\n"
    )

    assert [ligand["name"] for ligand in ligands] == ["Imatinib", "Dasatinib"]
    assert [ligand["smiles"] for ligand in ligands] == ["CCN", "CCC"]
    assert "PubChem" in ligands[0]["metadata"]


def test_rank_rows_uses_binder_probability_then_affinity_value(tmp_path):
    module = Boltz2BatchLigandRanker(work_dir=str(tmp_path))
    rows = [
        {"ligand": "low", "binder_probability": 0.2, "affinity_like_value": 10.0},
        {"ligand": "top", "binder_probability": 0.8, "affinity_like_value": 1.0},
        {"ligand": "tie-break", "binder_probability": 0.8, "affinity_like_value": 2.0},
    ]

    ranked = module._rank_rows(rows)

    assert [row["ligand"] for row in ranked] == ["tie-break", "top", "low"]


def test_missing_csv_surfaces_error_payload(tmp_path):
    module = Boltz2BatchLigandRanker(default_protein_sequence="MKT", work_dir=str(tmp_path))

    outputs = module.execute({}, context=ExecutionContext(policy=ExecutionPolicy.ONCE_BEFORE_RUN, run_start=0.0, run_end=0.1))

    metadata = module._output_payloads["run_metadata"]
    assert metadata["status"] == "error"
    assert "ligand_csv" in metadata["error"]
