"""Declared Boltz output checks; no potency, stereochemistry or accuracy claim."""

import math
from collections import Counter

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from biosimulant import BioModule, ExecutionPolicy, SignalSpec, unwrap_payload
from biosimulant.acceptance import SCHEMA_VERSION
from rdkit import Chem


def inspect_structure(
    path,
    sequence=None,
    smiles=None,
    *,
    protein_chains=None,
    ligand_chains=None,
    residue_names=None,
):
    """Validate declared chain contents, including multiple proteins and ligands.

    residue_names explicitly maps nonstandard residue names to one-letter codes.
    This checks sequence and heavy-element inventory, never ligand connectivity.
    """
    proteins = dict(protein_chains) if protein_chains is not None else {"A": sequence}
    ligands = dict(ligand_chains) if ligand_chains is not None else {"B": smiles}
    if not proteins or not ligands or set(proteins) & set(ligands):
        raise ValueError("Declare disjoint protein and ligand chain assignments")
    if any(
        not isinstance(k, str) or not k or not isinstance(v, str) or not v
        for k, v in {**proteins, **ligands}.items()
    ):
        raise ValueError(
            "Chain assignments and requested identities must be nonempty strings"
        )
    residue_names = residue_names or {}
    if any(
        not isinstance(k, str) or not isinstance(v, str) or len(v) != 1
        for k, v in residue_names.items()
    ):
        raise ValueError("Residue mappings require explicit single-letter identities")
    cif = MMCIF2Dict(path)
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
    values = [cif.get("_atom_site." + column) for column in columns]
    if any(not isinstance(column, list) or not column for column in values):
        raise ValueError("Structure is missing required atom-site columns")
    if len({len(column) for column in values}) != 1:
        raise ValueError("Atom-site columns have different lengths")
    if len(set(values[-1])) != 1 or len(set(values[-2])) != len(values[-2]):
        raise ValueError("Expected one structure model with unique atom IDs")
    residues = {chain: {} for chain in proteins}
    ligand_elements = {chain: Counter() for chain in ligands}
    for chain, residue_id, residue_name, element, x, y, z, _atom, _model in zip(
        *values
    ):
        if not all(math.isfinite(float(v)) for v in (x, y, z)):
            raise ValueError("Nonfinite structure coordinate")
        if chain in proteins:
            index = int(residue_id)
            if index in residues[chain] and residues[chain][index] != residue_name:
                raise ValueError("Conflicting residue identities")
            residues[chain][index] = residue_name
        elif chain in ligands:
            symbol = element.capitalize()
            if symbol not in {"H", "D"}:
                ligand_elements[chain][symbol] += 1
        else:
            raise ValueError("Unexpected chain outside declared chain assignments")
    for chain, expected_sequence in proteins.items():
        if sorted(residues[chain]) != list(range(1, len(expected_sequence) + 1)):
            raise ValueError(
                "Protein residue coverage does not match the requested sequence"
            )
        predicted = "".join(
            residue_names.get(name, seq1(name, undef_code="?"))
            for _, name in sorted(residues[chain].items())
        )
        if "?" in predicted or predicted != expected_sequence:
            raise ValueError(
                "Predicted protein sequence differs from the requested sequence"
            )
    for chain, expected_smiles in ligands.items():
        molecule = Chem.MolFromSmiles(expected_smiles)
        if molecule is None:
            raise ValueError("Requested ligand SMILES cannot be parsed")
        expected_elements = Counter(
            atom.GetSymbol() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1
        )
        if not expected_elements or ligand_elements[chain] != expected_elements:
            raise ValueError(
                "Predicted ligand heavy-element composition differs from the requested molecule"
            )
    return {
        "protein_residues": sum(len(v) for v in residues.values()),
        "ligand_heavy_atoms": sum(sum(v.values()) for v in ligand_elements.values()),
        "protein_chains": sorted(proteins),
        "ligand_chains": sorted(ligands),
        "atom_count": len(values[0]),
        "finite_coordinates": True,
        "protein_sequence_matches": True,
        "ligand_heavy_elements_match": True,
    }


class BoltzAnalysis(BioModule):
    execution_policy = ExecutionPolicy.ONCE_BEFORE_RUN

    def inputs(self):
        return {
            "requested_chains": SignalSpec.record(
                schema={"protein": "str", "ligand": "str"}
            ),
            "binding_probability": SignalSpec.scalar(
                dtype="float64",
                accepted_units=["1"],
                contract={"profile": "boltz.binding-probability/v1"},
            ),
            "affinity_log10_ic50_micromolar": SignalSpec.scalar(
                dtype="float64",
                accepted_units=["1"],
                contract={"profile": "boltz.log10-ic50-micromolar/v1"},
            ),
            "predicted_structure": SignalSpec.scalar(
                dtype="str",
                value_type="file",
                format="mmcif",
                contract={"profile": "protein-ligand.complex-structure-mmcif/v1"},
            ),
            "requested_protein_sequence": SignalSpec.scalar(
                dtype="str",
                format="sequence",
                contract={"profile": "protein.sequence/v1", "species": "any"},
            ),
            "requested_ligand_smiles": SignalSpec.scalar(
                dtype="str", format="smiles", contract={"profile": "chemical.smiles/v1"}
            ),
        }

    def outputs(self):
        return {
            "completed": SignalSpec.scalar(dtype="bool"),
            "summary": SignalSpec.record(schema={"payload": "json"}),
        }

    def execute(self, inputs, *, context):
        values = {name: unwrap_payload(inputs[name]) for name in self.inputs()}
        summary = {
            **values,
            "acceptance_schema": SCHEMA_VERSION,
            "interpretation": "Model predictions only. Ligand connectivity, stereochemistry, binding and calibrated Kd are not established.",
        }
        try:
            summary["structure_checks"] = inspect_structure(
                values["predicted_structure"],
                values["requested_protein_sequence"],
                values["requested_ligand_smiles"],
                protein_chains={
                    values["requested_chains"]["protein"]: values[
                        "requested_protein_sequence"
                    ]
                },
                ligand_chains={
                    values["requested_chains"]["ligand"]: values[
                        "requested_ligand_smiles"
                    ]
                },
            )
            completed = True
        except (ValueError, KeyError, TypeError, OSError) as exc:
            summary["validation_error"] = str(exc)
            completed = False
        return {"completed": completed, "summary": {"payload": summary}}
