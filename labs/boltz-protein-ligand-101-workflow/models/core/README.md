# Boltz-2 Protein-Ligand Workflow Runner

This model is the execution module for the `Boltz Workflow: Protein-Ligand 101` lab. It reuses the native Boltz-2 CLI wrapper for a curated teaching protein-ligand example and exposes a narrow set of workflow inputs:

- `protein_sequence`
- `ligand_smiles`
- `msa_path`
- `run_options`

The model emits parsed affinity, confidence, structure-artifact, and run-metadata records for the downstream visualisation model. Outputs are Boltz-2 predictions for learning and early hypothesis generation; they are not experimental binding measurements.
