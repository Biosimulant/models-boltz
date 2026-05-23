# Boltz-2 Antiviral Protease Workflow Runner

This model is the execution module for the `Boltz Workflow: Antiviral Protease` lab. It reuses the native Boltz-2 CLI wrapper for a curated SARS-CoV-2 3C-like protease sequence and nirmatrelvir reference ligand example, with these workflow inputs:

- `protein_sequence`
- `ligand_smiles`
- `msa_path`
- `run_options`

The model emits parsed affinity, confidence, structure-artifact, and run-metadata records for the downstream visualisation model. Outputs are Boltz-2 predictions for learning and early hypothesis generation; they are not clinical efficacy or experimental affinity measurements.
