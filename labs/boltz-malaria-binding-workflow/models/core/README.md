# Boltz-2 Malaria Binding Workflow Runner

This model is the execution module for the `Boltz Workflow: Malaria Binding` lab. It reuses the native Boltz-2 CLI wrapper for a curated Plasmodium falciparum falcipain-2 sequence and E-64 reference inhibitor-like ligand example, with these workflow inputs:

- `protein_sequence`
- `ligand_smiles`
- `msa_path`
- `run_options`

The model emits parsed affinity, confidence, structure-artifact, and run-metadata records for the downstream visualisation model. Outputs are Boltz-2 predictions for learning and early hypothesis generation; they are not antimalarial efficacy or experimental affinity measurements.
