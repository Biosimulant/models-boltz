# Boltz-2 Batch Ligand Ranking Workflow Runner

This model is the execution module for the `Boltz Workflow: Batch Ligand Ranking` lab. It reuses the native Boltz-2 CLI wrapper for a small, capped ligand set and exposes a narrow set of workflow inputs:

- `protein_sequence`
- `ligand_csv`
- `msa_path`
- `run_options`

The model emits a ranked batch summary plus the affinity, confidence, structure-artifact, and run-metadata records for the top-ranked completed ligand. Ranking is based on Boltz-2 binder probability and affinity-like value; it is an early comparison workflow, not experimental binding evidence.
