# Boltz-2 Batch Ligand Ranking Workflow Runner

This model is the execution module for the `Boltz Workflow: Batch Ligand Ranking` lab. It reuses the native Boltz-2 CLI wrapper for a small, capped ligand set and exposes a narrow set of workflow inputs:

- `protein_sequence`
- `ligand_csv`
- `msa_path`
- `run_options`

The model emits a ranked batch summary plus the affinity, confidence, structure-artifact, and run-metadata records for the top-ranked completed ligand. Active-affinity mode ranks lower affinity-like values first, using binder probability to break ties. Binder-screen mode ranks higher predicted binder probability first. Both are exploratory comparisons, not experimental binding evidence.

Every evaluated ligand records observed checkpoint hashes and the environment of
the executed managed Python runtime. Raw alignment bytes and processed alignment
arrays are retained in an inspectable ZIP with their hashes and the exact request.
The batch exposes alignment archives for all candidates, although predicted
structure delivery remains limited to the top candidate. Missing provenance is
reported explicitly as incomplete; it does not invalidate or qualify the numerical
prediction by itself. External executable mode cannot reliably identify its Python
interpreter and therefore reports that environment gap. Captured hashes identify
observed bytes, not upstream authenticity, deterministic reruns or experimental
accuracy. Archives exceeding the default 64 MiB artifact limit remain an explicit
delivery gap.
