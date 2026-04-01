# Expected Outputs

This minimal example should emit four BioSignals:

- `affinity_summary`
- `confidence_summary`
- `structure_artifacts`
- `run_metadata`

`structure_artifacts` should contain absolute paths to the top structure file
plus the Boltz confidence and affinity JSON files.

This example assumes:
- the `boltz` CLI is installed
- GPU-oriented execution is available
- server-side MSA generation is acceptable for the run
