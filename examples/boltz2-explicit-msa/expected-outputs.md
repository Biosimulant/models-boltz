# Expected Outputs

This example emits the same four BioSignals as the minimal example:

- `affinity_summary`
- `confidence_summary`
- `structure_artifacts`
- `run_metadata`

Differences from the server-backed example:
- `msa_path` must exist locally
- `use_msa_server` stays disabled
- the generated Boltz request includes `msa: /absolute/path/to/query.a3m`
- the first managed-runtime bootstrap still needs internet access
