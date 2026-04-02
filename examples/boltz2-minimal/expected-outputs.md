# Expected Outputs

This example is based on the upstream Boltz `examples/affinity.yaml` input.

Run it with:

```bash
python3 examples/run_example.py boltz2-minimal
```

It should emit four BioSignals:

- `affinity_summary`
- `confidence_summary`
- `structure_artifacts`
- `run_metadata`

`structure_artifacts` should contain absolute paths to the top structure file
plus the Boltz confidence and affinity JSON files.

Expected real Boltz outputs include:
- `predictions/request/request_model_0.cif`
- `predictions/request/confidence_request_model_0.json`
- `predictions/request/affinity_request.json`

This example assumes:
- internet access is available for the first managed-runtime bootstrap
- server-side MSA generation is acceptable for the run
- CPU execution is allowed, though it may be slow
