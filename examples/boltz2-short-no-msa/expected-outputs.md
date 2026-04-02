# Expected Outputs

This example combines the upstream Boltz `examples/prot_no_msa.yaml` protein
with the ligand from `examples/affinity.yaml`.

Run it with:

```bash
python3 examples/run_example.py boltz2-short-no-msa
```

It emits the standard four BioSignals:

- `affinity_summary`
- `confidence_summary`
- `structure_artifacts`
- `run_metadata`

Expected real Boltz outputs include:
- `predictions/request/request_model_0.cif`
- `predictions/request/confidence_request_model_0.json`
- `predictions/request/affinity_request.json`

This is the best local smoke-test example in the repo because:
- it uses real upstream inputs rather than placeholders
- it avoids MSA generation with `msa_path: empty`
- it uses reduced sampling and recycling settings so CPU validation is more practical
