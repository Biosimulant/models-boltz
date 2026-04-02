# Expected Outputs

This example uses the real upstream Boltz affinity sequence together with a
checked-in upstream-derived `seq1.a3m`.

Run it with:

```bash
python3 examples/run_example.py boltz2-explicit-msa
```

It emits the same four BioSignals as the minimal example:

- `affinity_summary`
- `confidence_summary`
- `structure_artifacts`
- `run_metadata`

Differences from the server-backed example:
- the checked-in `assets/seq1.a3m` is used directly
- `use_msa_server` stays disabled
- the generated Boltz request includes the resolved path to
  `examples/boltz2-explicit-msa/assets/seq1.a3m`
- the first managed-runtime bootstrap still needs internet access
