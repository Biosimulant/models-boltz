# Examples

This directory contains reference examples for the first Boltz-family BioSim
module in this repository.

## Included Examples

### `boltz2-minimal`
- real upstream Boltz affinity example
- uses `use_msa_server: true`
- best starting point when you do not have a precomputed MSA

### `boltz2-explicit-msa`
- same real upstream affinity example with a checked-in upstream-derived A3M
- useful for reproducible precomputed-MSA workflows

### `boltz2-wiring`
- a `space.yaml` with real constructor defaults so the example is not just a
  placeholder manifest reference

## Exact Commands

Run the server-side MSA example:

```bash
python3 examples/run_example.py boltz2-minimal
```

Run the explicit-MSA example:

```bash
python3 examples/run_example.py boltz2-explicit-msa
```

Write results to a JSON file:

```bash
python3 examples/run_example.py boltz2-minimal --output-json examples/boltz2-minimal/latest-output.json
```

Use isolated work and runtime directories:

```bash
python3 examples/run_example.py boltz2-minimal \
  --work-dir /tmp/models-boltz-runs \
  --runtime-dir /tmp/models-boltz-runtime
```

## Notes

- These examples use real upstream Boltz example inputs rather than placeholders.
- They stay aligned with the current `Boltz2AffinityPredictor` constructor,
  input ports, and output ports.
- The first real run bootstraps a managed local Boltz runtime unless you
  explicitly switch to `runtime_mode: external`.
- CPU is used in the checked-in example configs for portability; real runs may
  still take a long time.
