# Examples

This directory contains reference examples for the first Boltz-family BioSim
module in this repository.

## Included Examples

### `boltz2-minimal`
- smallest practical configuration
- uses `use_msa_server: true`
- good starting point when you do not have a precomputed MSA

### `boltz2-explicit-msa`
- same module with an explicit `msa_path`
- useful for reproducible offline or precomputed-MSA workflows

### `boltz2-wiring`
- a `space.yaml` showing how the model is referenced and wired in a BioSim
  composition

## Notes

- These examples intentionally use placeholder paths and short sequences.
- They are designed to stay aligned with the actual constructor parameters,
  input ports, and output ports of `Boltz2AffinityPredictor`.
- They do not ship large MSAs, structures, or Boltz outputs.
- The first real run will bootstrap a managed local Boltz runtime unless you
  explicitly switch to `runtime_mode: external`.
