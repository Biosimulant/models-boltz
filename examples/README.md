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

### `boltz2-short-no-msa`
- short real upstream protein plus real upstream affinity ligand
- uses `msa_path: empty` from upstream Boltz input conventions
- best local smoke-test path for CPU-only validation

### `boltz2-wiring`
- a fully owned checked-in `lab.yaml` source example for the portable
  Boltz-2 remote package, not just a placeholder manifest reference

## Notes

- These examples use real upstream Boltz example inputs rather than placeholders.
- They stay aligned with the current `Boltz2AffinityPredictor` constructor,
  input ports, and output ports.
- `boltz2-wiring` is the checked-in source-of-truth for the owned remote
  example space that the export helper packages into a `.bsilab`.
- The first real run bootstraps a managed local Boltz runtime unless you
  explicitly switch to `runtime_mode: external`.
- The standalone run-example configs default to CPU for portability, but the
  owned remote space example is configured for GPU-backed remote execution.
- `boltz2-short-no-msa` is the most practical end-to-end smoke example on a local machine.
- The two affinity-length examples may still take a long time on CPU.
