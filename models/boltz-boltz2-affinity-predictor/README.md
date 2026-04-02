# Boltz: Boltz2AffinityPredictor

Native `biosim.BioModule` wrapper around `boltz predict` for a focused
Boltz-2 protein-ligand affinity workflow.

## What It Does

- accepts a single protein sequence and single ligand SMILES
- supports either:
  - an explicit `msa_path` input, or
  - `use_msa_server: true`
- writes a Boltz YAML request into a run directory
- bootstraps a local managed Boltz runtime on first use by default
- executes the managed `boltz` CLI
- exposes compact BioSignals for:
  - `affinity_summary`
  - `confidence_summary`
  - `structure_artifacts`
  - `run_metadata`

## Input Ports

| Input | Meaning |
|---|---|
| `protein_sequence` | Target protein sequence as a string or `{sequence: ...}` payload |
| `ligand_smiles` | Ligand SMILES string or `{smiles: ...}` payload |
| `msa_path` | Optional path to a precomputed MSA file |
| `run_options` | Optional per-run overrides such as `use_msa_server`, `output_format`, or `template_path` |

## Output Ports

| Output | Meaning |
|---|---|
| `affinity_summary` | Parsed affinity JSON payload from Boltz |
| `confidence_summary` | Parsed top-model confidence JSON payload from Boltz |
| `structure_artifacts` | Absolute paths to generated structure and summary artifacts |
| `run_metadata` | Command, output directory, status, and captured logs |

## Runtime Behavior

By default this module is self-contained at the wrapper level:

- it creates a local virtual environment under `.runtime/boltz2`
- it installs Boltz there on first run
- it then reuses that managed runtime on later runs

Important constraints:
- internet access is still required on the first run so the managed runtime can
  install Boltz and its dependencies
- GPU-oriented execution is still the intended default
- large model weights, MSAs, and generated structures are not checked into this
  repository

You can still opt into external-runtime behavior with `runtime_mode: external`,
but the default is `managed`.

## Examples

See:
- [`examples/boltz2-minimal`](/Volumes/dem-ssd/imp/projects/Nitoons/Biosimulant/models/models-boltz/examples/boltz2-minimal)
- [`examples/boltz2-explicit-msa`](/Volumes/dem-ssd/imp/projects/Nitoons/Biosimulant/models/models-boltz/examples/boltz2-explicit-msa)
- [`examples/boltz2-wiring`](/Volumes/dem-ssd/imp/projects/Nitoons/Biosimulant/models/models-boltz/examples/boltz2-wiring)
