# models-boltz

Curated collection of **Boltz-family biomolecular interaction models** for the
**biosim** platform.

This repository hosts native Python `biosim.BioModule` wrappers around
upstream Boltz runtimes. It is intended for multiple Boltz-family models over
time; the first implemented model is a Boltz-2 affinity-focused predictor for
protein-ligand runs.

This repo is distinct from
[`models-onnx`](https://github.com/Biosimulant/models-onnx):
- `models-onnx` is for checked-in ONNX artifacts wrapped behind BioModules
- `models-boltz` is for native Boltz-family runtime wrappers and file-backed
  structural outputs

## What's Inside

### Models

| Model | Description |
|---|---|
| `boltz-boltz2-affinity-predictor` | Native Boltz-2 subprocess wrapper that predicts structure, confidence, and affinity for a single protein plus single ligand workflow. |

### Repo Scope

This repository is for:
- native Boltz-family wrappers that implement the `biosim.BioModule` contract
- file-backed structural-biology runs that emit summaries plus artifact paths
- reusable examples showing how to wire Boltz-family modules into BioSim

This repository is not for:
- ONNX exports of Boltz models
- checked-in large biological inputs such as full MSAs or structure outputs
- unrelated structural-biology runtimes that do not belong to the Boltz family

## How It Works

The initial Boltz-2 model:
- accepts a protein sequence and ligand SMILES as BioSim inputs
- supports either a provided MSA path or Boltz server-side MSA generation
- bootstraps a local managed Boltz runtime on first use by default
- builds a Boltz YAML request inside a run directory
- invokes `boltz predict`
- emits compact BioSignals for affinity, confidence, structure artifacts, and
  run metadata

## Examples

The repo includes three examples under [`examples/`](examples):
- `boltz2-minimal` for a server-side MSA workflow
- `boltz2-explicit-msa` for a provided MSA path workflow
- `boltz2-wiring` for a BioSim `space.yaml` wiring example

See [`examples/README.md`](examples/README.md) for the example inventory.

## Prerequisites

The wrapper module itself uses the Python standard library plus `biosim`.
Actual Boltz execution no longer depends on a system-wide Boltz install by
default: the model bootstraps a local managed runtime under this repository on
first use.

```bash
pip install "biosim @ git+https://github.com/BioSimulant/biosim.git@main"
```

For real runs, the first execution still needs:
- internet access to install Boltz and its Python dependencies
- a working Python environment with `venv`
- suitable hardware for the requested accelerator mode

## Validation

```bash
python scripts/validate_manifests.py
python scripts/check_entrypoints.py
bash scripts/check_public_boundary.sh
pytest
```

## License

Dual-licensed: Apache-2.0 (code), CC BY 4.0 (content).
