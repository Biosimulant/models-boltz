# models-boltz

> Publish-ready Boltz lab repo: each lab is self-contained under `labs/<slug>/`
> with a core Boltz runner, visualisation module, README, screenshots, and
> source-faithful caveats for generated structure and affinity-style outputs.

Curated collection of **Boltz-family biomolecular interaction models** for the
**biosim** platform.

This repository hosts native Python `biosim.BioModule` wrappers around
upstream Boltz runtimes. It is intended for multiple Boltz-family models over
time; the first implemented model is a Boltz-2 affinity-focused predictor for
protein-ligand runs.

## Runtime compatibility

The components shipped by the seven publish-ready Labs use
`BioModule.execute()` with `ExecutionPolicy.ONCE_BEFORE_RUN`. BioWorld invokes
each component once per run and drains the complete dependency chain in stable
order before temporal windows begin. Lab manifests remain unchanged for current
product compatibility; their short duration and settle fields no longer control
component invocation count. Biosimulant runtimes predating invocation-policy
support require the prior wrapper release.

This repo is distinct from
[`models-onnx`](https://github.com/Biosimulant/models-onnx):
- `models-onnx` is for checked-in ONNX artifacts wrapped behind BioModules
- `models-boltz` is for native Boltz-family runtime wrappers and file-backed
  structural outputs

## What's Inside

### Publish-Ready Labs

| Lab | Description |
|---|---|
| `boltz-boltz2-affinity-predictor` | Native Boltz-2 subprocess wrapper that predicts structure, confidence, and affinity for a single protein plus single ligand workflow. |
| `boltz-protein-ligand-101-workflow` | Guided protein-ligand walkthrough using a curated sequence and ligand example. |
| `boltz-cancer-kinase-workflow` | ABL1 kinase-domain workflow with an imatinib reference ligand example. |
| `boltz-antiviral-protease-workflow` | SARS-CoV-2 3C-like protease workflow with a nirmatrelvir reference ligand example. |
| `boltz-tuberculosis-target-workflow` | Mycobacterium tuberculosis InhA workflow with a triclosan reference ligand example. |
| `boltz-malaria-binding-workflow` | Plasmodium falciparum falcipain-2 workflow with an E-64 reference ligand example. |
| `boltz-batch-ligand-ranking-workflow` | Small-set ligand ranking workflow that compares three curated ABL1 ligands. |

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
- keeps a repo-local Boltz cache under `.runtime/boltz-cache` by default
- builds a Boltz YAML request inside a run directory
- invokes `boltz predict`
- retries once after purging known-corrupted cached Boltz assets
- emits compact BioSignals for affinity, confidence, structure artifacts, and
  run metadata
- emits a `structure3d` visualization payload for compatible BioSim web and
  desktop clients so the top-ranked complex can be inspected directly

## Prerequisites

The wrapper module itself uses the Python standard library plus `biosim`.
Actual Boltz execution no longer depends on a system-wide Boltz install by
default: the model bootstraps a local managed runtime under this repository on
first use and reuses a local cache directory on later runs.

```bash
pip install "biosim==0.0.7"
```

For real runs, the first execution still needs:
- internet access to install Boltz and its Python dependencies
- a working Python environment with `venv`
- suitable hardware for the requested accelerator mode

The default local paths are:
- managed runtime: `.runtime/boltz2`
- managed cache: `.runtime/boltz-cache`

If Boltz fails with a known corrupted-cache extraction error, the wrapper
purges the affected cached assets and retries once automatically.

## Remote Execution

For BioSim remote execution on Modal, the same package is used without a
Boltz-specific provider branch:

- the package manifest pins `boltz[cuda]==2.0.2` for GPU-backed sandbox
  dependency installs on Modal
- remote runs force `runtime_mode: external` through manifest-declared
  remote init overrides
- the Boltz cache is redirected to
  `${REMOTE_EXECUTION_MOUNT_ROOT}/runtime-cache/boltz`
- generated structure and JSON artifacts are staged back through the generic
  remote executor artifact flow so completed runs can still serve `structure3d`
  downloads after the sandbox exits

The supported validation target for real Boltz-2 runs is Linux + NVIDIA GPU on
Modal. Local macOS runs remain useful for smoke-testing the wrapper and example
inputs, but they are not the release-grade validation target.

## Validation

```bash
./.venv-check/bin/python scripts/validate_manifests.py
./.venv-check/bin/python scripts/check_entrypoints.py
bash scripts/check_public_boundary.sh
./.venv-check/bin/pytest -q models/boltz-boltz2-affinity-predictor/tests
```

Remote Modal smoke:

```bash
cd ../bsim-platform/backend
BIOSIM_MODAL_RUN_REAL_SMOKE=1 ./.venv/bin/pytest -q tests/test_executor_remote_overrides.py -k modal_real_boltz_smoke
```

## License

Dual-licensed: Apache-2.0 (code), CC BY 4.0 (content).
