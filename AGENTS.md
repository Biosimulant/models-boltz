# AGENTS.md

Instructions for AI agents working in the models repository.

## Repository Purpose

This is a **public** curated repository of Boltz-family biological model packs
for the [biosim](https://github.com/BioSimulant/biosim) platform. It provides
reusable, composable biomodules that wrap native Boltz runtimes for structural
biology and affinity workflows.

## Repository Structure

```
models/
├── models/          # Runnable model packages, each with a model.yaml manifest
├── examples/        # Example configs and wiring
├── templates/       # Starter template for new model packs
├── scripts/         # Validation and CI scripts
├── docs/            # Governance and policy documentation
└── .github/         # CI/CD workflows
```

## Key Concepts

- **Model Pack**: A self-contained simulation component with a `model.yaml` manifest, optional Python source, and dependencies. Each model defines inputs, outputs, and an `advance_to(t)` method.
- **Example Space**: A composed simulation defined by `space.yaml` that wires multiple models together with signal routing, runtime parameters, and initial conditions.
- **BioModule**: The base class (`biosim.BioModule`) that all custom models inherit from. It defines the `inputs()`, `outputs()`, `set_inputs()`, `advance_to()`, and `get_outputs()` interface.
- **Wiring**: YAML declarations in `space.yaml` that route signals from one model's output to another model's input (`from: module.output → to: module.input`).

## Manifest Schemas

### model.yaml (Version 2.0)

Required fields:
- `schema_version`: Must be `"2.0"`
- `title`, `description`: Human-readable metadata
- `standard`: Model standard (e.g., `"other"`)
- `biosim.entrypoint`: Python import path in `module.path:ClassName` or `module.path:function_name` format
- `biosim.init_kwargs`: Parameters passed to the model constructor
- `runtime.dependencies.packages`: All dependencies must be pinned with `==`

## Working with Models

### Adding a New Model

All new models and spaces must satisfy the acceptance criteria in [STANDARDS.md](STANDARDS.md). Use the checklist at the bottom of that file before submitting.

1. Copy `templates/model-pack/` to `models/<your-model-slug>/`
2. Edit `model.yaml` with your model's metadata, entrypoint, and dependencies
3. Implement your module by subclassing `biosim.BioModule`
4. Add tests in `tests/test_<module>.py` (see [STANDARDS.md § Tests](STANDARDS.md#tests))
5. Ensure all dependencies use exact version pinning (`==`)
6. Run validation: `python scripts/validate_manifests.py && python scripts/check_entrypoints.py`

### Model Implementation Notes

- Most curated models include Python source in `src/` (entrypoint `src.<module>:<ClassName>`).
- When making changes, keep `outputs()` accurate and keep example configs aligned with the real public interface.
- Boltz-family models should prefer file-backed structural outputs and compact summary signals.

## Validation and CI

Three validation scripts must pass before merging:

1. **`scripts/validate_manifests.py`** — Schema validation for all YAML manifests (required fields, dependency pinning, structure)
2. **`scripts/check_entrypoints.py`** — Verifies all Python entrypoints are importable and callable
3. **`scripts/check_public_boundary.sh`** — Scans markdown files for business-sensitive keywords that must not appear in this public repo

The CI pipeline (`.github/workflows/ci.yml`) runs these in order: secret scan → manifest validation → smoke sandbox (Docker).

## Public Boundary Rules

This repository is **public**. Do not add:

- Business strategy, pricing, GTM, or fundraising content
- Investor-facing materials or revenue forecasts
- Private operational runbooks or internal credentials
- Secrets or API keys (enforced by gitleaks and detect-secrets)

See [docs/PUBLIC_INTERNAL_BOUNDARY.md](docs/PUBLIC_INTERNAL_BOUNDARY.md) for the full policy.

## Code Style and Conventions

Full conventions, interface contracts, and acceptance checklists are in [STANDARDS.md](STANDARDS.md).

- All YAML manifests use `schema_version: "2.0"`
- Python dependencies must be pinned with exact versions (`==`)
- Custom Python modules follow the `biosim.BioModule` interface contract
- Model slugs use kebab-case with domain prefix (e.g., `neuro-`, `ecology-`, `virtualcell-`)
- Every model must include unit tests in `tests/test_<module>.py`
- Pre-commit hooks enforce trailing whitespace, EOF newlines, YAML syntax, and secret detection

## Domain Prefixes

| Prefix | Domain | Example |
|--------|--------|---------|
| `neuro-` | Neuroscience (spiking neurons, synapses, monitors) | `neuro-izhikevich-population` |
| `ecology-` | Ecosystem dynamics (populations, environment) | `ecology-predator-prey-interaction` |
| `virtualcell-` | Cellular and molecular biology | `virtualcell-grn-predictor`, `virtualcell-arc-state-predictor`, `virtualcell-perturbation-source` |

## Dependencies

- **Python 3.11+** required
- **biosim**: Core framework, installed from `git+https://github.com/BioSimulant/biosim.git@main`
- **pyyaml**: For manifest parsing in validation scripts
- An external **Boltz** installation is required for real structural runs
