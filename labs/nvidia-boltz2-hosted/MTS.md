# BIO-001: hosted adapter technical design

Status: implementation design under review; no scientific approval asserted.
Implements the user-authorized Boltz-2 work. MRS scientific gates stay open.

One direct BioModule uses `ExecutionPolicy.ONCE_BEFORE_RUN`. Its only model
computation is a real NVIDIA HTTPS invocation; local code validates requests,
observes execution, parses and preserves results. The CPU runner does not load
weights. BioWorld timestamps are orchestration boundaries, not predicted time.

## Interface

`molecular_system`: required record with exact keys `polymers`, `ligands`,
`constraints`, each a JSON list following the provider schema. Entity IDs are
required by this adapter for stable mapping, including IDs on ligands. No
scientific default. A3M is inline text, never an arbitrary path or URL.

`run_options`: optional record with exact key `parameters` (JSON object).
Only the explicit allowlist in model.py is forwarded. Omitted values use provider
defaults and are recorded as omitted, not invented resolved backend settings.
There is no seed input because a hosted seed contract has not been verified.
Full PAE/PDE filesystem flags are not offered by this hosted interface: the
provider-side paths are not available as returned artifacts. This missing
capability stays recorded; no zero matrices substitute for them.

Outputs have exact outer record schemas; nested objects are validated in code:

| Port | Payload | Scientific units |
|---|---|---|
| predicted_complexes | items with inline mmCIF, SHA-256, bytes, sample_index, native name/source and chain_ids; request_sha256; coordinate_unit | coordinates in angstrom |
| affinity_predictions | by_ligand; field_definitions; request_sha256 | native log-micromolar-like score, probability and separate provider pIC50; no conversion |
| quality_metrics | values; field_units; request_sha256 | PDE/ipDE in angstrom; other native score scales preserved |
| run_provenance | request, digests, request ID, endpoint, UTC capture time, duration and backend identity gaps | wall_seconds in seconds; no biological time |
| native_response | exact decoded UTF-8 JSON response | heterogeneous archival data; no scientific compatibility profile |

Record envelope unit `1` denotes mixed terminal reports; per-field units stay
explicit. No draft compatibility profile is attached. Sequence, SMILES and A3M
primitive profiles need explicit entity-extraction adapters, and a multi-sample
mixed complex must not claim the narrower single protein-ligand file profile.

## Runtime and failure semantics

Pin `biosimulant==0.0.34`, `httpx==0.28.1`, `biopython==1.84`,
`pyyaml==6.0.2` and `pytest==8.3.5`; Python 3.10. The last two support the
collectable manifest/contract acceptance suite in the same declared environment.
This is the installed test environment, distinct from local source 0.0.35 and
the unverified deployed Gateway runtime. Managed preflight must accept it.

The credential is read only from NVIDIA_API_KEY inside execution. It is never
accepted as a model parameter. HTTPS redirects on authenticated calls are
disabled. Unexpected redirect/download responses fail without forwarding the
credential or logging URLs. Such results require an explicit artifact adapter.
No automatic retry of a POST, including transport failures. Explicit resume
requires both the original provider UUID and canonical request SHA-256 from the
failed observation. The supplied request must reproduce that digest before GET.

Default observation budget: 1200 seconds; per-request timeout: at most 60 seconds
and bounded by remaining observation budget. Response limit: 20 MB decoded
bytes. These are client bounds, not NVIDIA GPU time/cost guarantees. Provider
execution may continue after an observation failure. Exceptions carry only safe
status/request identifiers; clients may resume the same UUID explicitly.

Raw JSON and mmCIF are inline portable outputs. Successful managed execution must
capture them in workspace-results with independently verified outer artifact
checksum. Failure paths expose a safe provider ID in the error rather than emit
plausible scientific outputs. Large artifact redirects are an unresolved feature.

## Verification versus qualification

Tests separately label transport mocks, captured-response parser fixtures and
actual live inference. No fixture is bundled as a runtime fallback. Exact raw
byte preservation is tested by digest; no stochastic cross-run numerical
tolerance is invented. Full chemical identity, scientific benchmarks, exact
backend versioning, secure Hub secret provisioning and public rights review
remain blockers to completion even if local technical tests pass.
