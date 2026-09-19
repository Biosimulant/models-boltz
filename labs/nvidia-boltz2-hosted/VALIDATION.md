# BIO-001 implementation evidence

Local technical suite: 34 passed, 0 failed (Python 3.10.20,
biosimulant 0.0.34, httpx 0.28.1, biopython 1.84, PyYAML 6.0.2,
pytest 8.3.5). Tests cover request validation, native nullable embeddings,
manifest/runtime/Lab port agreement, finite BioWorld execution, transport errors,
single POST followed by polling, bound resume, credential handling, response
limits and corrupt outputs. Synthetic transport tests are explicitly labeled.

Command: `python -m pytest -q tests/test_adapter.py` from this Lab directory.
Local `validate_lab_source` also returned valid with no errors or warnings.

The genuine prior native fixture was parsed independently of transport mocks.
Its request and response digests match the earlier receipt, and all native
affinity fields remain unchanged. This exposed nullable embedding fields and
led to a regression test. Fixture parsing is not fresh inference.

## Fresh adapter execution

One new NVIDIA request was submitted by this hand-authored model inside local
BioWorld. The single model committed all five typed outputs. This was not a
Hub run and did not download weights or run a local GPU model.

- Provider request: `07d9b06e-9c11-4ed6-b547-3f7e36a92a9c`
- UTC: 2026-09-19T23:15:10.720277 to 2026-09-19T23:15:22.783662
- Measured client observation: 12.009977708 seconds
- Model.py SHA-256: `73b6eb771007a2caca2340c89d772e21a4b75ee9ec198c16ebbd9a36c501392d`
- Request SHA-256: `d6fc730297845b004a649ccedca8cfb94c72396ab6ab94fefee1fd848915408e`
- Raw response: 269530 bytes, SHA-256 `8426321d5618ffc78cf8dfcd37c42193ecff6151cf08665c475949e086ed652e`
- Captured typed results: 556605 bytes, SHA-256 `d06f2a7f3fb0d923175897f020b70f76d62d9c001d6686218ee45e358d01de2c`
- Structure SHA-256: `6be2bd5a8a1bc85271affb102927b96e4d1fcbaa86ac29eeb855ee7784ea71a8`

An independent inspection recomputed those digests and sizes, verified exact
raw-response preservation and all five terminal ports, parsed finite coordinates
and unique atom IDs, matched the entire 384-residue protein and checked ligand
heavy-element inventory C9/N1/O3 against RDKit parsing of the input SMILES.
All signal timestamps were the zero orchestration boundary, not biological time.

Evidence bundle: `boltz-adapter-live-001/{execution,request,response,results,inspection}.json`
and `sample-0.cif` in this task's output directory. The bundle is local evidence;
it is not asserted to be pinned in the owned Hub workspace.

## Outstanding evidence

Managed secret availability, exact managed environment, fresh Hub inference,
captured workspace-results and its Passport, private release, and public release
approval remain absent. Full ligand graph/stereochemistry, native chain-score
index mapping, pIC50 calibration, optional-feature domain coverage, scientific
benchmarks and acceptance thresholds remain unresolved. These gaps prevent
claiming completion or research qualification.
