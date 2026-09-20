# Boltz with explicit scientific interfaces and output acceptance

A new, unpublished Lab based on the existing Boltz 2.0.2 adapter. The prediction
component sends probability, native affinity score and mmCIF structure through
five independently profiled connections (including the requested sequence and molecule) to an analysis component. The analysis
preserves each model-specific meaning; it does not convert scores to calibrated
Kd or experimental potency.

Requires the unreleased Biosimulant 0.0.36 implementation. Existing released Labs
and benchmark artifacts are unchanged. The analysis imports the acceptance API
at load time so an older runtime cannot silently ignore this contract.

The root manifest declares ten required profiled ports and seven output checks.
The seventeen resulting checks require exercised interfaces, finite/range-checked
scores, a basic mmCIF representation, predictor completion and downstream analysis
completion. Analysis parses atom-site records, requires finite coordinates and one
model with unique atom IDs, matches protein residue identity/coverage to the
requested sequence, and compares ligand heavy-element counts to RDKit-parsed
SMILES. These do not establish ligand connectivity, stereochemistry, binding or
structural accuracy. Same-formula isomers can pass the heavy-element check.

For scientific execution, install the released supporting runtime and the pinned
model dependencies, provide a matching protein/MSA and molecule, then run:

```bash
biosimulant labs validate ./labs/boltz-scientific-acceptance --json
biosimulant labs run ./labs/boltz-scientific-acceptance --require-acceptance --results-file result.json --json
```

The packaged default protein has a matching packaged A3M. Changing the protein
requires a matching MSA or explicitly selecting the MSA service. A real model run
requires an appropriate GPU. No new hardware compatibility qualification is
claimed by this Lab version.

`tests/test_workflow.py` executes actual package preparation, CLI dispatch, signal
routing and acceptance against an explicitly synthetic executable. It allocates
no GPU and must never be counted as a model prediction. Core adapter tests are
retained from the source adapter. No production publication has been performed.

For managed runs select `quality_profile_ref: research-ready@3`; deploy the
supporting runtime/backend before submitting. Review the acceptance result
separately from execution status.


## Explicit chain assignments

The predictor's `protein_chain_id` and `ligand_chain_id` parameters default to A
and B for existing inputs, but are configurable and must be distinct. The exact
assignments used in the Boltz request are emitted through `requested_chains` and
wired into analysis. The shared runtime has no chain-name rules.

`inspect_structure` also accepts explicit `protein_chains` and `ligand_chains`
mappings for multi-chain complexes, plus an explicit `residue_names` mapping for
nonstandard residues. Every declared chain must match; undeclared chains fail.
This analyzer still checks one coordinate model and uses mmCIF label sequence
numbering. The predictor's current single-protein/single-ligand affinity input
interface is preserved; a multi-chain analyzer is not a claim of GPU-qualified
multimer affinity inference. Heavy-element equality does not establish ligand
connectivity, stereochemistry, pose quality or experimental binding.
