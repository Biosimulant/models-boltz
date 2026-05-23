# Boltz Workflow: Malaria Binding

Malaria Binding Workflow is a curated BioSimulant Boltz workflow. It runs a single curated target sequence and a single reference ligand SMILES through Boltz-2, then presents the predicted complex, binder probability, affinity-like value, confidence metrics, run metadata, and report-ready caveats in the standard Biosimulant lab interface.

The workflow is designed for learning, comparison against known examples, and early biological hypothesis generation. It is not experimental validation, a clinical prediction, or a replacement for docking review, MD/FEP, assay design, or wet-lab confirmation.

## Workflow Status

This lab validates locally, exports as a portable `.bsilab` package, and has a successful private pre-publication GPU-backed run. It is not published yet. Publication should wait for explicit approval and then the Hub workflow card can be switched from `Coming soon` to the published lab id.

Publication checklist:

- manifest validation passes: complete
- strict package export passes: complete
- entrypoints import successfully: complete
- unit tests pass: complete
- at least one real GPU run completes: complete
- run results include structure, affinity, confidence, metadata, and visuals: complete
- screenshots/assets are captured from the real run: complete
- Hub workflow card is updated with the published lab id: pending


## Pre-Publication Run Evidence

The current assets and metrics are derived from this private staged remote run:

- run id: `bdda5ff7-e2ff-4621-8377-1e32ae3e7168`
- staged lab id: `ab1ef5e4-1940-457d-8bb3-50ec0a6d923a`
- run lab commit: `c3215e3019974f50b537b714de23b0c7afabb2a26d8882324a88e215ede00df8`
- remote size: GPU A10G
- status: completed
- duration: 407.1 seconds
- credits settled: 35

Key outputs from the run:

- `affinity_probability_binary`: `0.8718628883361816`
- `affinity_pred_value`: `-0.3072456419467926`
- `confidence_score`: `0.959920346736908`
- `complex_plddt`: `0.9600638747215272`
- `iptm`: `0.9593459367752076`
- `ligand_iptm`: `0.9593459367752076`

## Curated Known Example

The bundled known-example mode starts with:

- target: Plasmodium falciparum falcipain-2 cysteine protease
- target source: RCSB PDB `3BPF` sequence FASTA
- ligand: E-64 reference cysteine-protease inhibitor
- ligand source: PubChem CID `123985`
- protein sequence length: 241 amino acids
- MSA server usage enabled for the default example

The bundled E-64 example is a reference inhibitor-like teaching ligand, not a clinical antimalarial recommendation.

## What This Workflow Does

When run, the workflow:

1. Builds a Boltz-2 request from the curated or user-supplied protein and ligand inputs.
2. Runs the Boltz-2 CLI on a GPU-backed runtime.
3. Parses the top-ranked structure artifact.
4. Parses Boltz-2 affinity and confidence summaries.
5. Emits Biosimulant visuals and report-ready outputs.

## Inputs

- `protein_sequence`: amino-acid sequence for the target protein. If omitted, the bundled known example is used.
- `ligand_smiles`: SMILES string for the ligand. If omitted, the bundled known example ligand is used.
- `msa_path`: optional path to a precomputed `.a3m` MSA file.
- `run_options`: optional record for workflow/runtime options.

The known example mode works because `lab.yaml` defines defaults directly on the Boltz runner model. A new user can click Run without knowing YAML, SMILES formatting details, or Boltz CLI arguments.

## Outputs

- `structure_artifacts`: paths to the predicted complex structure files, usually mmCIF.
- `affinity_summary`: affinity-style outputs parsed from Boltz-2.
- `confidence_summary`: model confidence outputs for the top-ranked prediction.
- `run_metadata`: execution metadata, output paths, captured logs, and status.

The visualisation model turns these records into standard Biosimulant run visuals, including a structure viewer and summary metrics.

## Reading The Affinity Outputs

Boltz-2 distinguishes two affinity-oriented outputs that should not be collapsed into one meaning.

`affinity_probability_binary` is most useful as a binder-vs-decoy style signal. In product language, it is the binder probability. It is useful when comparing likely binders against unlikely binders under a hit-discovery framing.

`affinity_pred_value` is intended for ligand-optimization style use cases. In product language, it is an affinity-like value. It should be used cautiously and comparatively, not as a direct experimental measurement.

Both outputs are model predictions. They can help rank hypotheses for follow-up, but they do not prove binding, potency, selectivity, mechanism, or biological effect.

## Reading The Structure

The structure viewer shows the predicted protein-ligand complex from the top-ranked Boltz output. Use it as a sanity check:

- Is the ligand near a plausible pocket?
- Is the interface confidence reasonable?
- Does the pose look inconsistent with the score?
- Are there warnings in run metadata?

A high binder probability with an implausible pose should be treated as suspicious. A plausible pose with weak confidence should also be treated as uncertain.

## Safe Use Cases

This workflow is best for malaria-related teaching examples, parasite protease target context, and early hypothesis generation.

Safe uses include:

- Learn protein-ligand modeling workflows.
- Compare known ligand examples.
- Generate early hypotheses.
- Prepare candidate lists for deeper review.
- Create reproducible computational biology reports.
- Teach structure-based drug-discovery concepts.

## Do Not Claim

- Validated drug discovery.
- Clinical prediction.
- Medical diagnosis.
- Final compound selection.
- Wet-lab replacement.
- FEP replacement.
- Experimentally guaranteed binding-affinity certainty.

## Assets

The current screenshots were captured from the successful private pre-publication GPU run above, using its persisted mmCIF structure artifact and parsed run metrics.

![Boltz-2 predicted protein-ligand complex structure](assets/boltz2-affinity-structure-results.png)

![Boltz-2 affinity and confidence summary metrics](assets/boltz2-affinity-summary-metrics.png)

## Implementation Notes

This workflow intentionally reuses the existing Boltz-2 affinity predictor and visualisation modules. The product difference is in the lab packaging:

- guided disease or target context
- workflow tags
- known example defaults
- safe input names
- report-oriented caveats
- future Hub placement in the Boltz Workflows section
