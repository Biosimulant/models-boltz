# Boltz Workflow: Batch Ligand Ranking

Batch Ligand Ranking is a guided Biosimulant Boltz workflow for comparing a small ligand CSV against one protein target. It runs Boltz-2 once per ligand, extracts binder probability, affinity-like value, confidence metrics, and top-structure artifacts, then ranks the completed candidates into a report-ready table.

This is the first workflow here that is more than a repackaged single Boltz run. It adds CSV intake, repeated execution, result aggregation, ranking, flags, and a batch-specific visualisation table. It is still designed for learning, small-set comparison, and early biological hypothesis generation, not validated drug discovery or final compound selection.

## Execution behavior

Every graph component implements `BioModule.execute()` with
`ExecutionPolicy.ONCE_BEFORE_RUN`. BioWorld invokes each component once per
run and drains the dependency graph in stable layers; no artificial settle turn
is required to move data between these components. Existing temporal manifest
fields remain unchanged for compatibility with current Biosimulant products.

## Workflow Status

The Hub listing is private pending repair and validation. This working revision
fixes task-specific ranking, invalid-score handling, partial batches and downstream
reporting. Local tests exercise the five-stage graph using recorded GPU scores
and simulated failures. They do not validate a new Boltz prediction or a biological
rank order. GPU execution and publication of this repaired revision remain pending.

Replay provenance: run `df60d769-d0bb-492f-8219-915e81a2a1da`, artifact
`workspace-results`, independently verified SHA-256
`7ffb1ed1a5646e0ec4c9b9335e844e7335f9a747ce588b2a345a5c469ff6d8d9`.
Using its recorded affinity values gives Dasatinib, Imatinib, Nilotinib. The former
probability-based order was Dasatinib, Nilotinib, Imatinib. Both score types remain
available, with an explicit scientific task selector.

Remaining release work includes an end-to-end run in the declared GPU environment,
a run-wide execution budget, and verification that user overrides replace curated
provenance consistently. This Lab stays private until those checks are resolved.

## Historical Run Evidence (earlier revision)

The current assets and metrics are derived from this private staged remote run:

- run id: `c180d5ab-7d17-4037-8f4a-7b28aaf35a22`
- staged lab id: `3c6200af-bc4a-4041-9714-121954af1520`
- run lab commit: `f851d8490e3b725fdb685b1b998d1e96c0f6311191bef82067d60cfdacc55c89`
- remote size: GPU A10G
- status: completed
- duration: 873.2 seconds
- credits settled: 70

Key outputs from the run:

- completed ligands: `3`
- failed ligands: `0`
- top ligand: `Dasatinib`
- top `affinity_probability_binary`: `0.954538881778717`
- top `affinity_pred_value`: `-2.822425127029419`
- top `confidence_score`: `0.945051610469818`
- top `complex_plddt`: `0.9359697103500366`
- top `ligand_iptm`: `0.981379210948944`

Ranked table evidence:

- rank 1: Dasatinib - binder probability `0.954538881778717`, affinity-like value `-2.822425127029419`, confidence `0.945051610469818`
- rank 2: Nilotinib - binder probability `0.7802920341491699`, affinity-like value `-0.9997102618217468`, confidence `0.9546635746955872`
- rank 3: Imatinib - binder probability `0.6947591304779053`, affinity-like value `-0.9419443607330322`, confidence `0.9453792572021484`

## Curated Known Example

The bundled known-example mode starts with:

- target: Human ABL1 kinase domain
- target source: RCSB PDB `2HYY` sequence FASTA
- ligand CSV examples: Imatinib, Dasatinib, and Nilotinib
- ligand sources: PubChem CID `5291`, `3062316`, and `644241`
- protein sequence length: 273 amino acids
- maximum default ligands per run: 3
- MSA server usage enabled for the default example

The example is a kinase-inhibitor teaching set. It does not predict patient response, clinical efficacy, resistance, dosing, safety, or therapeutic suitability.

<!-- BIOSIMULANT_WORKFLOW_GRAPH_START -->
## Compose Workflow Graph

The published workflow is intentionally split into real Biosimulant modules:

1. `batch_target_context` emits source-backed target, ligand, disease/use-case, provenance, and caveat context.
2. `ligand_library_loader` resolves the public protein, ligand, MSA, and run-option inputs into the exact Boltz request.
3. `boltz_batch_ligand_ranker` runs the unchanged Boltz-2 scientific wrapper.
4. `ranking_interpreter` converts raw Boltz outputs into conservative evidence fields without adding new biological claims.
5. `visualisation` renders the 3D structure, confidence/affinity summaries, source context, request traceability, and Q/A caveat cards.

This makes the Compose view match the workflow promise while keeping Boltz-2 as the only predictive scientific model. The surrounding modules are provenance, request assembly, interpretation, and presentation stages.

<!-- BIOSIMULANT_WORKFLOW_GRAPH_END -->

## Inputs

- `protein_sequence`: amino-acid sequence for the shared target protein. If omitted, the bundled ABL1 kinase-domain example is used.
- `ligand_csv`: CSV text with `name` and `smiles` columns. Optional metadata columns are retained in the run context.
- `msa_path`: optional path to a precomputed `.a3m` MSA file.
- `run_options`: optional record for workflow/runtime options.

The known example mode works because `lab.yaml` defines the target and ligand CSV defaults on the input assembler. A new user can click Run without knowing YAML, SMILES formatting details, or Boltz CLI arguments.

## Outputs

- `batch_summary`: ranked ligand rows with binder probability, affinity-like value, confidence, status, and flags.
- `structure_artifacts`: paths to the top-ranked ligand predicted complex structure files, usually mmCIF.
- `affinity_summary`: affinity-style outputs for the top-ranked completed ligand.
- `confidence_summary`: model confidence outputs for the top-ranked completed ligand.
- `run_metadata`: batch execution metadata, per-ligand statuses, output paths, captured logs, and status.

The visualisation model turns these records into standard Biosimulant run visuals, including a top-ranked structure viewer and a batch ranking table.

## Ranking Semantics

The default `active_affinity` mode compares the bundled known ABL1 inhibitors by
`affinity_pred_value` ascending (lower predicts stronger binding), then binder
probability descending to break ties. It assumes the submitted compounds are
known active against the submitted target. For a mixed binder/decoy library, set
`run_options.ranking_mode` to `binder_probability`; this sorts by probability
descending, then affinity value ascending. These modes answer different questions.

Boltz reports affinity as `log10(IC50 / micromolar)`, not kcal/mol. Binder
probability is a separate classification score and is not a potency measurement.
See the [upstream score definitions](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).
Scores are model predictions, not experimental measurements or uncertainty bounds.
The default uses one diffusion sample; the displayed precision does not establish
ranking reproducibility.

Only completed runs with finite affinity values and probabilities in [0, 1] receive
ranks. Failed rows remain visible with an error and no rank. A batch with failures
is marked `partial`; an entirely failed batch has no top ligand. Submitted,
evaluated, completed and failed counts are reported. Blank SMILES or more than
`max_ligands` (default 3) reject the request before computation; no rows are silently
omitted. Pose review reminders use no uncalibrated binding/confidence thresholds.

Runtime setup, prediction and retry subprocesses share a 1,500-second batch
budget, leaving time within the managed 1,800-second limit to report results.
The budget starts at the batch stage; environment startup and platform scheduling
are outside this model's control. Remaining ligands are marked `not_started` when
the budget expires and receive no score or rank. `evaluated_count` counts started
ligands; `failed_count` includes every noncompleted row and `not_started_count`
identifies those never attempted. These are partial comparisons, not full-library
rankings. A cold runtime or slow external MSA service may consume the budget.

The assembled request records hashes of the actual protein and ligand library.
Changed inputs lose inherited example names and source claims. Reports use this
resolved context; the separate context-stage output is explicitly the packaged
example. User-supplied names and source metadata remain unverified. Explicitly
blank inputs stay blank and fail validation rather than selecting the example.

## Safe Use Cases

- Compare a small known ligand set against one target.
- Learn how Boltz-2 affinity outputs behave across ligands.
- Generate early candidate lists for deeper review.
- Prepare reproducible computational biology reports.
- Teach structure-based screening concepts.

## Do Not Claim

- Validated drug discovery.
- Clinical prediction.
- Medical diagnosis.
- Final compound selection.
- Wet-lab replacement.
- FEP replacement.
- Experimentally guaranteed binding-affinity certainty.

## Assets

The historical screenshots below were captured from the successful private pre-publication GPU run above, using its persisted mmCIF structure artifact, ligand-ranking output, and parsed run metrics.

![Boltz-2 predicted protein-ligand complex structure](assets/boltz2-affinity-structure-results.png)

![Boltz-2 affinity and confidence summary metrics](assets/boltz2-affinity-summary-metrics.png)

## Implementation Notes

This workflow contains a batch-specific core model:

- `models/core/src/boltz2_batch_ligand_ranker.py`: CSV parsing, repeated Boltz execution, ranking, flags, top-ligand output selection.
- `models/core/src/boltz2_affinity_predictor.py`: the reused single-ligand Boltz-2 runner.
- `models/visualisation/src/docking_visualisation.py`: batch-aware top-structure and ranked-table visualisation.
