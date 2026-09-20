# Working evaluation requirements, v1

Status: implementation authorized by the user benchmark request on 2026-09-19; no separate formal scientific signoff claimed.

# Preregistered Boltz three-arm pilot

Date: 2026-09-19. Written before inference results are available.

This is an execution, composition and reproducibility pilot, not a blinded evaluation of independent AI agents, a training benchmark, or an estimate of clinical/experimental efficacy. The same author prepares all adapters; setup knowledge is shared. No arm receives fabricated model outputs. Test fixtures are never counted as predictions.

## The three questions and expected answers

1. **Create and run a reusable structure/affinity prediction model.** Starting from human FKBP1A (UniProt P62942) and tacrolimus (PubChem CID 445643), create a source-pinned Boltz-2 inference adapter, declare its input/output compatibility contract, and run it. Return the complex mmCIF, confidence, native affinity score and binder probability, with provenance and runtime. **Expected answer:** one real, nonempty predicted protein–ligand structure; finite score; probability in [0,1]; correctly identified chains and ligand; no preassigned numerical affinity. A model here means an inference component around pretrained weights, not a newly trained foundation model.

2. **Compose models across sequence, chemical graph and 3-D structure modalities.** Compose a molecular-graph descriptor component (RDKit), the same Boltz protein/ligand predictor, a geometric contact analysis component, and a ranking component for tacrolimus, sirolimus (CID 5284616), and ascomycin (CID 5282071). Run the composition, preserving stereochemistry and molecular identity. **Expected answer:** exactly three candidate rows with molecular mass (g/mol), heavy-atom count, predicted affinity score, probability, confidence, and protein–ligand heavy-atom contacts within 4 Å. Rank native affinity ascending, retaining probabilities independently. The experimentally best ligand is unknown here; do not invent a ground-truth ranking. A geometric contact is not a demonstrated hydrogen bond. These are heterogeneous data modalities, not a validated multiscale kinetic model.

3. **Run a reproducibility and structural-consistency experiment.** Repeat FKBP1A–tacrolimus prediction three times, using seeds 101, 202 and 303 wherever the runtime exposes them. Compose the predictor with sequence-aware alignment to experimental complex PDB 1FKJ and a statistical summary component. Return the three raw outputs, mean/sample SD/range of predicted affinity and probability, and Cα RMSD after superposition of exactly matched residues. **Expected answer:** three independently executed predictions, calculated statistics (SD undefined if fewer than two succeed), explicit alignment coverage, and an RMSD calculated from coordinates. No threshold is asserted for biological validity. Report uncontrolled seeds for a hosted runtime. PDB 1FKJ predates training; this is retrospective structural consistency, not held-out generalization. Ligand RMSD is not reported without verified atom mapping/symmetry handling.

## Arms

| Arm | Execution |
|---|---|
| A | Upstream Boltz + ordinary Python, without BioNeMo or BioSimulant |
| B | NVIDIA BioNeMo Agent Toolkit Boltz-2 NIM workflow, without BioSimulant |
| C | Boltz wrapped as legitimate BioSimulant Labs, saved remotely and executed through authenticated MCP, without BioNeMo |

The BioNeMo skill is installed only in arm B. BioNeMo Inference Runtime cannot silently replace NIM: its published pipeline does not support affinity. If credentials/hardware are unavailable, B is blocked, not assigned synthetic scores or treated as a model failure.

## Shared controls

- Identical downloaded sequence and isomeric SMILES, source bytes and SHA-256 recorded before inference. Full P62942 sequence used; no sequence editing to improve a result.
- A query-only A3M is explicit and identical in every arm; no live MSA service or template input. This deliberately controls the input and limits prediction quality.
- 3 recycling steps; 200 structure diffusion steps; 1 structure sample; 200 affinity steps; 3 affinity diffusion samples where supported. Potentials off where exposed; molecular-weight correction off. Native pinned version 2.0.2. NIM version/weight opacity must be disclosed; numerical equivalence cannot be assumed.
- A and C use seeds 42 (Q1/Q2), 101/202/303 (Q3). If an arm cannot expose a setting, record it as uncontrolled rather than imply equality.
- Seven inference calls per complete arm: Q1=1, Q2=3, Q3=3. No cross-question reuse counts as a fresh run. All failures and retries remain in the ledger.
- The primary comparison is completion and artifact integrity. Numerical deltas are descriptive; heterogeneous GPUs, versions, weights and service queueing prevent causal speed/accuracy claims.
- Separate cold setup/download, inference wall time, and analysis time. Do not infer monetary charges without provider evidence.
- No model-derived affinity score is wired into a Kd, kinetic rate or mechanistic occupancy model. The log-scale score mixes affinity/activity training labels and is not a measured equilibrium constant.

## Resource limits

One benchmark GPU job at a time. Up to 30 minutes per individual prediction and 90 minutes per arm, with explicit termination and persisted failure records. Cache downloaded weights. Do not allocate a persistent machine merely for an idle service. Stop benchmark-created jobs/services after retrieval. Do not stop unrelated pre-existing Brev jobs. No broad parameter searches, retraining, hidden retries or output selection.

## Acceptance and reporting

Record status for all nine question×arm cells: completed, failed, partial or blocked. Only a completed real execution with raw artifacts can satisfy an inference question. Local wrapper tests and schema validation are separate. Preserve immutable remote workspace revision, run ID, verified result checksum and Passport for C. Scientific limitations remain even if a Passport is READY.

For confidence/affinity, finite/range tests are technical consistency checks, not calibration evidence. Three repeats describe stochastic variability only; do not claim robust statistical significance or experimentally superior drug activity.
