# Boltz benchmark model compatibility contract v1

This contract supplements Biosimulant model-compatibility standard v0. Its batch-record schemas are local benchmark interfaces, not new globally registered profiles. Declaring compatibility does not prove biological correctness.

The authoritative inference implementation is upstream Boltz v2.0.2, git 3755c35ed6b9cd764435ec9379bc03a0b72b5e22. The Hub adapter source is demi/boltz-boltz2-affinity-predictor@1.1.1, package SHA-256 ef6c5b3838f960065c9544a1d75e08275cffed7905b13775b4c62cf25cdd8233. The new Lab retains source provenance but uses a finite batch adapter to expose seeds and raw structural artifacts.

| Boundary | Required representation and scientific meaning |
|---|---|
| Sequence | Unmodified UniProt P62942 sequence, 108 standard amino acids; protein chain A. No identification solely by filename. |
| MSA | A3M query exactly equals sequence. Query-only, zero evolutionary homologs. This limits quality. |
| Ligand | One isomeric SMILES for chain B; preserve stereochemical annotations, PubChem CID and canonical isomeric SMILES. RDKit must parse it; emitted CIF heavy-atom count must agree. This count check does not establish stereochemical coordinate accuracy. |
| Native score | Finite real affinity_pred_value. Conventional log10 concentration scale in micromolar; lower means stronger predicted score. Mixed Kd/Ki/IC50 training labels. Never measured IC50, calibrated Kd or free energy. Related scalar profile: boltz.log10-ic50-micromolar/v1. |
| Derived p-scale | Explicit adapter p = 6 − native score. Algebraic scale conversion only; it does not change the underlying mixed-label prediction into an experimental pIC50. |
| Binder probability | Finite [0,1], dimensionless; boltz.binding-probability/v1. Not an empirical success rate or a confidence interval. |
| Structure | Nonempty mmCIF, Cartesian coordinates in angstrom, protein A and ligand B present. Profile protein-ligand.complex-structure-mmcif/v1. Compare only after sequence-aware rigid alignment. |
| Confidence | Preserve provider field names and scale; no equivalence between NIM aggregate confidence and native subfields without documented mapping. |
| Descriptors | RDKit 2025.3.3 average molecular mass in g/mol; heavy-atom count dimensionless. Do not confuse mass with affinity correction. |
| Geometry | Protein–ligand non-hydrogen atom pair distances ≤4 angstrom; multiple atom pairs per residue possible. Count is not hydrogen-bond energy or ligand pose accuracy. |
| Statistics | Three independently executed replicates; mean, sample SD (n−1), min/max; missing results remain missing. |

The finite graph is CandidateModel → PredictionModel → AnalysisModel. Each model declares typed Python SignalSpecs and matching model.yaml record ports. Records are heterogeneous envelopes with explicit field units, not interchangeable scalar quantities. A record containing coordinates cannot be connected to a concentration input merely because both serialize as JSON. Runtime validation must inspect the fields and molecular identities in addition to schema equality.

All modules execute once before the temporal run in topological order. The duration is bookkeeping, not physical time or an MD trajectory. Every dependency is exactly pinned. The native Boltz call includes random seed and both structure/affinity sampling controls. NIM controls are separately audited; an unavailable seed, checkpoint hash, version or profile must be marked uncontrolled.

Reject missing input, invalid SMILES, query mismatch, nonfinite outputs, out-of-range probability, missing structures, unexpected ligand atom counts and changed reference checksum. Abort a batch on inference failure to conserve resources; report remaining candidates as not run. Never replace failures with fixture predictions.

The C arm is only complete after this Lab is saved as a remote immutable workspace revision, run through MCP, and its nonempty workspace-results bytes independently match server size/SHA-256. A local archive or READY Passport alone does not meet that condition.

Sources: https://github.com/jwohlwend/boltz/tree/v2.0.2 ; https://github.com/jwohlwend/boltz/blob/v2.0.2/docs/prediction.md ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12262699/ ; https://www.rcsb.org/structure/1FKJ ; https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit
