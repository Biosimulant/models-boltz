# BioNeMo Boltz result compatibility contract v1

The Lab consumes captured NVIDIA Boltz-2 requests/responses, not a native Boltz checkpoint.
NVIDIA inference executes in the authorized agent environment. Its API key never enters the Lab.
The source node preserves exact response text, request content, HTTP status, execution ID and hashes.
An unmodified, identified 108-residue human FKBP1A sequence and matching query-only MSA are required.
One protein chain A and one ligand chain B are supported by this demo adapter.
Protein sequence must match exactly. Ligand heavy-element inventory must match RDKit-parsed input SMILES.
Heavy-element inventory is not a connectivity, stereochemistry or pose-accuracy check.
Coordinates are angstrom. The native affinity score retains its mixed-label log-micromolar convention.
Binder probability is dimensionless in [0,1]. Provider affinity_pic50 is preserved separately and never relabelled as measured affinity.
Record ports are mixed-field envelopes. Their schema does not itself establish scientific compatibility.
The acceptance component inspects scientific fields and emits executed checks. The ranking component
must publish no ranking if acceptance fails, including partial batches. A process exit is not acceptance.
Hashes establish captured-byte identity, not an NVIDIA digital signature or independent provider attestation.
Replay and reused results are explicitly labelled, with original execution IDs. Neither counts as new inference.
All modules run once before temporal bookkeeping. No biological temporal simulation is implied.
This is research workflow verification, not experimental validation or an independent-agent quality guarantee.

Upstream model: https://github.com/jwohlwend/boltz/tree/v2.0.2
BioNeMo skill source: NVIDIA-BioNeMo/bionemo-agent-toolkit commit 0e67a612e4045f007e38fa77adc8f3ebfc5616b6.
