# BIO-001: NVIDIA-hosted Boltz-2 requirements

Status: implementation requirements drafted from the user's authorized goal.
Scientific acceptance thresholds and independent review remain outstanding.
This document does not record an approval or a qualified release.

## Intended use

Run NVIDIA's listed Boltz-2 hosted service from one Biosimulant Lab for research
prediction of biomolecular complex coordinates and native confidence/affinity
outputs. Preserve protein, DNA, RNA, ligand and constraint inputs. Preserve all
returned samples, native field names and missing values. A batch of predictions
does not represent biological time evolution. Do not substitute upstream Boltz
or replay saved responses as inference.

## Required behavior

- R01: Accept an explicit molecular-system record containing polymers, ligands
  and constraints. Require unique entity IDs so results can be traced to inputs.
  Keep modifications, MSA and templates with their parent polymer. Reject
  malformed identities, conflicting chemical representations and nonfinite JSON.
- R02: Accept only declared inference controls. Record the exact submitted
  request and its SHA-256. Do not silently truncate sequences, choose a pose,
  convert score scales or supply invented scientific data.
- R03: Use the documented HTTPS prediction endpoint and an environment secret
  `NVIDIA_API_KEY`. Never expose the credential in inputs, manifests, errors,
  output records or source control. No configurable arbitrary endpoint.
- R04: Submit once. Poll an accepted provider request using its ID. A timeout
  is not proof that provider execution stopped; retain the ID and never retry
  the POST automatically. Allow explicit resume of that ID with the same input.
- R05: Expose every returned structure as inline mmCIF with its checksum, byte
  size, native name/source and sample index. Inline payloads are portable inside
  the captured result; no worker-local path is an artifact reference.
- R06: Validate parseability, finite coordinates, unique atom IDs and exact chain
  set. Retain native response bytes as UTF-8 text and their digest. Structural
  integrity does not establish pose accuracy or experimental validity.
- R07: Keep ligand affinity, structure quality and operational provenance in
  separate typed outputs. Preserve `affinity_pic50` independently from
  `affinity_pred_value`; their conversion is unresolved in the captured reference.
- R08: Canonical model IO, Python SignalSpec ports and Lab mappings must agree.
  Execute once before a run, with pinned Python dependencies and bounded HTTP
  observation/response sizes. Mock transport tests cannot satisfy real inference.
- R09: Publish privately only after fresh managed inference has durable nonempty
  results, independent checksum/content inspection and reviewed Passport findings.
  Public qualification needs separate approval of the exact immutable release.

## Acceptance evidence

Technical acceptance requires executable malformed-input, secret-handling,
poll/resume, output-integrity, manifest/runtime parity and BioWorld execution
tests. A captured native request/response may be a parser fixture only. A fresh
adapter request must yield traceable outputs; a subsequent actual Hub run is
required to establish managed execution and credential access.

Scientific acceptance remains BLOCKED pending approved held-out fixtures,
training-overlap review, structural metrics and affinity assay-specific metrics,
numerical tolerance rationale and uncertainty/reproducibility requirements.
Do not invent thresholds or equate confidence with observed accuracy. Protein,
nucleic-acid, modified/cyclic and constrained complexes require suitable separate
coverage. Full ligand connectivity/stereochemistry needs independent checking.

## Identity and rights

The Build card lists model v2.2.1; NIM documentation identifies container 1.9.0.
The hosted endpoint does not provide a pinned checkpoint/container digest in the
observed result. These version namespaces are not interchangeable. Execution of
the hosted service does not prove equivalence to an exact self-hosted release.
Hosted trial-service terms, model terms and redistribution permissions need
separate review before a public deployment. No model weights are redistributed.

Sources reviewed 2026-09-20:
- https://build.nvidia.com/mit/boltz2/modelcard
- https://docs.api.nvidia.com/nim/reference/mit-boltz2-infer
- https://docs.nvidia.com/nim/bionemo/boltz2/latest/inference.html
- https://docs.api.nvidia.com/cloud-functions/reference/getfunctioninvocationresult

The existing native reference has provider ID
`60235ebd-8d17-4134-b89d-5397dd2cc53d` and raw response SHA-256
`23f24ca05f300a6dcce349090c64a78353d9c9ff3ed0b342a7f7c81debf04ac1`.
It is not a managed Lab run or scientific benchmark.
