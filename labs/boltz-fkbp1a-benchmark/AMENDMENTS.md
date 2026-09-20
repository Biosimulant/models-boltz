# Transparent post-registration amendments

## 2026-09-19: adapter transport and provenance repairs

The first C attempt failed before inference because a typed record's payload envelope was not unwrapped. Revision 3 uses the SDK's `unwrap_payload` helper. A real BioWorld transport regression with an explicitly simulated CLI failure verifies propagation without creating benchmark predictions.

A second C attempt failed at provenance collection because the managed Python interpreter had no pip module. Revision 4 uses standard-library package metadata and nonfatal bounded hardware metadata collection. No durable scientific output was available from that attempt, so none is counted.

## 2026-09-19, after verified run 308fe4ee-b9a1-4549-913f-83c8688b705b: T4 attention fallback

The managed service chose a Tesla T4. The verified 20,612-byte result artifact (SHA-256 29794db4621b15efef1fcebc46219151ecc6cad141820279861abf4dc2e4688c) reports Boltz failure in trifast/Triton kernel compilation (`IndexError: map::at`), zero completed predictions, and exactly the baseline's model-checkpoint hashes. The gateway exposed no per-run GPU-type selector; it selects the available default GPU size.

Revision 5 adds upstream Boltz 2.0.2's documented `--no_trifast` flag to C, using its PyTorch attention implementation. This is a hardware-compatibility amendment, not score-based model tuning. Baseline A used the default accelerated attention kernel. NIM's internal kernel choices remain undisclosed. The fallback therefore adds a disclosed numerical/performance confound to A-versus-C comparison; no exact numerical parity or causal software speedup is claimed.

Scientific identities, query-only MSA, Boltz version and weights, sampling counts, seeds, potentials, score interpretation, analysis, questions and acceptance criteria are unchanged. The original PROTOCOL.md is preserved. Failed attempts remain in the evidence ledger. No synthetic result, reused prediction or selected best run replaces them.

Source: pinned upstream CLI defines `--no_trifast` and selects the same model's reference attention path: https://github.com/jwohlwend/boltz/blob/v2.0.2/src/boltz/main.py ; https://github.com/jwohlwend/boltz/blob/v2.0.2/docs/prediction.md .
