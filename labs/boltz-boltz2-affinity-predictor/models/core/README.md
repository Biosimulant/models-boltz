# Boltz-2 affinity model

This `biosimulant.BioModule` wraps `boltz predict` from
`boltz[cuda]==2.0.2`. It accepts one protein, one molecular SMILES value, and
either a matching A3M file or explicit MSA-server mode.

The precise scientific ports are:

- `protein_sequence` → `protein.sequence/v1`
- `ligand_smiles` → `chemical.smiles/v1`
- `msa_path` → `protein.multiple-sequence-alignment/v1`
- `binding_probability` → `boltz.binding-probability/v1`
- `affinity_log10_ic50_micromolar` → `boltz.log10-ic50-micromolar/v1`
- `predicted_structure` → `protein-ligand.complex-structure-mmcif/v1`

`run_options`, the three aggregate summaries, and `run_metadata` are
operational or multi-field records and intentionally have no profile.

Boltz's `affinity_pred_value` is `log10(IC50)` with IC50 expressed in
micromolar. It is not pIC50. The model accepts mmCIF output only, checks an A3M
query against the protein before compute, validates atomic values, and does not
emit placeholder atomic outputs after failure.

Managed mode creates a cached local Boltz environment and retains the existing
single cache-repair retry. Remote mode uses the GPU runtime image and stores the
Boltz cache under the remote execution mount. See the Lab-level README for the
complete interface, defaults, validation commands, and scientific limits.
