# Quickstart: `mlmm all`

## Goal

Run the end-to-end ML/MM ONIOM workflow once from a reactant and product PDB pair: ML-region selection → MM topology/layer preparation → MEP search (single-pass `path-opt` by default) → optional TS optimization, thermochemical correction, IRC, and DFT single-point.

## Prerequisites

- **Two complex PDBs** with matching atom order and identical residues (R/P; only coordinates differ — in PyMOL tick *Original atom order* on export).
- **Explicit hydrogens present** (`mlmm all` does not auto-protonate); match `-l RES:CHARGE` to the H count in the file (e.g. SAM 23 H = `SAM:1`, 22 H = `SAM:0`; a mismatch triggers an antechamber sqm odd-electron failure).
- **GPU recommended.** Default backend `uma` (also `orb` / `mace` / `aimnet2` via `-b`).

## Minimal command

Download `r_complex.pdb` and `p_complex.pdb` from the [prepared example](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/toy_system) into one working directory. Run this MEP-only command from that directory:

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 --out-dir ./result_all
```

With TS optimization, thermochemistry, and DFT single-point:

```bash
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -r 6.0 \
  --ligand-charge 'PRE:0' -q -1 -m 1 \
  --tsopt --thermo --dft --out-dir ./result_all
```

(`-c` = ligand residue list for ML-region centring; `-r` = ML-region radius in Å; `-q`/`-m` = ML-region charge / multiplicity.)

## Result

Check [result status and reasons](json-output.md#execution-and-requested-stage-completion) in `summary.json` before interpreting the segment energies and bond changes in `summary.log`. The merged path is `mep.pdb` (plus `mep.cif` for bridged input) with `energy_diagram_MEP.png` at the output root. See [all](all.md) for the full output tree and per-segment deliverables, and [output-layout](output-layout.md) for the filename reference.

## Next step

- Single-structure scan route: [Quickstart: `mlmm scan` with `-s` (YAML spec)](quickstart-scan-spec.md)
- TS validation route: [Quickstart: `mlmm tsopt` → `mlmm freq`](quickstart-tsopt-freq.md)
- Common failures: [recipes-common-errors](recipes-common-errors.md) · [troubleshooting](troubleshooting.md)
