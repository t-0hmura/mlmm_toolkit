# bond-summary

Detect and report covalent bond changes between consecutive molecular structures.

## Synopsis

```bash
mlmm bond-summary -i R.xyz -i P.xyz
mlmm bond-summary -i R.xyz -i TS.xyz -i P.xyz
mlmm bond-summary -i R.pdb -i IM1.pdb -i IM2.pdb -i P.pdb
```

## Description

`bond-summary` compares consecutive pairs of input structures and reports
bonds that are formed or broken. For *N* input files it produces *N − 1*
comparison blocks (A→B, B→C, …).

Bond perception uses element-specific covalent radii with configurable
tolerances. Distances are reported in Ångström.

Supported formats: **XYZ**, **PDB**, **GJF** (auto-detected by extension).

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input FILE` | Input structure file (repeat for each file, ≥ 2 required) | — |
| `--device TEXT` | Compute device (`cpu`, `cuda`) | `cpu` |
| `--bond-factor FLOAT` | Scaling factor for covalent radii sum | `1.20` |
| `--one-based / --zero-based` | Atom index convention in output | `--one-based` |

## Examples

### Two-structure comparison

```bash
mlmm bond-summary -i 1.R.xyz -i 3.P.xyz
```

Output:
```
============================================================
  1.R.xyz  →  3.P.xyz
============================================================
Bond formed (2):
  - O14-H106 : 1.502 Å --> 1.011 Å
  - P95-O107 : 3.477 Å --> 1.523 Å
Bond broken (2):
  - P95-O97 : 1.585 Å --> 3.270 Å
  - H106-O107 : 1.034 Å --> 1.673 Å
```

### Multi-structure (reaction pathway)

```bash
mlmm bond-summary -i 1.R.xyz -i 3.IM1.xyz -i 5.IM2.xyz -i 7.P.xyz
```

Produces three comparison blocks: R→IM1, IM1→IM2, IM2→P.

## Notes

- All input structures must have **identical atom counts and element ordering**.
- Bond detection uses the same algorithm as the internal `bond_changes` module
  used by the `all` workflow for IRC endpoint validation.
- To adjust sensitivity to borderline bonds (e.g., metal coordination at 2.0–2.4 Å),
  increase `--bond-factor` (e.g., `1.30`).
