# bond-summary

Detect and report covalent bond changes between consecutive molecular structures. `bond-summary` compares consecutive pairs of input structures and reports bonds that are formed or broken, which is useful for inspecting bond formation / breaking between reactant, intermediate, and product structures along a reaction pathway, validating IRC endpoint connectivity, or checking borderline coordination (e.g. metal coordination at 2.0–2.4 Å) by tuning `--bond-factor`. For *N* input files it produces *N − 1* comparison blocks (A→B, B→C, …). Bond perception uses element-specific covalent radii with configurable tolerances, and distances are reported in Ångström.

## Examples

Compare a reactant and a product:

```bash
mlmm bond-summary -i 1.R.xyz -i 3.P.xyz
```

Compare a multi-structure pathway (each consecutive pair is reported):

```bash
mlmm bond-summary -i 1.R.xyz -i 3.IM1.xyz -i 5.IM2.xyz -i 7.P.xyz
```

Bare positional files are also accepted (or mix `-i` with positionals):

```bash
mlmm bond-summary A.xyz B.xyz
mlmm bond-summary -i R.xyz TS.xyz P.xyz
```

Supported formats: **XYZ**, **PDB**, **GJF** (auto-detected by extension).

## Outputs

`bond-summary` writes no files. It prints one comparison block per consecutive pair to stdout, each listing the bonds formed and broken with their before/after distances in Ångström. Redirect stdout to persist the report; with `--json` the report is printed as machine-readable JSON instead.

```text
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

A multi-structure run such as `-i 1.R.xyz -i 3.IM1.xyz -i 5.IM2.xyz -i 7.P.xyz` produces three comparison blocks: R→IM1, IM1→IM2, IM2→P.

## CLI options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --input FILE` | Input structure file (repeat for each file, ≥ 2 required); consecutive pairs are compared in order. Bare positional files are also accepted, and may be mixed with `-i`. | — |
| `--device TEXT` | Compute device (`cpu`, `cuda`) | `cpu` |
| `--bond-factor FLOAT` | Scaling factor for covalent radii sum | `1.20` |
| `--one-based / --zero-based` | Atom index convention in output | `--one-based` |
| `--json / --no-json` | Print machine-readable JSON to stdout instead of the text report (no file is written; redirect stdout to persist it). | `--no-json` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes

- All input structures must have **identical atom counts and element ordering**.
- Bond detection uses the same algorithm as the internal `bond_changes` module used by the `all` workflow for IRC endpoint validation.
- To make bond detection more permissive for borderline bonds, increase `--bond-factor` (e.g., `1.30`).

## See Also

- [`all`](all.md)
- [`energy-diagram`](energy-diagram.md)
