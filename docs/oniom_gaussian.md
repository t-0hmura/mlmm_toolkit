# `oniom-gaussian`

## Overview

> **Summary:** Export an ML/MM system to Gaussian ONIOM (`.com`/`.gjf`) using an Amber parm7 topology.

### Quick reference
- **Use when:** You want to run Gaussian ONIOM on the same system prepared for ML/MM.
- **Inputs:** `--parm7` (required), optional `-i/--input` (PDB/XYZ), optional `--model-pdb`.
- **Output:** Gaussian ONIOM input with layer assignments and connectivity.
- **Validation:** `--element-check` (default) compares coordinate-file element order against parm7.

## Minimal example

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

## Output checklist

- `system.com` or `system.gjf`
- Console summary of detected QM atoms and boundary handling

## Common examples

1. Basic export with explicit method.

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --method "B3LYP/6-31G(d,p)"
```

2. Disable element-sequence validation when atom order is already trusted.

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.xyz --model-pdb ml_region.pdb \
 -o system.gjf --no-element-check
```

3. Tune optimization environment metadata.

```bash
mlmm oniom-gaussian --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --nproc 16 --mem 32GB --near 5.0
```

## Usage

```bash
mlmm oniom-gaussian --parm7 real.parm7 [-i coords.pdb] [--model-pdb ml_region.pdb] \
 -o output.com [-q CHARGE] [-m MULT] [--method "B3LYP/6-31G(d,p)"] \
 [--near FLOAT] [--nproc INT] [--mem TEXT] \
 [--element-check|--no-element-check]
```

## Description

`oniom-gaussian` reads topology data from `--parm7` (via ParmEd), optional coordinates from `-i/--input`, then writes a Gaussian ONIOM input file.

Workflow:

1. Load atom/bond/charge data from parm7.
2. Optionally load coordinates and run element-order validation (`--element-check`).
3. If `--model-pdb` is provided, map those atoms to define the QM layer.
4. Detect QM/MM boundaries and annotate link atoms.
5. Write Gaussian input with method, `%nproc`, `%mem`, coordinates, layer flags, and connectivity.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `--parm7 PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (PDB/XYZ); atom order must match parm7. | _None_ |
| `--element-check / --no-element-check` | Validate element sequence between input and parm7. | `True` |
| `--model-pdb PATH` | PDB file defining QM region atoms. | _None_ |
| `-o, --output PATH` | Output Gaussian input file (`.com` or `.gjf`). | Required |
| `--method TEXT` | QM method and basis set. | `B3LYP/6-31G(d,p)` |
| `-q, --charge INT` | Charge of QM region. | Required |
| `-m, --multiplicity INT` | Multiplicity of QM region. | `1` |
| `--near FLOAT` | Distance cutoff for movable atoms (angstrom). | `6.0` |
| `--nproc INT` | Number of processors. | `8` |
| `--mem TEXT` | Memory allocation. | `16GB` |

## Outputs

```text
<output>.com / <output>.gjf
```

## Notes

- Requires `parmed`.
- Element validation is best-effort and assumes invariant atom ordering.
- `--near` controls which MM atoms are marked movable around the QM region.

---

## See Also

- [oniom-orca](oniom_orca.md) -- ORCA QM/MM exporter
- [oniom_export](oniom_export.md) -- Export overview and chooser
- [mm_parm](mm_parm.md) -- Build Amber topology (`parm7`/`rst7`)
- [define_layer](define_layer.md) -- Build/check layer annotations
