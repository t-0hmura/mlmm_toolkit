# Gaussian ONIOM Mode (`oniom-export --mode g16`)

## Overview

> **Summary:** Export an ML/MM system to Gaussian ONIOM (`.com`/`.gjf`) using an Amber parm7 topology.

### At a glance
- **Use when:** You want to run Gaussian ONIOM on the same system prepared for ML/MM.
- **Method:** Reads Amber parm7 topology via ParmEd; maps model PDB to define QM layer; detects QM/MM boundaries and annotates link atoms.
- **Outputs:** Gaussian `.com` or `.gjf` with method, coordinates, layer flags, and connectivity.
- **Defaults:** `--method "B3LYP/6-31G(d,p)"`; `--nproc 8`; `--mem 16GB`; `--near 6.0`.
- **Next step:** Run the exported input in Gaussian.

## Minimal example

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

## Output checklist

- `system.com` or `system.gjf`
- Console summary of detected QM atoms and boundary handling

## Common examples

1. Basic export with explicit method.

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --method "B3LYP/6-31G(d,p)"
```

2. Disable element-sequence validation when atom order is already trusted.

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.xyz --model-pdb ml_region.pdb \
 -o system.gjf --no-element-check
```

3. Tune optimization environment metadata.

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --nproc 16 --mem 32GB --near 5.0
```

## Workflow

Gaussian mode (`mlmm oniom-export --mode g16`) reads topology data from `--parm7` (via ParmEd), optional coordinates from `-i/--input`, then writes a Gaussian ONIOM input file.

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

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Requires `parmed`.
- Element validation is best-effort and assumes invariant atom ordering.
- `--near` controls which MM atoms are marked movable around the QM region.

---

## See Also

- [oniom_orca](oniom_orca.md) -- ORCA mode guide (`--mode orca`)
- [oniom_export](oniom_export.md) -- Export overview and chooser
- [mm_parm](mm_parm.md) -- Build Amber topology (`parm7`/`rst7`)
- [define_layer](define_layer.md) -- Build/check layer annotations
