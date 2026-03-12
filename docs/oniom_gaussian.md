# Gaussian ONIOM Mode (`oniom-export --mode g16`)

## Overview

> **Summary:** Export an ML/MM system to Gaussian ONIOM (`.com`/`.gjf`) using an Amber parm7 topology.

## Minimal example

```bash
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com -q 0 -m 1
```

## Output checklist

- `system.com` or `system.gjf`
- Console summary of detected QM atoms and boundary handling

## Common examples

1. Basic export with explicit method.

```bash
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --method "wB97XD/def2-TZVPD"
```

2. Disable element-sequence validation when atom order is already trusted.

```bash
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket.xyz --model-pdb ml_region.pdb \
 -o system.gjf --no-element-check
```

3. Tune optimization environment metadata.

```bash
mlmm oniom-export --mode g16 --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.com --nproc 16 --mem 32GB --near 5.0
```

## Workflow

Gaussian mode (`mlmm oniom-export --mode g16`) reads topology data from `--parm` (via ParmEd), optional coordinates from `-i/--input`, then writes a Gaussian ONIOM input file.

1. Load atom/bond/charge data from parm7.
2. Optionally load coordinates and run element-order validation (`--element-check`).
3. If `--model-pdb` is provided, map those atoms to define the QM layer.
4. Detect QM/MM boundaries and annotate link atoms.
5. Write Gaussian input with method, `%nproc`, `%mem`, coordinates, layer flags, and connectivity.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (PDB/XYZ); atom order must match parm7. | _None_ |
| `--element-check / --no-element-check` | Validate element sequence between input and parm7. | `True` |
| `--model-pdb PATH` | PDB file defining QM region atoms. | _None_ |
| `-o, --output PATH` | Output Gaussian input file (`.com` or `.gjf`). | Required |
| `--method TEXT` | QM method and basis set. | `wB97XD/def2-TZVPD` |
| `-q, --charge INT` | Charge of QM region. | Required |
| `-m, --multiplicity INT` | Multiplicity of QM region. | `1` |
| `--near FLOAT` | Distance cutoff for movable atoms (Å). | `6.0` |
| `--nproc INT` | Number of processors. | `8` |
| `--mem TEXT` | Memory allocation. | `16GB` |

## See Also

- [oniom_orca](oniom_orca.md) -- ORCA mode guide (`--mode orca`)
- [oniom_export](oniom_export.md) -- Export overview and chooser
- [mm_parm](mm_parm.md) -- Build Amber topology (`parm7`/`rst7`)
- [define_layer](define_layer.md) -- Build/check layer annotations
