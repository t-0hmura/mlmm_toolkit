# ORCA QM/MM Mode (`oniom-export --mode orca`)

## Overview

> **Summary:** Export an ML/MM system to ORCA QM/MM (`.inp`) using an Amber parm7 topology.

## Minimal example

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

## Output checklist

- `system.inp`
- `ORCAFF.prms` (either reused, provided, or auto-generated)

## Common examples

1. Basic export.

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

2. Set explicit total charge/multiplicity for full QM+MM system.

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

3. Use an explicit ORCAFF path and disable auto conversion.

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp --orcaff ./ORCAFF.prms --no-convert-orcaff
```

## Workflow

ORCA mode (`mlmm oniom-export --mode orca`) reads topology information from `--parm` and writes ORCA QM/MM input.

1. Load atom/bond/charge data from parm7.
2. Optionally load coordinates and validate element ordering (`--element-check`).
3. Map `--model-pdb` atoms to the full system to define QM atoms.
4. Resolve ORCAFF parameters from `--orcaff` or auto-generated default.
5. Write ORCA `.inp` with method, parallel settings, QM set, and total charge/multiplicity fields.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (PDB/XYZ); atom order must match parm7. | _None_ |
| `--element-check / --no-element-check` | Validate element sequence between input and parm7. | `True` |
| `--model-pdb PATH` | PDB file defining QM region atoms. | _None_ |
| `-o, --output PATH` | Output ORCA input file (`.inp`). | Required |
| `--method TEXT` | QM method and basis set. | `B3LYP D3BJ def2-SVP` |
| `-q, --charge INT` | Charge of QM region. | Required |
| `-m, --multiplicity INT` | Multiplicity of QM region. | `1` |
| `--total-charge INT` | Total charge of full QM+MM system (`Charge_Total`). | topology-derived |
| `--total-mult INT` | Total multiplicity of full QM+MM system (`Mult_Total`). | same as `--multiplicity` |
| `--nproc INT` | Number of processors. | `8` |
| `--near FLOAT` | Distance cutoff used to define ActiveAtoms when layer tags are absent. | `6.0` |
| `--orcaff PATH` | Path to ORCAFF.prms. | auto |
| `--convert-orcaff/--no-convert-orcaff` | Attempt `orca_mm -convff -AMBER` when ORCAFF is missing. | `True` |

## See Also

- [oniom_gaussian](oniom_gaussian.md) -- Gaussian mode guide (`--mode g16`)
- [oniom_export](oniom_export.md) -- Export overview and chooser
- [mm_parm](mm_parm.md) -- Build Amber topology (`parm7`/`rst7`)
- [define_layer](define_layer.md) -- Build/check layer annotations
- ORCA 6.0 Manual (QM/MM): <https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html>
