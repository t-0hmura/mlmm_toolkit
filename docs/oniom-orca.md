# `oniom-orca`

Export an ML/MM system to ORCA QM/MM (`.inp`) using an Amber parm7 topology. This is the ORCA-specific detail page for `oniom-export` (`mlmm oniom-export --mode orca`). It reads topology information from an Amber parm7 file, maps a model-region PDB to the QM atoms, resolves ORCAFF parameters, and writes a single ORCA QM/MM input file ready for an ORCA 6.0 run.

## Examples

Basic export:

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1
```

Set explicit total charge/multiplicity for the full QM+MM system:

```bash
mlmm oniom-export --mode orca --parm real.parm7 -i pocket.pdb --model-pdb ml_region.pdb \
 -o system.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

Use an explicit ORCAFF path and disable auto-conversion:

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

## Outputs

- `system.inp`
- `ORCAFF.prms` (either reused, provided, or auto-generated)

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

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
| `--near FLOAT` | Distance cutoff (Å) used to define ActiveAtoms when layer tags are absent. | `6.0` |
| `--orcaff PATH` | Path to ORCAFF.prms. | _None_ (auto-resolved to `<parm7_stem>.ORCAFF.prms`) |
| `--convert-orcaff/--no-convert-orcaff` | Attempt `orca_mm -convff -AMBER` when ORCAFF is missing. | `True` |

## See Also

- [oniom_gaussian](oniom-gaussian.md) — Gaussian mode guide (`--mode g16`)
- [oniom_export](oniom-export.md) — Export overview and chooser
- [mm_parm](mm-parm.md) — Build Amber topology (`parm7`/`rst7`)
- [define_layer](define-layer.md) — Build/check layer annotations
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- ORCA 6.0 Manual (QM/MM): <https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html>
