# `oniom-export`

Export an Amber-topology ML/MM system into an external QM/MM input file — Gaussian ONIOM (`--mode g16`, with link-atom annotations) or ORCA QM/MM (`--mode orca`, with ORCAFF handling). Use it to turn an Amber `parm7` topology plus a coordinate file and the ML-region (QM) definition into a single ready-to-run QM/MM input file: the QM region is taken from `--model-pdb`, and the surrounding MM environment is emitted in the target program's native format with link-atom annotations at the QM/MM cut.

## Examples

```bash
# Gaussian ONIOM input
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 -q 0 -m 1
```

```bash
# ORCA QM/MM input (mode inferred from the .inp suffix)
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.inp -q 0 -m 1
```

```bash
# Gaussian input with a custom method/basis and resources
mlmm oniom-export --parm real.parm7 -i pocket.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB -q 0 -m 1
```

## Workflow

1. **Topology + coordinates** — read the `parm7` and the `-i` coordinate file (atom order must match the topology; `--element-check` validates the element sequence).
2. **QM region** — `--model-pdb` defines the QM (ML-region) atoms; `--near` sets the movable/active MM cutoff (Å).
3. **Link atoms** — placed at each severed QM/MM bond. `--link-atom-method scaled` (default) uses the Morokuma/Dapprich g-factor (matches the `MLMMCore` runtime); `fixed` uses fixed 1.09/1.01 Å bond lengths.
4. **Write** — emit the target-format input file at `-o`. ORCA mode additionally resolves `ORCAFF.prms` (auto-converting from Amber via `orca_mm -convff -AMBER` when `--convert-orcaff` is on).

## Outputs

- `<output>.{gjf,com}` (g16) or `<output>.inp` (ORCA) — the QM/MM input file
- ORCA mode also reads/creates `<parm7_stem>.ORCAFF.prms` in the output directory (force-field parameters)

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | Coordinate file (`.pdb` / `.xyz`) for the current structure. | _None_ |
| `--model-pdb PATH` | PDB defining the QM-region atoms. | _None_ |
| `-o, --output PATH` | Output file path (`.gjf` / `.com` for g16, `.inp` for ORCA). | Required |
| `--mode [g16\|orca]` | Export mode; inferred from the `-o` suffix when omitted. | _inferred_ |
| `--method TEXT` | QM method and basis set. | mode-dependent |
| `-q, --charge INT` | Charge of the QM region. | Required |
| `-m, --multiplicity INT` | Multiplicity of the QM region. | `1` |
| `--near FLOAT` | Distance cutoff (Å) for movable/active MM atoms. | `6.0` |
| `--nproc INT` | Number of processors. | `8` |
| `--mem TEXT` | Memory allocation (g16 mode). | `16GB` |
| `--total-charge INT` / `--total-mult INT` | Total charge / multiplicity of the full QM+MM system (ORCA `Charge_Total` / `Mult_Total`). | _None_ |
| `--orcaff PATH` | Path to `ORCAFF.prms` (ORCA mode); created in the output directory if omitted. | _None_ |
| `--convert-orcaff / --no-convert-orcaff` | Auto-convert a missing `ORCAFF.prms` via `orca_mm -convff -AMBER` (ORCA mode). | `True` |
| `--element-check / --no-element-check` | Validate the `--input` element sequence against the parm7 topology. | `True` |
| `--link-atom-method [scaled\|fixed]` | Link-H placement: `scaled` (g-factor, matches runtime) or `fixed` (1.09/1.01 Å). | `scaled` |

`mlmm oniom-export --help` shows core options; `mlmm oniom-export --help-advanced` shows the full list.

## Notes

- Mode selection: `--mode` is highest priority. If `--mode` is omitted, the mode is inferred from `-o`:
  - `.gjf` / `.com` → `g16`
  - `.inp` → `orca`
- If `--mode` is omitted and the `-o` suffix is unknown, the command fails.

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [oniom_gaussian](oniom-gaussian.md) — Gaussian-mode details (`--mode g16`)
- [oniom_orca](oniom-orca.md) — ORCA-mode details (`--mode orca`)
- [oniom_import](oniom-import.md) — Reconstruct XYZ/layered PDB from ONIOM inputs
- [mm_parm](mm-parm.md) — Build Amber topology
- [define_layer](define-layer.md) — Build/check layer annotations
