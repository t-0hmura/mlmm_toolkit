# `oniom-export`

Export an Amber-topology ML/MM system into an external QM/MM input file — Gaussian ONIOM (`--mode g16`, with link-atom annotations) or ORCA QM/MM (`--mode orca`, with ORCAFF handling). It combines an Amber `parm7` topology and an MLMM layered PDB into a single ready-to-run input file. The layered PDB B-factors define the movable/frozen partition; `--model-pdb`, when supplied, overrides its QM-region membership.

Both export modes require a CMAP-free `parm7`. Gaussian ONIOM cannot
represent these terms faithfully, and ORCA's MM engine does not apply them;
the exporter therefore fails before writing when the topology contains CMAP.
This is an export-format limitation—normal mlmm calculations may keep CMAP
enabled in both MM layers. See [CMAP-free preparation](mm-parm.md#cmap-free-topology-for-oniom-export).

## Examples

```bash
# Gaussian ONIOM input
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 -q 0 -m 1
```

```bash
# ORCA QM/MM input (mode inferred from the .inp suffix)
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.inp -q 0 -m 1
```

```bash
# Gaussian input with a custom method/basis and resources
mlmm oniom-export --parm real.parm7 -i pocket_layered.pdb --model-pdb ml.pdb \
 -o out.gjf --mode g16 --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB -q 0 -m 1
```

## Workflow

1. **Topology + layers** — read the `parm7` and the layered PDB passed to `-i` (atom order must match the topology; `--element-check` validates the element sequence). B-factors near 0/10/20 define ML, movable MM, and frozen MM atoms.
2. **QM region** — use the B-factor ML layer unless `--model-pdb` explicitly defines the QM atoms. Movable/frozen membership always comes from the layered PDB.
3. **QM/MM boundary** — Gaussian uses `--link-atom-method scaled` (the default Morokuma/Dapprich g-factor) or `fixed` (1.09/1.01 Å) to place link H atoms. ORCA uses `QMAtoms`/`ORCAFF` for capping; exported link coordinates are diagnostic comments only.
4. **Write** — emit the target-format input file at `-o`. ORCA mode additionally resolves `ORCAFF.prms`. With `--convert-orcaff`, conversion is attempted through `orca_mm -convff -AMBER`; if conversion is disabled or unavailable, the `.inp` is still written and reports the parameter file that must be supplied before ORCA is run.

## Outputs

- `<output>.{gjf,com}` (g16) or `<output>.inp` (ORCA) — the QM/MM input file
- ORCA mode references `<parm7_stem>.ORCAFF.prms`; it reuses an existing file or creates one only when automatic conversion is enabled and available

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `--parm PATH` | Amber parm7 topology file. | Required |
| `-i, --input PATH` | MLMM layered PDB; atom order must match the parm7 and B-factors define movable/frozen atoms. | Required |
| `--model-pdb PATH` | PDB defining the QM-region atoms. | _None_ |
| `-o, --output PATH` | Output file path (`.gjf` / `.com` for g16, `.inp` for ORCA). | Required |
| `--mode [g16\|orca]` | Export mode; inferred from the `-o` suffix when omitted. | _inferred_ |
| `--method TEXT` | QM method and basis set. | mode-dependent |
| `-q, --charge INT` | Charge of the QM region. | Required |
| `-m, --multiplicity INT` | Multiplicity of the QM region. | `1` |
| `--nproc INT` | Number of processors. | `8` |
| `--mem TEXT` | Memory allocation (g16 mode). | `16GB` |
| `--total-charge INT` / `--total-mult INT` | Total charge / multiplicity of the full QM+MM system (ORCA `Charge_Total` / `Mult_Total`). | topology-derived / same as `--multiplicity` |
| `--orcaff PATH` | Path to `ORCAFF.prms` (ORCA mode). If omitted, a derived path is referenced and automatic creation is attempted conditionally. | _None_ |
| `--convert-orcaff / --no-convert-orcaff` | Auto-convert a missing `ORCAFF.prms` via `orca_mm -convff -AMBER` (ORCA mode). | `True` |
| `--element-check / --no-element-check` | Validate the `--input` element sequence against the parm7 topology. | `True` |
| `--link-atom-method [scaled\|fixed]` | Gaussian link-H placement; ORCA records the corresponding coordinates only as diagnostics and creates caps from `QMAtoms`/`ORCAFF`. | `scaled` |

`mlmm oniom-export --help` shows core options; `mlmm oniom-export --help-advanced` shows the full list.

## Notes

- Mode selection: `--mode` is highest priority. If `--mode` is omitted, the mode is inferred from `-o`:
  - `.gjf` / `.com` → `g16`
  - `.inp` → `orca`
- If `--mode` is omitted and the `-o` suffix is unknown, the command fails.
- For PDB/ENT input, the exported file embeds
  `MLMM_REF_PDB_ORDER_V1_SHA256=<digest>`. Coordinates, occupancy, and B-factor
  are excluded, while fixed atom/residue/chain/insertion/element identity is
  covered. `oniom-import --ref-pdb` verifies this marker before positional
  metadata restoration.

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [oniom-gaussian](oniom-gaussian.md) — Gaussian-mode details (`--mode g16`)
- [oniom-orca](oniom-orca.md) — ORCA-mode details (`--mode orca`)
- [oniom-import](oniom-import.md) — Reconstruct XYZ/layered PDB from ONIOM inputs
- [mm-parm](mm-parm.md) — Build Amber topology
- [define-layer](define-layer.md) — Build/check layer annotations
