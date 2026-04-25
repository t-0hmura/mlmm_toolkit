# `mlmm mm-parm`

## Purpose

Generate Amber `parm7` (topology) + `rst7` (coordinates) files from a
protein-ligand PDB by driving AmberTools' `tleap` and (when needed)
`antechamber`. The output pair is consumed by every downstream
subcommand that needs MM gradients (`opt`, `tsopt`, `freq`, …).

## Synopsis

```bash
mlmm mm-parm -i complex.pdb \
    [--ligand 'RES1:Q1,RES2:Q2,...'] \
    [--force-field amber14sb|amber99sb-ildn|amber19sb] \
    [--water tip3p|tip4p-ew|none] \
    [-o complex] \
    [--ph 7.0] \
    [--frcmod custom.frcmod ...]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Protein-ligand PDB (with or without B-factor layer labels) |
| `--ligand` | str | none | Per-ligand charges, e.g. `'SAM:1,GPP:-3'`. Required for non-standard residues. |
| `--force-field` | str | `amber14sb` | Protein force field: `amber14sb`, `amber99sb-ildn`, `amber19sb`, … |
| `--water` | str | `tip3p` | Water model; `none` skips solvation |
| `--ph` | float | 7.0 | Used by `pdb4amber` to assign protonation states |
| `--frcmod` | path(s) | none | Extra parameter files (one per non-standard ligand) |
| `--charge-method` | str | `bcc` | Antechamber charge method (`bcc` / `gas` / etc.) for non-std residues |
| `-o, --out-dir` | path | `./mm_parm/` | Output directory |
| `-o-stem` | str | `<input-basename>` | Basename for the produced `.parm7`/`.rst7` |
| `--show-config` / `--dry-run` | flag | off | Print resolved config; useful to inspect tleap input |

## Examples

### Minimal — standard residues only

```bash
mlmm mm-parm -i complex.pdb -o complex
# produces complex.parm7 / complex.rst7 in ./mm_parm/
```

### With ligand needing antechamber

```bash
mlmm mm-parm -i complex.pdb \
    --ligand 'SAM:1,GPP:-3' \
    --force-field amber14sb \
    --water tip3p \
    -o complex
```

`mm-parm` runs `antechamber -nc <Q>` and `parmchk2` automatically for
each `--ligand` entry that is not in the Amber library, producing
`complex.GPP.frcmod` etc. alongside the topology.

### Reuse pre-built ligand parameters

```bash
mlmm mm-parm -i complex.pdb \
    --ligand 'GPP:-3' \
    --frcmod gpp_gaff2.frcmod \
    -o complex
```

## Output

```
mm_parm/
├── complex.parm7              # Amber topology
├── complex.rst7               # Coordinates (and box if water)
├── complex.frcmod             # Combined extra parameters (if any ligand needed antechamber)
├── tleap.in                   # Generated tleap input (for reproducibility)
├── tleap.out                  # tleap log
└── result.json
```

`result.json` keys:

```python
import json
d = json.load(open("mm_parm/result.json"))
print(d["status"])              # "completed" / "error"
print(d["n_atoms"], d["charge"])
print(d["force_field"], d["water_model"])
print(d["ligands"])             # {"SAM": {"charge": 1, "src": "amber"}, "GPP": {"charge": -3, "src": "antechamber"}}
```

## Caveats

- AmberTools must be installed and `tleap` on PATH; see
  `../mlmm-install-backends/ambertools.md`.
- The PDB must be tleap-clean: H-atoms named per Amber conventions,
  no altloc, no missing element column. Run `mlmm fix-altloc` and
  `mlmm add-elem-info` first if you hit "Atom 'XXX' not in residue
  template" errors.
- The B-factor layer encoding (ML / movable-MM / frozen) is **not**
  carried into `parm7`; downstream subcommands pair `parm7` + `rst7`
  with the original PDB via `--ref-pdb` to recover layer info.
- For small MM regions (< 100 atoms), running `mm-parm` may be
  overkill; consider `pdb2reaction` (no MM at all) instead.

## See also

- `../mlmm-structure-io/parm7.md` — what `parm7`/`rst7` actually contain.
- `../mlmm-install-backends/ambertools.md` — install AmberTools.
- `define-layer.md` — assign / refine ML / movable-MM / frozen labels
  on the PDB *after* `mm-parm` (or before, if you want layer-aware
  parameter generation).
- `extract.md` — alternative front-door that combines extraction +
  layer assignment in one shot.
- Defaults: `import mlmm.defaults as d; print(d.MM_PARM_KW)` (if exported)