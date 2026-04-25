# `mlmm scan`

## Purpose

1D bond-length-driven scan with staged harmonic restraints and inter-stage
relaxation. Drive one or more bonds toward target distances, producing a
trajectory you can reuse to seed `path-search`. For most workflows prefer
`mlmm all --scan-lists` (see `all-scan-list.md`); standalone
`scan` is for one-off exploration.

## Synopsis

```bash
mlmm scan -i input.pdb \
    -s '[(idx_a, idx_b, target_A), ...]' \
    [-l 'RES:Q,...'] [-q / -m] \
    [-b uma|orb|mace|aimnet2] [-o ./result_scan/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

Beyond the cluster-style flags below (inherited from `pdb2reaction`),
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); overrides `--model-pdb` |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Reactant `.pdb` / `.xyz` / `.gjf` |
| `-s, --scan-lists` | str | required | Inline Python literal `'[(a,b,target),...]'`, or YAML/JSON spec path. Repeat to add stages. |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | `none` | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF inputs |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

The tuple grammar in `-s` accepts atom-index ints (`(1, 5, 1.4)`) or atom
specs (`("CS1 SAM 320", "C7 GPP 321", 1.60)`). Multiple `-s` flags
chain stages sequentially: each stage starts from the previous stage's
final geometry.

## Examples

### Single stage by atom name

```bash
mlmm scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
    -b uma -o result_scan
```

### Two sequential stages

```bash
mlmm scan -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' \
    -s '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
    -b uma -o result_scan_staged
```

## Output

```
result_scan/
├── result.json
├── stage_01/                # per-stage relaxed snapshots
│   └── scan_*.xyz
├── stage_02/
├── mep.xyz                  # stitched scan trajectory
└── scan.log
```

`result.json` lists per-stage status, target distances, final energies,
and the stitched MEP path. Plot with `trj2fig.md`.

## Caveats

- `-s` is Python literal-eval. Quote with single quotes outside,
  double quotes inside. Atom-name strings use `"NAME RESNAME RESID"`
  with single spaces.
- Stage *k+1* starts from stage *k*'s final geometry; a diverged
  stage poisons downstream stages.
- For coupled multi-bond drives in one stage, put multiple tuples in
  one `-s` argument: `'[(a,b,1.6),(c,d,3.0)]'`.

## See also

- `scan2d.md`, `scan3d.md` — higher-dim analogs.
- `all-scan-list.md` — wraps `scan` inside the full pipeline.
- Defaults: `import mlmm.defaults as d; print(d.BIAS_KW, d.BOND_KW, d.OUT_DIR_SCAN)`