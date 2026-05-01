# `mlmm scan2d`

## Purpose

2D distance scan with harmonic restraints. Drives two bonds toward
target distances simultaneously and produces a 2D grid of relaxed
geometries. Useful for mapping concerted-vs-stepwise reaction
surfaces (e.g. SN2 attack + leaving-group departure).

## Synopsis

```bash
mlmm scan2d -i input.pdb \
    -s '[(a1, b1, target_A), (a2, b2, target_B)]' \
    [-l 'RES:Q,...'] \
    [-b uma|orb|mace|aimnet2] [-o ./result_scan2d/]
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
| `-s, --scan-lists` | str | required | Inline Python literal containing **two** tuples (one per axis), or a YAML/JSON spec file. |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_scan2d/` | Output directory |
| `--ref-pdb` | path | none | Residue context for XYZ/GJF |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

The two tuples in the `-s` literal define the two scan axes; both bonds
are driven simultaneously, generating the grid.

## Examples

```bash
mlmm scan2d -i 1.R.pdb -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60), ("GPP 321 H11","GLU 186 OE2",0.90)]' \
    -b uma -o result_scan2d
```

## Output

```
result_scan2d/
├── result.json
├── grid_NN_MM/                 # per grid-point relaxed geometries
│   └── final.xyz
├── surface.csv                 # 2D energy surface (axis_1, axis_2, energy)
└── scan2d.log
```

`result.json` stores grid metadata and energy values; `surface.csv` is
ready for downstream contour plotting.

## Caveats

- `-s` literal must contain **exactly two** tuples for `scan2d`.
- Output volume scales as `n_steps1 × n_steps2`. Keep grids modest;
  10×10 = 100 single-point optimizations.
- Atom-spec syntax is identical to `scan.md`.

## See also

- `scan.md`, `scan3d.md` — 1D / 3D analogs.
- `all-scan-list.md` — sequential scans (different from coupled grid).
- Defaults: `import mlmm.defaults as d; print(d.BIAS_KW, d.OUT_DIR_SCAN2D)`