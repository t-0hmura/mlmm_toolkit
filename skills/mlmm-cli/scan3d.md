# `mlmm scan3d`

## Purpose

3D distance scan with harmonic restraints. Drives three bonds toward
target distances and produces a 3D grid of relaxed geometries. Rare in
practice; usually a sequence of 1D / 2D scans (or `all-scan-list.md`)
captures the chemistry with less compute. Provided for inherently
3D-coupled mechanisms (e.g. proton-coupled electron transfer with two
protons + one redox-donor distance).

## Synopsis

```bash
mlmm scan3d -i input.pdb \
    -s '[(a1,b1,low1,high1), (a2,b2,low2,high2), (a3,b3,low3,high3)]' \
    [--csv surface.csv] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] [-o ./result_scan3d/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); used only when `--model-pdb` is omitted (`--model-pdb` takes precedence) |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required (unless `--csv`) | Reactant `.pdb` / `.xyz` |
| `-s, --scan-lists` | str | required (unless `--csv`) | Python literal with **three** quadruples `(i,j,low,high)` |
| `--csv` | path | none | Skip the scan; load a precomputed `surface.csv` for downstream plotting |
| `-q` / `-l` / `-m` | — | — | Charge / spin |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_scan3d/` | Output directory |
| `--ref-pdb` / `--config` / `--help-advanced` | — | — | Standard |
| `--print-parsed` | — | — | Validate the parsed scan spec without GPU compute |

## Caveats

- Output volume scales as `n_steps1 × n_steps2 × n_steps3`. Keep grid
  sizes modest (5×5×5 = 125 points already).
- For inherently 1D or 2D mechanisms, `scan.md` / `scan2d.md` are
  drastically cheaper.
- `--csv` is the re-plotting entry point: rerun energy diagrams or
  contour plots without redoing the scan.

## See also

- `scan.md`, `scan2d.md` — lower-dim analogs.
- `all-scan-list.md` — staged sequential scans inside the full
  pipeline (avoids the 3D grid blow-up when stages are decoupled).
- Defaults: `import mlmm.core.defaults as d; print(d.OUT_DIR_SCAN3D)`