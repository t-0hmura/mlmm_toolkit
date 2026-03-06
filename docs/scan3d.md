# `scan3d`

## Overview

> **Summary:** Perform a three-distance (d1, d2, d3) grid scan with harmonic restraints and ML/MM relaxations. Use `--spec` (YAML/JSON, recommended) or `--scan-lists`.

### At a glance
- **Input:** One full enzyme PDB + `--spec scan3d.yaml` (recommended), or one `--scan-lists` literal (three quadruples). Use `--csv` for plot-only mode.
- **Grid ordering:** d1 is scanned first, then d2 for each d1 value (both restraints active), then d3 for each (d1, d2) with all three restraints active.
- **Energies:** Recorded energies are evaluated **without bias**, so grid points are directly comparable.
- **Outputs:** `surface.csv`, per-point geometries under `grid/`, and an HTML isosurface plot (`scan3d_density.html`).
- **Caution:** 3D grids grow very quickly; consider coarser `--max-step-size` or smaller ranges first.

`mlmm scan3d` nests loops over d1, d2, and d3, relaxing each point with the appropriate restraints active using the ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`). The ML region comes from `--model-pdb`; Amber parameters are read from `--parm`; the MLIP backend is selected via `--backend` (default: `uma`); the optimizer is PySisyphus LBFGS.


## Minimal example
```bash
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan3d.yaml --print-parsed --out-dir ./result_scan3d/
```

## Output checklist
- `result_scan3d/surface.csv`
- `result_scan3d/grid/point_i000_j000_k000.xyz`
- `result_scan3d/scan3d_density.html`

## Common examples
1. Validate parsed `pairs` from a YAML spec before running a full grid.
2. Run with `--scan-lists`.
3. Enable `--dump` to keep inner d3 trajectories for each `(d1,d2)` slice.

> **Note:** Add `--print-parsed` when you want to verify parsed pair targets from `--spec` / `--scan-lists`.

## Usage
```bash
mlmm scan3d -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m MULT] \
 [--csv precomputed_surface.csv] \
 [--spec scan3d.yaml | --scan-lists "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2),(I3,J3,LOW3,HIGH3)]"] \
 [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
 [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
 [--dump/--no-dump] [--out-dir DIR] \
 [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

### Examples
```bash
# Recommended: YAML/JSON spec
cat > scan3d.yaml << 'YAML'
one_based: true
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
 - [15, 60, 1.10, 3.00]
YAML
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan3d.yaml --print-parsed

# Alternative: Python literal
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]"

# With pre-optimization and custom output directory
mlmm scan3d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20),(15,60,1.10,3.00)]" \
 --max-step-size 0.20 --dump --out-dir ./result_scan3d/ \
 --preopt --baseline min
```

## `--spec` format (recommended)

```yaml
one_based: true # optional; defaults to CLI --one-based/--zero-based
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
 - [15, 60, 1.10, 3.00]
```

- `pairs` is required and must contain exactly 3 quadruples.
- Each quadruple is `(i, j, low_A, high_A)`.
- Indices may be integers or PDB selectors (same as `--scan-lists`).

## `--scan-lists` format

`--scan-lists` is the advanced input mode. It accepts a **single Python literal** string. Shell quoting matters.

### Basic structure

The literal is a Python list of exactly **three** quadruples `(atom1, atom2, low_A, high_A)`:

```
--scan-lists '[(atom1, atom2, low_A, high_A), (atom3, atom4, low_A, high_A), (atom5, atom6, low_A, high_A)]'
```

- Wrap the entire literal in **single quotes** so the shell does not interpret parentheses or spaces.
- Each quadruple defines one scan axis: the distance between `atom1`--`atom2` is scanned from `low_A` to `high_A`.
- Unlike `scan`, only **one literal** is accepted (no multi-stage support).

### Specifying atoms

Atoms can be given as **integer indices** or **PDB selector strings**:

| Method | Example | Notes |
| --- | --- | --- |
| Integer index | `(1, 5, 1.30, 3.10)` | 1-based by default (`--one-based`) |
| PDB selector | `("TYR,285,CA", "MMT,309,C10", 1.30, 3.10)` | Residue name, residue number, atom name |

PDB selector tokens can be separated by any of: comma `,`, space, slash `/`, backtick `` ` ``, or backslash `\`. Token order is flexible.

```bash
# All of these specify the same atom:
"TYR,285,CA"
"TYR 285 CA"
"TYR/285/CA"
"285,TYR,CA" # order is flexible
```

### Quoting rules

```bash
# Correct: single-quote the list, double-quote selector strings inside
--scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20),("TYR,285,CG","MMT,309,C12",1.10,3.00)]'

# Correct: integer indices need no inner quotes
--scan-lists '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20), (3, 12, 1.10, 3.00)]'

# Avoid: double-quoting the outer literal requires escaping inner quotes
--scan-lists "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10),...]"
```

## Workflow
1. Load the structure through `geom_loader`, resolve charge/spin from CLI, and
 optionally run an unbiased preoptimization when `--preopt`.
2. Parse targets from `--spec` (recommended) or `--scan-lists` (default 1-based indices unless
 `--zero-based` is passed) into three quadruples. For PDB inputs, each atom
 entry can be an integer index or a selector string like `"TYR,285,CA"`;
 delimiters may be spaces, commas, slashes, backticks, or backslashes.
3. Outer loop over `d1[i]`: relax with only the d1 restraint active, starting
 from the previously scanned geometry whose d1 value is closest.
4. Middle loop over `d2[j]`: relax with d1 and d2 restraints, starting from the
 closest (d1, d2) geometry.
5. Inner loop over `d3[k]`: relax with all three restraints, measure the
 unbiased energy (bias removed for evaluation), and write the constrained
 geometry and convergence flag.
6. After the scan completes, assemble `surface.csv`, apply the kcal/mol
 baseline shift (`--baseline {min|first}`), and generate a 3D RBF-interpolated
 isosurface plot (`scan3d_density.html`) honoring `--zmin/--zmax`.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme PDB (no link atoms). | Required unless `--csv` |
| `--parm PATH` | Amber parm7 topology for the full enzyme. | Required unless `--csv` |
| `--model-pdb PATH` | PDB defining the ML region. | _None_ |
| `--model-indices TEXT` | Explicit ML-region atom indices (alternative to `--model-pdb`). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Indexing convention for `--model-indices`. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors. | `True` |
| `-q, --charge INT` | ML-region total charge. | Required unless `--csv` |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | 1-based comma-separated frozen atom indices. | _None_ |
| `--hess-cutoff FLOAT` | Cutoff distance for Hessian-MM layer. | _None_ |
| `--movable-cutoff FLOAT` | Cutoff distance for movable-MM layer. | _None_ |
| `--spec FILE` | YAML/JSON spec with `pairs` (3 quadruples) and optional `one_based`. | Recommended |
| `--csv FILE` | Load precomputed `surface.csv` and generate plot without running a scan. | _None_ |
| `--scan-lists TEXT` | single Python literal with three quadruples `(i,j,low,high)`. `i`/`j` can be integer indices or PDB atom selectors. | Alternative to `--spec` |
| `--one-based / --zero-based` | Interpret `(i, j)` indices as 1- or 0-based. | `True` (1-based) |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `--spec`/`--scan-lists` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum distance increment per step (Å). Controls grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic well strength k (eV/Å²). | `300.0` |
| `--relax-max-cycles INT` | Maximum optimizer cycles during each biased relaxation. | `10000` |
| `--dump/--no-dump` | Write inner d3 scan TRJs per (d1, d2) slice. | `False` |
| `--out-dir TEXT` | Output directory root for grids and plots. | `./result_scan3d/` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `--ref-pdb FILE` | Reference PDB topology for non-PDB inputs. | _None_ |
| `--preopt/--no-preopt` | Run an unbiased optimization before scanning. | `True` |
| `--baseline {min,first}` | Shift kcal/mol energies so the global min or `(i,j,k)=(0,0,0)` is zero. | `min` |
| `--zmin FLOAT` | Manual lower limit for the isosurface color bands (kcal/mol). | Autoscaled |
| `--zmax FLOAT` | Manual upper limit for the isosurface color bands (kcal/mol). | Autoscaled |
| `--backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## Outputs
```
out_dir/ (default:./result_scan3d/)
 surface.csv # Grid metadata (d1, d2, d3, energy, convergence)
 scan3d_density.html # 3D energy isosurface visualization
 grid/point_i###_j###_k###.xyz # Relaxed geometry for each grid point
 grid/point_i###_j###_k###.pdb # PDB companions (B-factors: ML=0, Movable-MM=10, Frozen=20)
 grid/inner_path_d1_###_d2_###_trj.xyz # Present only when --dump is True
```

## YAML configuration

```yaml
geom:
 coord_type: cart
 freeze_atoms: []
calc:
 charge: 0
 spin: 1
mlmm:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
opt:
 thresh: baker
 max_cycles: 10000
 dump: false
 out_dir: ./result_scan3d/
lbfgs:
 max_step: 0.3
 out_dir: ./result_scan3d/
bias:
 k: 300.0
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`) reuses the same
 `HarmonicBiasCalculator` as the 1D/2D scans.
- `--baseline` defaults to the global minimum; `--baseline first` anchors the
 `(i,j,k)=(0,0,0)` grid point when present.
- 3D visualization uses RBF interpolation on a 50x50x50 grid with
 semi-transparent step-colored isosurfaces.
- When the input is a PDB, each grid-point XYZ and (if present) inner-path TRJ are also
 converted to PDB files, using the input PDB as a template. B-factors are annotated
 as: ML=0, MovableMM=10, FrozenMM=20.
- Plot color scales can be clamped with `--zmin/--zmax` to compare scans consistently.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [scan](scan.md) -- 1D bond-length driven scan
- [scan2d](scan2d.md) -- 2D distance grid scan
- [opt](opt.md) -- Geometry optimization (often precedes scan)
- [all](all.md) -- End-to-end workflow
- [YAML Reference](yaml_reference.md) -- Full scan configuration options
