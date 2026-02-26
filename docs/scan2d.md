# `scan2d`

## Overview

> **Summary:** Perform a two-distance (d1, d2) grid scan with harmonic restraints and ML/MM relaxations. Use `--spec` (YAML/JSON, recommended) or `--scan-lists`.

### At a glance
- **Input:** One structure + `--spec scan2d.yaml` (recommended), or one `--scan-lists` literal containing exactly two quadruples.
- **Grid ordering:** Each axis is reordered so the point closest to the (pre)optimized structure is visited first.
- **Energies:** Values written to `surface.csv` are always evaluated **without bias**, so grid points are directly comparable.
- **Outputs:** `surface.csv` plus `scan2d_map.png` and `scan2d_landscape.html`, and per-point structures under `grid/`.
- **Caution:** Grid size grows quickly as `(high - low) / --max-step-size` increases.

`mlmm scan2d` constructs linear grids for two bond distances using `--max-step-size`, relaxes each grid point with the appropriate restraints active, and records unbiased ML/MM energies for visualization. The scan iterates d1 first, relaxing the structure with only the d1 restraint active, then iterates d2 for each d1 value with both restraints applied.

Energies at each grid point are re-evaluated without the bias to populate a PES grid and contour plot. Outputs include per-point XYZ snapshots, `surface.csv` summarizing the PES, a 2D contour map (`scan2d_map.png`), and a 3D landscape with bottom projection (`scan2d_landscape.html`).


## Minimal example
```bash
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan2d.yaml --print-parsed --out-dir ./result_scan2d/
```

## Output checklist
- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` and `result_scan2d/scan2d_landscape.html`

## Common examples
1. Validate parsed scan targets from a YAML spec.
2. Run with the `--scan-lists` literal.
3. Enable `--dump` to store inner trajectories per outer d1 step.

> **Note:** Add `--print-parsed` when you want to verify parsed pair targets from `--spec` / `--scan-lists`.

## Usage
```bash
mlmm scan2d -i INPUT.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m SPIN] \
 [--spec scan2d.yaml | --scan-lists "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]"] \
 [--one-based|--zero-based] [--max-step-size FLOAT] [--bias-k FLOAT] \
 [--freeze-atoms "1,3,5"] [--relax-max-cycles INT] [--thresh PRESET] \
 [--dump/--no-dump] [--out-dir DIR] \
 [--preopt/--no-preopt] [--baseline {min|first}] [--zmin FLOAT] [--zmax FLOAT]
```

### Examples
```bash
# Recommended: YAML/JSON spec
cat > scan2d.yaml << 'YAML'
one_based: true
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
YAML
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --spec scan2d.yaml --print-parsed

# Alternative: Python literal
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

# LBFGS scan with TRJ dumps and fixed color scale for the contour plot
mlmm scan2d -i input.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --scan-lists "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
 --max-step-size 0.20 --dump --out-dir ./result_scan2d/ --preopt --baseline min \
 --zmin 0.0 --zmax 40.0
```

## `--spec` format (recommended)

```yaml
one_based: true # optional; defaults to CLI --one-based/--zero-based
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
```

- `pairs` is required and must contain exactly 2 quadruples.
- Each quadruple is `(i, j, low_A, high_A)`.
- Indices may be integers or PDB selectors (same as `--scan-lists`).

## `--scan-lists` format

`--scan-lists` is the advanced input mode. It accepts a **single Python literal** string. Shell quoting matters.

### Basic structure

The literal is a Python list of exactly **two** quadruples `(atom1, atom2, low_A, high_A)`:

```
--scan-lists '[(atom1, atom2, low_A, high_A), (atom3, atom4, low_A, high_A)]'
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
--scan-lists '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20)]'

# Correct: integer indices need no inner quotes
--scan-lists '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20)]'

# Avoid: double-quoting the outer literal requires escaping inner quotes
--scan-lists "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10),...]"
```

## Workflow
1. **Input & preoptimization** -- Load the enzyme PDB, resolve charge/spin, build the ML/MM calculator (FAIR-Chem UMA + hessian_ff), and optionally run an unbiased pre-optimization when `--preopt`.
2. **Grid construction** -- Parse targets from `--spec` (recommended) or `--scan-lists` into two quadruples, normalize indices (1-based by default or PDB atom selectors like `"TYR,285,CA"`). Build linear grids with `ceil(|high - low| / h) + 1` points where `h = --max-step-size`.
3. **Outer loop (d1)** -- For each d1 value, relax the system with **only the d1 restraint** active.
4. **Inner loop (d2)** -- For each d2 value at the current d1, relax with **both restraints** active starting from the nearest previously converged structure.
5. **Energy evaluation** -- At each (i, j) pair, evaluate the ML/MM energy without bias and record to `surface.csv`.
6. **Visualization** -- Write `scan2d_map.png` (2D contour) and `scan2d_landscape.html` (3D surface). Use `--zmin/--zmax` to clamp the color scale. Baselines: `--baseline min` zeroes the minimum energy; `--baseline first` zeroes the (i=0, j=0) grid point.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input enzyme complex PDB (required). | Required |
| `--real-parm7 PATH` | Amber parm7 topology for the enzyme (required). | Required |
| `--model-pdb PATH` | PDB defining the ML region. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | ML-region total charge. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze. | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (A) for MM Hessian atoms. Providing cutoffs disables `--detect-layer`. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (A) for movable MM atoms. | _None_ |
| `--spec FILE` | YAML/JSON spec with `pairs` (2 quadruples) and optional `one_based`. | Recommended |
| `--scan-lists TEXT` | Python literal with two quadruples: `"[(i1,j1,low1,high1),(i2,j2,low2,high2)]"`. Indices can be integers or PDB atom selectors. | Alternative to `--spec` |
| `--one-based / --zero-based` | Interpret `(i,j)` indices in `--scan-lists` as 1-based or 0-based. | `True` (1-based) |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `--spec`/`--scan-lists` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum distance increment per step (A). Determines grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic well strength k (eV/A^2). | `100.0` |
| `--relax-max-cycles INT` | Maximum LBFGS cycles per biased relaxation. | `10000` |
| `--dump/--no-dump` | Write inner d2 scan TRJs per d1 slice. | `False` |
| `--out-dir TEXT` | Base output directory. | `./result_scan2d/` |
| `--thresh TEXT` | Convergence preset (`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`). | _None_ |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `--ref-pdb FILE` | Reference PDB topology when `--input` is XYZ. | _None_ |
| `--preopt/--no-preopt` | Run an unbiased pre-optimization before scanning. | `True` |
| `--baseline {min,first}` | Reference for relative energy (kcal/mol). | `min` |
| `--zmin FLOAT` | Lower bound of the contour color scale (kcal/mol). | Autoscaled |
| `--zmax FLOAT` | Upper bound of the contour color scale (kcal/mol). | Autoscaled |

## Outputs
```
out_dir/ (default:./result_scan2d/)
├── surface.csv # PES grid: i, j, d1_A, d2_A, energy_hartree, energy_kcal, bias_converged
├── scan2d_map.png # 2D contour map
├── scan2d_landscape.html # 3D surface visualization (Plotly)
├── grid/
│ ├── point_i###_j###.xyz # Relaxed geometry for every (i, j) pair
│ ├── point_i###_j###.pdb # PDB companion (when input is PDB)
│ ├── preopt_i###_j###.xyz # Pre-optimized structure (when --preopt)
│ └── inner_path_d1_###_trj.xyz # Inner d2 trajectory per d1 slice (when --dump)
└── (stdout) # Progress and energy summaries
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
 out_dir:./result_scan2d/
lbfgs:
 max_step: 0.3
 out_dir:./result_scan2d/
bias:
 k: 100.0
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`) keeps the entire enzyme complex. The ML region comes from `--model-pdb`; Amber parameters are read from `--real-parm7`.
- The bias is always removed before final energies are recorded, so `surface.csv` is directly comparable across grid points.
- When the input is a PDB, each grid-point XYZ and (if present) inner-path TRJ are also converted to PDB files with B-factor annotations: ML-region atoms = 100.00, frozen atoms = 50.00, both = 150.00.
- `i`/`j` entries in `--scan-lists` may be integer indices (1-based by default) or PDB atom selectors like `"TYR,285,CA"`.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [scan](scan.md) -- 1D bond-length driven scan
- [scan3d](scan3d.md) -- 3D distance grid scan
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow
