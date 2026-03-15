# `scan2d`

## Overview

> **Summary:** Perform a two-distance (d1, d2) grid scan with harmonic restraints and ML/MM relaxations. Use `-s/--scan-lists` with a YAML/JSON spec file (recommended) or an inline Python literal.

`mlmm scan2d` constructs linear grids for two bond distances using `--max-step-size`, relaxes each grid point with the appropriate restraints active, and records unbiased ML/MM energies for visualization. The scan iterates d1 first, relaxing the structure with only the d1 restraint active, then iterates d2 for each d1 value with both restraints applied.

Energies at each grid point are re-evaluated without the bias to populate a PES grid and contour plot. Outputs include per-point XYZ snapshots, `surface.csv` summarizing the PES, a 2D contour map (`scan2d_map.png`), and a 3D landscape with bottom projection (`scan2d_landscape.html`).


## Minimal example
```bash
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan2d.yaml --print-parsed -o ./result_scan2d/
```

## Output checklist
- `result_scan2d/surface.csv`
- `result_scan2d/grid/point_i000_j000.xyz`
- `result_scan2d/scan2d_map.png` and `result_scan2d/scan2d_landscape.html`

## Common examples
1. Validate parsed scan targets from a YAML spec.
2. Run with an inline `-s` literal.
3. Enable `--dump` to store inner trajectories per outer d1 step.

> **Note:** Add `--print-parsed` when you want to verify parsed pair targets from `-s/--scan-lists`.

## Usage
```bash
mlmm scan2d -i INPUT.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q CHARGE [-m SPIN] \
 [-s scan2d.yaml | -s "[(I1,J1,LOW1,HIGH1),(I2,J2,LOW2,HIGH2)]"] \
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
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan2d.yaml --print-parsed

# Alternative: inline Python literal
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,1.30,3.10),(10,55,1.20,3.20)]"

# LBFGS scan with TRJ dumps and fixed color scale for the contour plot
mlmm scan2d -i input.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s "[(12,45,1.30,3.10),(10,55,1.20,3.20)]" \
 --max-step-size 0.20 --dump -o ./result_scan2d/ --preopt --baseline min \
 --zmin 0.0 --zmax 40.0
```

## YAML/JSON spec format (recommended)

`-s/--scan-lists` auto-detects YAML/JSON files. Pass a file path to use the spec format:

```yaml
one_based: true # optional; defaults to CLI --one-based/--zero-based
pairs:
 - [12, 45, 1.30, 3.10]
 - [10, 55, 1.20, 3.20]
```

- `pairs` is required and must contain exactly 2 quadruples.
- Each quadruple is `(i, j, low_A, high_A)`.
- Indices may be integers or PDB selectors (same as inline literals).

## Inline literal format

When `-s/--scan-lists` receives a value that is not a file path, it is treated as a **single Python literal** string. Shell quoting matters.

### Basic structure

The literal is a Python list of exactly **two** quadruples `(atom1, atom2, low_A, high_A)`:

```
-s '[(atom1, atom2, low_A, high_A), (atom3, atom4, low_A, high_A)]'
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
-s '[("TYR,285,CA","MMT,309,C10",1.30,3.10),("TYR,285,CB","MMT,309,C11",1.20,3.20)]'

# Correct: integer indices need no inner quotes
-s '[(1, 5, 1.30, 3.10), (2, 8, 1.20, 3.20)]'

# Avoid: double-quoting the outer literal requires escaping inner quotes
-s "[(\"TYR,285,CA\",\"MMT,309,C10\",1.30,3.10),...]"
```

## Workflow
1. **Input & preoptimization** -- Load the enzyme PDB, resolve charge/spin, build the ML/MM calculator (MLIP backend + hessian_ff; backend selected via `-b/--backend`, default `uma`), and optionally run an unbiased pre-optimization when `--preopt`. When `--embedcharge` is enabled, xTB point-charge embedding is applied for MM-to-ML environmental corrections.
2. **Grid construction** -- Parse targets from `-s/--scan-lists` (YAML/JSON spec file or inline literal) into two quadruples, normalize indices (1-based by default or PDB atom selectors like `"TYR,285,CA"`). Build linear grids with `ceil(|high - low| / h) + 1` points where `h = --max-step-size`.
3. **Outer loop (d1)** -- For each d1 value, relax the system with **only the d1 restraint** active.
4. **Inner loop (d2)** -- For each d2 value at the current d1, relax with **both restraints** active starting from the nearest previously converged structure.
5. **Energy evaluation** -- At each (i, j) pair, evaluate the ML/MM energy without bias and record to `surface.csv`.
6. **Visualization** -- Write `scan2d_map.png` (2D contour) and `scan2d_landscape.html` (3D surface). Use `--zmin/--zmax` to clamp the color scale. Baselines: `--baseline min` zeroes the minimum energy; `--baseline first` zeroes the (i=0, j=0) grid point.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input enzyme complex PDB (required). | Required |
| `--parm PATH` | Amber parm7 topology for the enzyme (required). | Required |
| `--model-pdb PATH` | PDB defining the ML region. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | ML-region total charge. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g., `GPP:-3,SAM:1`). Derives total charge when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze. | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. Can be combined with `--detect-layer`. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. Providing this disables `--detect-layer`. | _None_ |
| `-s, --scan-lists TEXT` | Scan targets: a YAML/JSON spec file path (auto-detected, with `pairs` containing 2 quadruples) or an inline Python literal `"[(i1,j1,low1,high1),(i2,j2,low2,high2)]"`. Indices can be integers or PDB atom selectors. | Required |
| `--one-based / --zero-based` | Interpret `(i,j)` indices in `-s/--scan-lists` as 1-based or 0-based. | `True` (1-based) |
| `--print-parsed/--no-print-parsed` | Print parsed pair tuples after `-s/--scan-lists` resolution. | `False` |
| `--max-step-size FLOAT` | Maximum distance increment per step (Å). Determines grid density. | `0.20` |
| `--bias-k FLOAT` | Harmonic well strength k (eV/Å²). | `300.0` |
| `--relax-max-cycles INT` | Maximum LBFGS cycles per biased relaxation. | `10000` |
| `--dump/--no-dump` | Write inner d2 scan TRJs per d1 slice. | `False` |
| `-o, --out-dir TEXT` | Base output directory. | `./result_scan2d/` |
| `--thresh TEXT` | Convergence preset (`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`). | `baker` |
| `--config FILE` | Base YAML configuration file (applied first). | _None_ |
| `--ref-pdb FILE` | Reference PDB topology when `--input` is XYZ. | _None_ |
| `--preopt/--no-preopt` | Run an unbiased pre-optimization before scanning. | `True` |
| `--baseline {min,first}` | Reference for relative energy (kcal/mol). | `min` |
| `--zmin FLOAT` | Lower bound of the contour color scale (kcal/mol). | Autoscaled |
| `--zmax FLOAT` | Upper bound of the contour color scale (kcal/mol). | Autoscaled |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | _None_ (`uma` applied internally) |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## Outputs
```
out_dir/ (default: ./result_scan2d/)
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
 out_dir: ./result_scan2d/
lbfgs:
 max_step: 0.3
 out_dir: ./result_scan2d/
bias:
 k: 300.0
```

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [scan](scan.md) -- 1D bond-length driven scan
- [scan3d](scan3d.md) -- 3D distance grid scan
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow
