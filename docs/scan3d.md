# `scan3d` (3D restrained grid scan)

`scan3d` relaxes every point of a grid over three coordinates of a layered enzyme structure with harmonic restraints and the ML/MM calculator, records the energy without the restraints, and draws the energy volume as isosurfaces in an HTML page. Each axis is a range: a distance `(i,j,low,high)` in Å, an angle `(i,j,k,low,high)`, or a dihedral `(i,j,k,l,low,high)` in degrees.

---

## What it is for

* **Reactions that involve three coordinates at once**: see the energy landscape when, for example, a bond forms, another breaks, and a proton moves in the same step.
* **Redrawing a finished grid**: plot an existing `surface.csv` again over another energy range (`--csv`).
* **Judging the result**: `surface.csv` marks every point as converged or not, and `result.json` gives `scientific_status` (see *Reading surface.csv* below).

The ML region uses **UMA** (Meta) by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or DFT (`dft`). For a single path driven by one or more coordinates, use [`scan`](scan.md); for two coordinates, use [`scan2d`](scan2d.md).

---

## Examples

Here `pocket.pdb` contains the full system matching `real.parm7`, and `ml_region.pdb` selects its ML region without link hydrogens.

The atoms are those of the bundled enzyme example (`examples/beza/1.R.pdb`), whose PDB has an empty chain column. Each atom is therefore written with three fields, residue name, residue number, and atom name, in any order and separated by commas or spaces (`"SAM,320,CS1"`).

### 1. From a YAML spec

Write the three ranges under `pairs:` and add `--out-json` to also get `result.json`. With the default step of 0.2 Å, this file gives a 9 × 9 × 7 grid (567 points).

```yaml
# scan3d.yaml
pairs:
  - ["SAM,320,CS1", "GPP,321,C7", 1.50, 3.00]
  - ["GPP,321,H11", "GLU,186,OE2", 0.90, 2.50]
  - ["SAM,320,SD", "SAM,320,CS1", 1.80, 3.00]
```

```bash
mlmm scan3d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan3d.yaml --out-json -o ./result_scan3d/
```

Open `result_scan3d/scan3d_density.html` in a browser for the isosurfaces; `result.json` gives `scientific_status` and the number of usable points (`n_points_usable`).

### 2. Inline literal

The same three ranges can be written on the command line as one literal.

```bash
mlmm scan3d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50),("SAM,320,SD","SAM,320,CS1",1.80,3.00)]'
```

### 3. L-BFGS, dump, and pre-optimization

Optimize the input before the scan, relax each point with L-BFGS, keep the inner-loop trajectories, and measure the relative energies from the lowest usable point. `--max-step-size 0.20`, `--opt-mode grad`, and `--baseline min` are the defaults, written out here.

```bash
mlmm scan3d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50),("SAM,320,SD","SAM,320,CS1",1.80,3.00)]' \
    --max-step-size 0.20 --dump -o ./result_scan3d/ --opt-mode grad \
    --preopt --baseline min
```

### 4. Re-plot an existing surface.csv

Redraw the isosurfaces of a finished grid over −10 to 40 kcal/mol; no energy is computed. A separate `-o` keeps the files of the original scan (see *Notes*).

```bash
mlmm scan3d --csv ./result_scan3d/surface.csv --zmin -10 --zmax 40 -o ./result_scan3d_replot/
```

---

## How it works

1. **Starting structure and grid**:
The ML-region charge comes from `-q` or `-l` (see {ref}`Charge specification <charge-specification>`). With `--preopt`, the input is first optimized without restraints; if that does not converge, the input geometry is used. Each axis gets ceil(|high − low| / h) + 1 evenly spaced values, both ends included, where h is `--max-step-size` (Å) for a distance and `--max-angle-step-size` or `--max-dihedral-step-size` (degrees) for an angle or dihedral. The values are visited from the one closest to the starting structure.
2. **Three nested loops**:
For each d₁ value, the structure is relaxed with only the d₁ restraint; for each d₂ value, with the d₁ and d₂ restraints; the inner loop then scans d₃ with all three restraints. Each relaxation starts from the nearest structure that has already converged in the same loop; until one has, it starts from the structure the enclosing loop produced (the starting structure for d₁).
3. **Relaxation at each point**:
A harmonic restraint E = ½ k (q − q_target)² holds each coordinate q at its target (k is `--restraint-k`), and L-BFGS (`--opt-mode grad`, default) or RFO (rational function optimization, `--opt-mode hess`) relaxes the rest with the ML/MM calculator; the atoms of the frozen MM layer stay fixed. The energy is then computed without the restraints, and the structure is written under `grid/`.
4. **Table and figure**:
After the last point, all points go into `surface.csv`. The usable points are interpolated with radial basis functions (RBF) on a 50 × 50 × 50 grid, and eight semi-transparent isosurfaces with banded colors are drawn in `scan3d_density.html`. With `--csv`, only this step runs, on the given table.

---

## Reading surface.csv

`surface.csv` has one row per grid point and one reference row.

| Column | Meaning |
| --- | --- |
| `i`, `j`, `k` | Grid indices. Index 0 is the value closest to the starting structure, so the indices follow the visiting order, not ascending values |
| `d1_A`, `d2_A`, `d3_A` (also `q1`, `q2`, `q3`) | Coordinate values measured after the relaxation. The `_A` names are kept for every axis, so an angle axis holds degrees; `q1_unit`, `q2_unit`, and `q3_unit` give the unit (`angstrom` or `degree`) |
| `target_d1_A`, `target_d2_A`, `target_d3_A` (also `target_q1`, `target_q2`, `target_q3`) | Restraint targets of the point |
| `energy_hartree` | Energy without the restraints (Hartree) |
| `bias_converged` | Whether the restrained relaxation converged |
| `is_preopt` | `true` only for the reference row |
| `energy_kcal` | Energy relative to the baseline (kcal/mol) |
| `d1_label`, `d2_label`, `d3_label` | Axis labels used in the figure |

* **Reference row**: `i = j = k = -1` and `is_preopt = true` hold the starting structure. The row stays in the table but is never a grid point, a baseline, or a plotted point.
* **Usable points**: a point is usable when its relaxation converged, its energy is finite, and its structure was written. Only usable points set the baseline and enter the figure.
* **Verdict**: in `result.json` (`--out-json`), `scientific_status` is `success` when every grid point is usable, `partial` when only some are (exit status 0), and `failed` when none is (exit status 1); `n_points_attempted` and `n_points_usable` give the counts (see {ref}`Exit codes <exit-codes>`).
* **Next step**: the isosurfaces are interpolated, so take a structure from a computed point near the saddle, `grid/point_*.pdb`, and pass it to [`tsopt`](tsopt.md); points in the reactant and product basins can be inputs for [`path-search`](path-search.md).
* **Re-plotting a table (`--csv`)**: the table needs `d1_A`, `d2_A`, `d3_A`, and `energy_hartree` or `energy_kcal`. The reference row and rows with `bias_converged = false` or a non-finite energy are left out.

---

## Output files

After a run, `--out-dir` (default: `./result_scan3d/`) contains the following files.

```text
result_scan3d/
├─ surface.csv                          # grid table with the reference row
├─ scan3d_density.html                  # 3D isosurfaces (open in a browser)
├─ grid/
│  ├─ point_i150_j090_k180.xyz          # relaxed structure of each grid point
│  ├─ preopt_iDDD_jDDD_kDDD.xyz         # starting structure (the reference row)
│  └─ inner_path_d1_000_d2_000_trj.xyz  # inner-loop trajectory for each (d₁, d₂) pair (--dump)
└─ result.json                          # summary (--out-json); summary.json holds the same content
```

Start with `scan3d_density.html` and `surface.csv`; the structure of each point is under `grid/`. In `result.json`, `grid_points[]` maps each grid index to its values, targets, energy, convergence, and structure file.

* **File names**: the number after `i`, `j`, or `k` (the tag `DDD`) is the target × 100 (Å, or degrees for an angle), padded to at least three digits, not the grid index of `surface.csv`: `d1 = 1.50 Å, d2 = 0.90 Å, d3 = 1.80 Å` gives `point_i150_j090_k180.xyz`, and an angle of 120° gives `12000`. When two points round to the same tag, the later file name gets `_grid_III_JJJ_KKK` with the zero-based indices; the numbers in `inner_path_d1_000_d2_000` are zero-based indices too.
* **Other formats**: each structure is also written as `.pdb` (`--no-convert-files` turns this off), with the layer of each atom in the B-factor column: 0 for the ML region, 10 for movable MM, and 20 for frozen MM. mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files with the original identifiers (see {ref}`mmCIF input <mmcif-input>`).
* **With `--csv`**: only `scan3d_density.html` is written, plus `result.json` (without `grid_points`) with `--out-json`. `scientific_status` is `success` when the table has a usable point and `failed` otherwise.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `scan3d`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | `None` | Full-system structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`). Required unless `--csv` is given |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` or `--csv` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Total charge of the unknown ligand residues (e.g. `-1`) or per-residue charges (e.g. `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-s, --scan-lists` | text | `None` | Three ranges as a YAML/JSON spec file or one inline literal: distance `(i,j,low,high)`, angle `(i,j,k,low,high)`, or dihedral `(i,j,k,l,low,high)`. Required unless `--csv` is given |
| `-o, --out-dir` | path | `./result_scan3d/` | Output directory |
| `--max-step-size` | float | `0.2` | Largest grid spacing of a distance axis (Å) |
| `--max-angle-step-size` | float | `5.0` | Largest grid spacing of an angle axis (degrees) |
| `--max-dihedral-step-size` | float | `10.0` | Largest grid spacing of a dihedral axis (degrees) |
| `--restraint-k` | float | `300.0` | Restraint strength k (eV/Å² for distances, eV/rad² for angles); alias `--bias-k`. When omitted, YAML `bias.k` applies |
| `--opt-mode` | `grad` / `hess` | `grad` | Relaxation of each point: L-BFGS / RFO |
| `--preopt/--no-preopt` | flag | `False` | Optimize the input without restraints before the scan |
| `--dump/--no-dump` | flag | `False` | Write the inner-loop (d₃) trajectory for each (d₁, d₂) pair under `grid/` |
| `--baseline` | `min` / `first` | `min` | Zero of `energy_kcal`: lowest usable point, or point `(0, 0, 0)` |
| `--zmin`, `--zmax` | float | interpolated min / max | Lower and upper ends of the energy range over which the eight isosurfaces are placed (kcal/mol) |
| `--thresh` | text | `baker` | Convergence preset of each relaxation (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`) |
| `--csv` | path | `None` | Read a finished `surface.csv` and only draw the figure; `-i`, `--parm7`, `-s`, and `-q` are not needed |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/scan3d.md) for every option.

> **Note:** In YAML (`--config`), the `lbfgs` (or `opt.lbfgs`) section applies with `--opt-mode grad` and the `rfo` (or `opt.rfo`) section with `--opt-mode hess`. Trajectories are written only with `--dump`; YAML `opt.dump` is not used.

---

## Notes

* **Three ranges in one literal**: `-s` takes exactly three ranges, in one inline literal or under `pairs:` in a YAML/JSON file. For staged scans, use [`scan`](scan.md).
* **For a PDB with chains**, the positional form `A:SAM:320:CS1` picks one atom without ambiguity (see {ref}`Scan-list spec <scan-list-spec>`).
* **Grid size**: the number of relaxations is the product of the three axis lengths and grows quickly (567 for example 1). Start with a larger `--max-step-size` or narrower ranges.
* **`--baseline first`**: zero is put at point `(i, j, k) = (0, 0, 0)` when it is usable; otherwise the run prints `[baseline] 'first' requested but usable (i=0,j=0,k=0) is missing; using the usable minimum instead.` and uses the lowest usable point.
* **Frozen atoms**: the atoms given by `--freeze-atoms` or YAML `geom.freeze_atoms`, and the atoms of the frozen MM layer, stay fixed in every relaxation. A scanned coordinate whose atoms are all frozen is an error (see {ref}`Frozen atoms and restraints <freeze-atoms-and-restraints>`).
* **Check the spec without computing**: `--dry-run` reads the input, the charge and spin, and `-s`, prints the plan, and exits without any optimization. With `--csv`, it checks only the options.
* **Cycle limit**: `--relax-max-cycles` (default `100000`) limits each relaxation; an explicit value overrides YAML `opt.max_cycles`, and without it the YAML value applies.
* **Too few usable points**: with fewer than four usable points, or with all of them in one plane, only the figure is skipped. The run prints `[plot] NOTE: Volume plot skipped: …`, still writes `surface.csv` (and `result.json` with `--out-json`, listing no figure), and exits with status 0. With no usable point it prints `[plot] No finite data for plotting.` and exits with status 1.
* **Re-plot warning**: a table without the `bias_converged` or `artifact_written` column is still read, with `[plot] WARNING: CSV lacks complete point provenance; …`; the `surface.csv` written by `scan3d` has no `artifact_written` column, so re-plotting it prints this warning, and `result.json` then has no `n_points_usable`.
* **Re-plotting into the scan directory**: with or without `--csv`, a run without `--out-json` removes `result.json` and `summary.json` from `--out-dir`, a run with it overwrites them, and every run removes the old `scan3d_density.html` before drawing a new one. To keep the files of the original scan, give `--csv` another `-o`, as in example 4.

---

## See also

* [scan](scan.md) — staged scans of one or more coordinates from one structure
* [scan2d](scan2d.md) — energy map over two coordinates
* [opt](opt.md) — single-structure optimization before or after a scan
* [tsopt](tsopt.md) — optimize a structure near the saddle into a TS
* [path-search](path-search.md) — MEP through structures taken from the grid
* [all](all.md) — end-to-end workflow
* [Troubleshooting](troubleshooting.md) — diagnosing failed runs
