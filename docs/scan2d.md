# `scan2d` (2D restrained grid scan)

## Overview

`scan2d` relaxes every point of a grid over two coordinates of a layered enzyme structure with harmonic restraints and the ML/MM calculator, and records the energy without the restraints, giving a 2D energy map of the reaction. Each axis is a range: a distance `(i,j,low,high)` in Å, an angle `(i,j,k,low,high)`, or a dihedral `(i,j,k,l,low,high)` in degrees.

### What it is for

* **Locating the TS region**: see where the saddle between the reactant and product basins lies before a path search or TS optimization.
* **Viewing the landscape before the MEP**: check whether two events, such as bond formation and proton transfer, happen together or one after the other, before refining the minimum-energy path (MEP).
* **Judging the result**: `surface.csv` marks every point as converged or not, and `result.json` gives `scientific_status` (see *Reading surface.csv* below).

The ML region uses **UMA** (Meta) by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or DFT (`dft`). For a single path driven by one or more coordinates, use [`scan`](scan.md); for three coordinates, use [`scan3d`](scan3d.md).

---

## Examples

Here `pocket.pdb` contains the full system matching `real.parm7`, and `ml_region.pdb` selects its ML region without link hydrogens.

The atoms are those of the bundled enzyme example (`examples/beza/1.R.pdb`), whose PDB has an empty chain column. Each atom is therefore written with three fields, residue name, residue number, and atom name, in any order and separated by commas or spaces (`"SAM,320,CS1"`).

### 1. From a YAML spec

Write the two ranges under `pairs:`. With the default step of 0.2 Å, this file gives a 9 × 9 grid.

```yaml
# scan2d.yaml
pairs:
  - ["SAM,320,CS1", "GPP,321,C7", 1.50, 3.00]
  - ["GPP,321,H11", "GLU,186,OE2", 0.90, 2.50]
```

Give the structure, the charge, and the file; add `--out-json` to also get `result.json`.

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan2d.yaml --out-json -o ./result_scan2d/
```

Open `result_scan2d/scan2d_map.png` for the contour map and `scan2d_landscape.html` for the 3D surface; `result.json` gives `scientific_status` and the number of usable points.

### 2. Inline literal

The same two ranges can be written on the command line as one literal.

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50)]'
```

### 3. Pre-optimize, dump, and set the baseline

Optimize the input before the scan, keep the inner-loop trajectories, and measure the relative energies from the lowest usable point. `--max-step-size 0.20`, `--opt-mode grad`, and `--baseline min` are the defaults, written out here.

```bash
mlmm scan2d -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50)]' \
    --max-step-size 0.20 --dump -o ./result_scan2d/ --opt-mode grad \
    --preopt --baseline min
```

---

## How it works

1. **Starting structure and grid**:
The ML-region charge comes from `-q` or `-l` (see {ref}`Charge specification <charge-specification>`). With `--preopt`, the input is first optimized without restraints; if that does not converge, the input geometry is used. Each axis gets ceil(|high − low| / h) + 1 evenly spaced values, both ends included, where h is `--max-step-size` (Å) for a distance and `--max-angle-step-size` or `--max-dihedral-step-size` (degrees) for an angle or dihedral. The values are visited from the one closest to the starting structure.
2. **Outer and inner loops**:
For each d₁ value, the structure is relaxed with only the d₁ restraint. The inner loop then scans d₂ with both restraints, starting each point from the nearest point that has already converged.
3. **Relaxation at each point**:
A harmonic restraint E = ½ k (q − q_target)² holds each coordinate q at its target (k is `--restraint-k`), and L-BFGS (`--opt-mode grad`, default) or RFO (rational function optimization, `--opt-mode hess`) relaxes the rest with the ML/MM calculator; the atoms of the frozen MM layer stay fixed. The energy is then computed without the restraints, and the structure is written under `grid/`.
4. **Table and plots**:
After the last point, all points go into `surface.csv`. The usable points are interpolated with radial basis functions (RBF) on a 50 × 50 grid and drawn as a contour map (`scan2d_map.png`) and a 3D surface (`scan2d_landscape.html`).

---

## Reading surface.csv

`surface.csv` has one row per grid point and one reference row.

| Column | Meaning |
| --- | --- |
| `i`, `j` | Grid indices. Index 0 is the value closest to the starting structure, so the indices follow the visiting order, not ascending values |
| `d1_A`, `d2_A` (also `q1`, `q2`) | Coordinate values measured after the relaxation. The `_A` names are kept for every axis, so an angle axis holds degrees; `q1_unit` and `q2_unit` give the unit (`angstrom` or `degree`) |
| `target_d1_A`, `target_d2_A` (also `target_q1`, `target_q2`) | Restraint targets of the point |
| `energy_hartree` | Energy without the restraints (Hartree) |
| `energy_kcal` | Energy relative to the baseline (kcal/mol) |
| `bias_converged` | Whether the restrained relaxation converged |
| `is_preopt` | `true` only for the reference row |
| `d1_label`, `d2_label` | Axis labels used in the plots |

* **Reference row**: `i = j = -1` and `is_preopt = true` hold the starting structure. The row stays in the table but is never a grid point, a baseline, or a plotted point.
* **Usable points**: a point is usable when its relaxation converged, its energy is finite, and its structure was written. Only usable points set the baseline and enter the plots.
* **Verdict**: in `result.json` (`--out-json`), `scientific_status` is `success` when every grid point is usable, `partial` when only some are (exit status 0), and `failed` when none is (exit status 1); `n_points_attempted` and `n_points_usable` give the counts (see {ref}`Exit codes <exit-codes>`).
* **Next step**: the plots are interpolated, so take a structure from a computed point near the saddle, `grid/point_*.pdb`, and pass it to [`tsopt`](tsopt.md); points in the two basins can be inputs for [`path-search`](path-search.md).

---

## Output files

After a run, `--out-dir` (default: `./result_scan2d/`) contains the following files.

```text
result_scan2d/
├─ surface.csv                   # grid table with the reference row
├─ scan2d_map.png                # 2D contour map
├─ scan2d_landscape.html         # 3D surface (open in a browser)
├─ grid/
│  ├─ point_i150_j090.xyz        # relaxed structure of each grid point
│  ├─ preopt_iDDD_jDDD.xyz       # starting structure (the reference row)
│  └─ inner_path_d1_000_trj.xyz  # inner-loop trajectory for each d₁ value (--dump)
└─ result.json                   # summary (--out-json); summary.json holds the same content
```

Start with `surface.csv` and the two plots; the structure of each point is under `grid/`. In `result.json`, `grid_points[]` maps each grid index to its values, targets, energy, convergence, and structure file; use it instead of reading values back from file names.

* **File names**: the number after `i` or `j` (the tag `DDD`) is the target × 100 (Å, or degrees for an angle), padded to at least three digits, not the grid index of `surface.csv`: `d1 = 1.50 Å, d2 = 0.90 Å` gives `point_i150_j090.xyz`, and an angle of 120° gives `12000`. When two points round to the same tag, the later file name gets `_grid_III_JJJ` with the zero-based indices.
* **Other formats**: each structure is also written as `.pdb` (`--no-convert-files` turns this off), with the layer of each atom in the B-factor column: 0 for the ML region, 10 for movable MM, and 20 for frozen MM. mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files with the original identifiers (see {ref}`mmCIF input <mmcif-input>`).

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `scan2d`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Total charge of the unknown ligand residues (e.g. `-1`) or per-residue charges (e.g. `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-s, --scan-lists` | text | (required) | Two ranges as a YAML/JSON spec file or one inline literal: distance `(i,j,low,high)`, angle `(i,j,k,low,high)`, or dihedral `(i,j,k,l,low,high)` |
| `-o, --out-dir` | path | `./result_scan2d/` | Output directory |
| `--max-step-size` | float | `0.2` | Largest grid spacing of a distance axis (Å) |
| `--max-angle-step-size` | float | `5.0` | Largest grid spacing of an angle axis (degrees) |
| `--max-dihedral-step-size` | float | `10.0` | Largest grid spacing of a dihedral axis (degrees) |
| `--restraint-k` | float | `300.0` | Restraint strength k (eV/Å² for distances, eV/rad² for angles); alias `--bias-k` |
| `--preopt/--no-preopt` | flag | `False` | Optimize the input without restraints before the scan |
| `--opt-mode` | `grad` / `hess` | `grad` | Relaxation of each point: L-BFGS / RFO |
| `--dump/--no-dump` | flag | `False` | Write the inner-loop trajectory for each d₁ value under `grid/` |
| `--baseline` | `min` / `first` | `min` | Zero of `energy_kcal`: lowest usable point, or point `(0, 0)` |
| `--zmin`, `--zmax` | float | surface min / max | Lower and upper ends of the color scale (kcal/mol) |
| `--thresh` | text | `baker` | Convergence preset of each relaxation (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`) |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/scan2d.md) for every option.

> **Note:** In YAML (`--config`), the `lbfgs` (or `opt.lbfgs`) section applies with `--opt-mode grad` and the `rfo` (or `opt.rfo`) section with `--opt-mode hess`. Trajectories are written only with `--dump`; YAML `opt.dump` is not used.

---

## Notes

* **Two ranges in one literal**: `-s` takes exactly two ranges, in one inline literal or under `pairs:` in a YAML/JSON file. For staged scans, use [`scan`](scan.md).
* **For a PDB with chains**, the positional form `A:SAM:320:CS1` picks one atom without ambiguity (see {ref}`Scan-list spec <scan-list-spec>`).
* **Restraint units**: `--restraint-k` is in eV/Å² for distances and eV/rad² for angles and dihedrals. When the option is omitted, YAML `bias.k` applies.
* **Energies without restraints**: every recorded energy is computed after the restraints are removed, so `surface.csv` can go straight into fitting or plotting scripts.
* **Baseline**: `--baseline min` (default) puts zero at the lowest usable point; `--baseline first` puts it at point `(i, j) = (0, 0)`, or at the lowest usable point when `(0, 0)` is not usable.
* **Frozen atoms**: the atoms given by `--freeze-atoms` or YAML `geom.freeze_atoms`, and the atoms of the frozen MM layer, stay fixed in every relaxation. A scanned coordinate whose atoms are all frozen is an error (see {ref}`Frozen atoms and restraints <freeze-atoms-and-restraints>`).
* **Check the spec without computing**: `--dry-run` reads the input, the charge and spin, and `-s`, prints the plan, and exits without any optimization.
* **Cycle limit**: `--relax-max-cycles` (default `100000`) limits each relaxation; an explicit value overrides YAML `opt.max_cycles`.
* **Too few usable points**: with fewer than three usable points, or with all of them on one line, only the plots are skipped. The run prints `[plot] NOTE: Plots skipped: …`, still writes `surface.csv` (and `result.json` with `--out-json`), and exits with status 0. With no usable point it prints `[plot] No finite data for plotting.` and exits with status 1.
* **PNG export**: the PNG is written with Plotly and Kaleido. If the export fails, the run prints `[plot] NOTE: PNG export skipped: …`, still writes the HTML surface, and leaves the PNG out of `result.json`.

---

## See also

* [scan](scan.md) — staged scans of one or more coordinates from one structure
* [scan3d](scan3d.md) — energy grid over three coordinates
* [opt](opt.md) — single-structure optimization before or after a scan
* [tsopt](tsopt.md) — optimize a structure near the saddle into a TS
* [path-search](path-search.md) — MEP through structures taken from the map
* [all](all.md) — end-to-end workflow
* [Troubleshooting](troubleshooting.md) — diagnosing failed runs
