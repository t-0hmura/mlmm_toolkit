# `mlmm scan`, `mlmm scan2d`, and `mlmm scan3d`

Drive distances, angles, or dihedrals with harmonic restraints, relaxing the
rest of the ML/MM system at every step. Run
`mlmm scan -i 1.R.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' -s scan.yaml --out-json`.
Success is `[stage k] Covalent-bond changes (start vs final): Yes` for the
bonds you meant, `====== Scan finished ======`, and `scientific_status`
`success` in `result.json` (for grids, enough usable points in
`surface.csv`).

## When to use

- `scan`: staged 1D scans with relaxation between stages. It makes a
  candidate path from one structure to seed `path-search`, or a TS
  candidate. Inside a full run, `mlmm all -s` does this
  ([all-scan-list.md](all-scan-list.md)); run `scan` alone for one-off
  exploration or other step sizes and restraints.
- `scan2d`: a grid over two coordinates driven together, to map concerted
  versus stepwise surfaces (for example nucleophilic attack and
  leaving-group departure).
- `scan3d`: a grid over three coordinates. Rare: a sequence of 1D or 2D
  scans usually captures the chemistry at much lower cost. Use it for
  inherently 3D-coupled mechanisms, such as two proton transfers coupled to
  one donor distance.

## Writing -s

- `-s` takes inline Python literals or a YAML/JSON file (`stages:` for
  `scan`, `pairs:` for `scan2d` and `scan3d`). Wrap an inline literal in
  single quotes and the atom selectors in double quotes.
- An atom is a 1-based index or a selector. With chain IDs, use
  `"A:SAM:320:CS1"` (chain, residue name, number, atom). Without a chain,
  as in the bundled example, give three fields (residue name, number, atom
  name) in any order, separated by spaces, commas, colons, slashes,
  backticks, or backslashes: `"SAM 320 CS1"` and `"CS1,SAM,320"` are the
  same atom.
- `scan`: `(i, j, target)` drives a distance to a target. Tuples in one
  literal move together as one stage; several literals after one `-s` run as
  stages in sequence, each starting from the previous stage's final
  geometry. A range `(i, j, low, high)` scans both ways from the input and
  gives two stages, toward `low` and then `high`; angles
  `(i, j, k, low, high)` and dihedrals `(i, j, k, l, low, high)` work the
  same way. Inline, targets and ranges cannot be mixed; use YAML `stages:`.
- `scan2d` and `scan3d`: exactly two or three ranges in one literal (or
  under `pairs:`). The tuples are the grid axes.
- Staged versus concerted:
  [Staged vs concerted scan](../mlmm-overview/ts-strategy.md#5-staged-vs-concerted-scan).

## Minimal run

```bash
mlmm scan -i 1.R.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' -b uma -o result_scan
mlmm scan -i 1.R.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60)]' '[("GPP 321 H11","GLU 186 OE2",0.90)]' \
    -b uma -o result_scan_staged
```

```bash
mlmm scan2d -i 1.R.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.60,3.10), ("GPP 321 H11","GLU 186 OE2",0.90,1.80)]' \
    -b uma -o result_scan2d
mlmm scan3d -i 1.R.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' \
    -s '[("SAM,320,CS1","GPP,321,C7",1.50,3.00),("GPP,321,H11","GLU,186,OE2",0.90,2.50),("SAM,320,SD","SAM,320,CS1",1.80,3.00)]' \
    -o result_scan3d
mlmm scan3d --csv result_scan3d/surface.csv -o result_scan3d_plot
```

The last line redraws a finished 3D grid without a structure or topology.
`-b` selects `uma`, `orb`, `mace`, `aimnet2`, or `dft`. The step is
`--max-step-size` (0.2 Å) for distances, `--max-angle-step-size` (5°), and
`--max-dihedral-step-size` (10°). `--dry-run` reads the input, the charge,
and `-s`, prints the plan, and stops.

## Judge success

**scan.** Each stage prints
`[stage k] Covalent-bond changes (start vs final): Yes` (or `No`), and the
run ends with a `Summary` and `====== Scan finished ======`. With `--out-json`, `result.json` gives
`scientific_status`: `success` when every step converged (and `--preopt`
and `--endopt`, if requested), `partial` when some did (exit 0), `failed`
when none did (exit 1). Per stage it has `stages[].converged`,
`stages[].bond_changes`, `stages[].final_energy_hartree`, and
`stages[].energies_hartree`. Files: `stage_NN/result.{xyz,pdb}` (end of each
stage), `stage_NN/scan_trj.xyz`, and `scan_trj.xyz` and `scan.pdb` joining
all stages; `preopt/` with `--preopt`. Each frame's comment line holds the
energy without restraints; plot it with [trj2fig](utilities.md#trj2fig).
The highest step is a TS candidate, and for `path-search` give the start
followed by the `stage_NN/result.*` files in order.

**scan2d and scan3d.** `surface.csv` has one row per grid point and one
reference row (`i = j = -1`, `is_preopt` true) for the start. A point is
usable when `bias_converged` is true, its energy is finite, and its
structure was written. `scan2d` also writes `scan2d_map.png` (contour) and
`scan2d_landscape.html` (3D surface); `scan3d` writes
`scan3d_density.html` (isosurfaces). With `--out-json`, `result.json` has
`scientific_status` (`success` when every point is usable, `partial` when
some are, exit 0; `failed` when none is, exit 1), `n_points_attempted`, and
`n_points_usable`. Structures are `grid/point_i150_j090.xyz` (and `.pdb`):
the tag is the target × 100, not the grid index; `grid/preopt_*.xyz` is the
start, and `--dump` adds `grid/inner_path_d1_*_trj.xyz`. The plots are
interpolated, so take a computed `grid/point_*.pdb` near the saddle for
`tsopt`.

## Pitfalls and recovery

- Stage k+1 starts from stage k, so one diverged stage derails every later
  stage. Give all stage literals after a single `-s`; do not repeat the
  flag.
- Target distances must be positive, and a coordinate may appear only once
  per stage. A coordinate whose atoms are all frozen is an error.
- A scan from the product gives the reverse barrier, E(TS) − E(P); the
  forward barrier needs E(R) from an optimized reactant:
  [Reading the barrier](../mlmm-overview/ts-strategy.md#4-reading-the-barrier-when-the-scan-started-from-p).
- Standalone `scan` does not pre-optimize unless you pass `--preopt`;
  inside `all`, `--preopt` is on by default.
- Grid cost is the product of the axis lengths: 10 × 10 = 100
  relaxations, 5 × 5 × 5 = 125, and 9 × 9 × 7 = 567. Start with a larger
  `--max-step-size` or narrower ranges; for 1D or 2D chemistry, `scan` and
  `scan2d` are far cheaper than `scan3d`.
- Too few usable points (under three for `scan2d`, four for `scan3d`, or
  all on one line or plane): only the figure is skipped, with
  `[plot] NOTE: Plots skipped: …` or `[plot] NOTE: Volume plot skipped: …`;
  `surface.csv` is still written and the exit code is 0.
- `--baseline min` (default) puts zero at the lowest usable point;
  `--baseline first` at grid point 0, or at the lowest usable point when
  that one is not usable.
- `scan3d --csv` needs the columns `d1_A`, `d2_A`, `d3_A`, and
  `energy_hartree` or `energy_kcal`. A `surface.csv` from `scan3d` lacks
  `artifact_written`, so re-plotting prints
  `[plot] WARNING: CSV lacks complete point provenance; …`. Redrawing into
  the scan's own directory replaces or removes its `result.json` and old
  `scan3d_density.html`; give another `-o`.

## Next step

- [path.md](path.md): MEP search from the scan results.
- [tsopt.md](tsopt.md): optimize the highest step or a grid point near the
  saddle.
- [all-scan-list.md](all-scan-list.md): staged scans inside the full
  workflow.
- [trj2fig](utilities.md#trj2fig): plot `scan_trj.xyz`.
