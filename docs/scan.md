# `scan` (restrained coordinate scan)

`scan` drives chosen distances, angles, or dihedrals of a layered enzyme structure step by step with harmonic restraints, relaxing every other degree of freedom with the ML/MM calculator at each step, and so builds a candidate reaction path from a single structure. The coordinates in one literal (or the distance targets of one YAML stage) move together as one **stage**; several literals run as stages in sequence, each starting from the relaxed end of the previous one.

---

## What it is for

* **A path from one structure**: drive the reacting bonds of a reactant to get intermediate- and product-like structures for [`path-search`](path-search.md).
* **Testing the order of events**: drive bond formation and proton transfer in one stage or in separate stages, and compare the energy profiles.
* **Running the scan step of `all` on its own**: repeat the scan that [`all`](all.md) runs for `-s`, with other step sizes or restraints.
* **Judging the result**: each stage reports whether covalent bonds formed or broke, and `result.json` gives `scientific_status`. See {ref}`Reading the result <scan-checking-result>`.

The ML region uses **UMA** (Meta) by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or DFT (`dft`). For an energy grid over two or three independent coordinates, use [`scan2d`](scan2d.md) or [`scan3d`](scan3d.md).

---

## Examples

Here `complex.pdb` contains the full system matching `real.parm7`, and `ml_region.pdb` selects its ML region without link hydrogens.

The atoms are those of the bundled enzyme example (`examples/beza/1.R.pdb`), whose PDB has an empty chain column. Each atom is therefore written with three fields, residue name, residue number, and atom name, in any order and separated by commas or spaces.

### 1. From a YAML spec

Write the stages in a file and add `--out-json` for a summary.

```yaml
# scan.yaml
stages:
  - [["SAM,320,CS1", "GPP,321,C7", 1.60]]
  - [["GPP,321,H11", "GLU,186,OE2", 0.90]]
```

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan.yaml --out-json -o ./result_scan
```

For each stage the console prints `[stage k] Covalent-bond changes (start vs final): Yes` (or `No`), and the run ends with a `Summary` of every stage and `====== Scan finished ======`. In `result_scan/result.json`, `scientific_status` is `success` when every stage converged.

### 2. Inline literal

A short single-stage scan can be given on the command line.

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]'
```

### 3. Two coordinates in one stage

Coordinates in the same literal move together (a concerted step).

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("CS1 SAM 320","GPP 321 C7",1.60),("GPP 321 H11","GLU 186 OE2",0.90)]' -o ./result_concerted
```

### 4. Two stages in sequence

Give several literals after one `-s`; stage 2 starts from the relaxed result of stage 1.

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[("SAM,320,CS1","GPP,321,C7",1.60)]' '[("GPP,321,H11","GLU,186,OE2",0.90)]' -o ./result_staged
```

### 5. Bidirectional scan

A [4-tuple](#bidirectional-scan-4-tuple) scans one distance in both directions from the input geometry.

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s '[(12, 45, 1.35, 2.50)]'
```

### 6. Dump trajectories

Add `--dump` to keep the optimizer trajectory of every step.

```bash
mlmm scan -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    -s scan.yaml --dump -o ./result_scan_dump
```

---

## How it works

1. **Reading the structure**: the ML-region {ref}`charge <charge-specification>` comes from `-q` or `-l`. With `--preopt`, the structure is first optimized without restraints; if that does not converge, the input geometry is used.
2. **Splitting each stage into steps**: for every coordinate, `scan` takes the change Δ = target − current and divides the stage into N = ceil(max(|Δ| / h)) steps, where h is `--max-step-size` (Å) for distances, `--max-angle-step-size` for angles, and `--max-dihedral-step-size` for dihedrals (degrees). Each coordinate moves by Δ / N per step, so all coordinates of a stage arrive together.
3. **Restrained relaxation**: at each step, a harmonic restraint E = ½ k (q − q_target)² holds every scanned coordinate q at its step target (k from `--restraint-k`), and the rest of the structure is relaxed with the ML/MM calculator by L-BFGS (`--opt-mode grad`, default) or RFO (rational function optimization, `--opt-mode hess`); the atoms of the frozen MM layer stay fixed. The energy written for each step is the ML/MM energy computed with the restraints removed. Cartesian coordinates (`geom.coord_type: cart`) are the default and recommended for ML/MM scans; `dlc` can be selected in YAML but can take much longer to converge.
4. **End of the stage**: with `--endopt`, the last structure of the stage is optimized once more without restraints. `scan` then compares the first and last structures of the stage for covalent-bond changes and writes the stage result.
5. **Next stage**: the next stage starts from this result. After the last stage, the trajectories of all stages are joined into one file.

### Bidirectional scan (4-tuple)

A range `(i, j, low, high)` instead of a target `(i, j, target)` scans in both directions from the input geometry. It expands into two stages:

1. **Pass 1**: drive `i`–`j` from the current distance toward `low`.
2. **Pass 2**: restore the input geometry and drive `i`–`j` toward `high`.

The joined trajectory runs `low → input geometry → high`, a continuous path through the starting structure. Angle ranges `(i, j, k, low, high)` and dihedral ranges `(i, j, k, l, low, high)` are scanned the same way.

(section-bond)=
### Bond-change detection

Let T be the sum of the covalent radii of two atoms scaled by `bond_factor` (default `1.20`). The atoms count as bonded when their distance is at most T − `margin_fraction` × T (`margin_fraction` defaults to `0.05`). A pair is reported as formed or broken only when its distance changed by at least `delta_fraction` × T (`delta_fraction` defaults to `0.05`). `path-search` uses the same rules; the keys are in the YAML [`bond`](yaml-reference.md#bond) section.

---

(scan-checking-result)=
## Reading the result

| Where | What to check |
| --- | --- |
| Console, each stage | `[stage k] Covalent-bond changes (start vs final): Yes` with the formed and broken bonds listed, or `No` with `(no covalent changes detected)` |
| Console, end of run | `Summary`: targets, initial values, per-coordinate step, number of steps, and bond changes of each stage, followed by `====== Scan finished ======` |
| `result.json` (`--out-json`) | `scientific_status`: `success` when every step of every stage converged (and the `--preopt` and `--endopt` optimizations, when requested) with a finite energy; `partial` when only some of them did; `failed` when none did |
| `result.json` (`--out-json`) | `stages[].converged`, `stages[].bond_changes.changed`, `stages[].final_energy_hartree`, and the energy of every step in `stages[].energies_hartree` |

A `partial` run exits with 0 and a `failed` run with 1. See {ref}`Exit codes <exit-codes>`. A converged scan with the intended bond changes gives a candidate path; the highest-energy step is a TS candidate for [`tsopt`](tsopt.md), which you can {ref}`extract <trajectory-one-frame>` from `scan_trj.xyz`.

---

(scan-direction-and-barrier-sign)=
## Barrier sign

`scan` records energies but does not report a barrier. If you read a barrier off a scan (or a path, or a TS candidate made from one) that **started from the product**, the difference from the starting structure is the **reverse** barrier, `E(TS) − E(product)`. The forward barrier is computed from the reactant:

| You ran | Forward barrier |
| --- | --- |
| A scan from the reactant | `E(TS) − E(reactant)`, the difference from the starting structure |
| A scan from the product | `E(TS) − E(reactant)`, **not** the difference from the starting structure; E(reactant) comes from an optimized reactant, for example an IRC endpoint optimized with [`opt`](opt.md) |

No option changes this. Before quoting a barrier, check which endpoint the scan started from, especially when the starting structure was a crystallographic product complex.

---

## Output files

`scan` writes these files to `--out-dir`:

```text
result_scan/
├─ preopt/
│  └─ result.xyz                    # Pre-optimized structure (--preopt)
├─ stage_01/                        # One directory per stage (stage_NN)
│  ├─ result.xyz                    # Final geometry of the stage
│  ├─ scan_trj.xyz                  # Structure and energy of every step in the stage
│  └─ scan_s0001_optimization_trj.xyz  # Optimizer trajectory of each step (--dump)
├─ scan_trj.xyz                     # All stages joined
└─ result.json                      # Summary (--out-json); summary.json has the same content
```

The structures and trajectories are also written as PDB under the same names (`result.pdb`, `scan.pdb`); `--no-convert-files` turns this off. {ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Stage results**: `stage_NN/result.*` is the structure at the end of stage NN. For [`path-search`](path-search.md), give the starting structure followed by the `stage_NN/result.*` files in stage order.
* **Energy profile**: the comment line of each frame in `scan_trj.xyz` holds the energy without restraints (Hartree); plot it with [`trj2fig`](trj2fig.md).

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `scan`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Total charge of the unknown ligand residues (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-s, --scan-lists` | text | (required) | A YAML/JSON spec file, or one or more inline literals (one per stage): distance targets `(i,j,target)`, or ranges for a distance `(i,j,low,high)`, an angle `(i,j,k,low,high)`, or a dihedral `(i,j,k,l,low,high)` |
| `-o, --out-dir` | path | `./result_scan/` | Output directory |
| `--one-based/--zero-based` | flag | `--one-based` | Read atom indices in `-s` as 1-based or 0-based |
| `--max-step-size` | float | `0.2` | Largest change of a distance per step (Å) |
| `--max-angle-step-size` | float | `5.0` | Largest change of an angle per step (degrees) |
| `--max-dihedral-step-size` | float | `10.0` | Largest change of a dihedral per step (degrees) |
| `--restraint-k` | float | `300.0` | Restraint strength k (eV/Å² for distances, eV/rad² for angles); alias `--bias-k` |
| `--preopt/--no-preopt` | flag | `False` | Optimize the input structure without restraints before the scan |
| `--endopt/--no-endopt` | flag | `False` | Optimize the result of each stage without restraints |
| `--dump/--no-dump` | flag | `False` | Write the optimizer trajectory of every step |
| `--opt-mode` | `grad` / `hess` | `grad` | Relaxation: L-BFGS / RFO (on `tsopt` the same words select other optimizers; see {ref}`--opt-mode by command <opt-mode-semantics>`) |
| `--freeze-atoms` | text | `None` | Comma-separated 1-based atom indices to freeze, added to YAML `geom.freeze_atoms` and the frozen MM layer |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/scan.md) for every option.

> **Note:** In YAML (`--config`), [`bias.k`](yaml-reference.md#bias) sets the restraint strength when `--restraint-k` is not given, and the [`bond`](yaml-reference.md#bond) section sets the bond-change thresholds `bond_factor`, `margin_fraction`, and `delta_fraction`.

---

## Notes

* **`--preopt` depends on the caller**: run on its own, `scan` does not pre-optimize unless you pass `--preopt`. Inside `all`, the scan pre-optimizes when `all --preopt` is on (the default), and [`all --scan-preopt/--no-scan-preopt`](reference/commands/all.md) overrides it.
* **Targets and ranges are not mixed inline**: one inline literal, and all literals of one run, hold either targets `(i,j,target)` or ranges. To combine them, list them under `stages:` in a YAML/JSON spec.
* **Stage numbers with ranges**: a range gives two stages, toward `low` and then toward `high` (a single 4-tuple gives `stage_01/` and `stage_02/`). Inline, all ranges of one literal move together in these two stages; in a YAML `stages:` list, each entry of a stage that holds a range becomes its own stage, one for a target and two for a range.
* **Target distances must be positive**, and one coordinate may appear only once per stage.
* **Check the spec without computing**: `--dry-run` reads the input, the charge and spin, and `-s`, prints the number of stages, and exits without any optimization.
* **Frozen atoms**: the atoms given by `--freeze-atoms` or YAML `geom.freeze_atoms`, and the atoms of the frozen MM layer, stay fixed in every relaxation. A scanned coordinate whose atoms are all {ref}`frozen <freeze-atoms-and-restraints>` is an error.
* **Cycle limit**: `--relax-max-cycles` (default `100000`) limits each relaxation; when given, it overrides YAML `opt.max_cycles`.

---

## See also

* {ref}`Scan-list spec <scan-list-spec>` — YAML/JSON spec files, inline literals, and atom selectors
* [scan2d](scan2d.md) — energy map over two coordinates
* [scan3d](scan3d.md) — energy grid over three coordinates
* [path-search](path-search.md) — minimum energy path (MEP) search from the scan results
* [all](all.md) — the full workflow, including a scan from one structure with `-s`
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
