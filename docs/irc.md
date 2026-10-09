# `irc` (intrinsic reaction coordinate)

`irc` traces the intrinsic reaction coordinate (IRC) of an ML/MM system from an optimized transition state (TS) in both directions with EulerPC (an Euler predictor–corrector integrator). It writes the trajectory of each branch and the two endpoint candidates. Optimizing those endpoints with [`opt`](opt.md) shows which reactant (R) and product (P) the TS connects.

---

## What it is for

* **Checking a TS**: after [`tsopt`](tsopt.md) and [`freq`](freq.md) (n_imag = 1), confirm that the TS connects the intended R and P.
* **Getting R and P**: optimize the endpoints with [`opt`](opt.md) to obtain the R and P structures of this TS.
* **Rerunning the IRC step of `all`**: trace the IRC of an [`all`](all.md) run again on its own, with different settings.

The ML region uses **UMA** (Meta) by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or **DFT**. The MM atoms use the Amber force field of `--parm7`.

---

## Examples

### 1. Both branches

Trace both directions from the TS `ts.pdb`, and write a summary with `--out-json`.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --out-json --out-dir ./result_irc
```

The endpoint candidates are `forward_first.xyz` and `backward_last.xyz`.

### 2. Forward branch only

Trace only the forward branch.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --no-backward --out-dir ./result_irc_forward
```

### 3. Analytical Hessian

Have the ML backend compute the starting Hessian analytically instead of by finite differences.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

### 4. Retry with a smaller step

When a branch stops after only a few frames, retry with a maximum step of 0.05 bohr.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --step-size 0.05 --out-dir ./result_irc_small_step
```

### 5. Trace to the cycle limit

Add `--never-stop` to ignore the gradient and energy stop criteria and trace each branch until `--max-cycles`.

```bash
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --step-size 0.05 --never-stop --max-cycles 250 \
    --out-dir ./result_irc_continue
```

---

## How it works

1. **Building the ML/MM system**: `irc` reads the TS structure from `-i`, the Amber topology from `--parm7`, and the {ref}`ML region <mlmm-options>` from `--model-pdb`. `-q` and `-m` are the charge and the spin multiplicity of the ML region.
2. **Starting direction**: `irc` computes the Hessian at the TS (or reads it with `--read-hess`), removes rigid motions as [`freq`](freq.md#rigid-modes-with-frozen-boundaries) does, and takes the eigenvector `--root` (default `0`, the lowest eigenvalue) as the reaction mode. If that mode is not imaginary, the run stops with an error.
3. **EulerPC integration**: each branch (forward, then backward) starts from the TS. Every step is an Euler predictor along the mass-weighted steepest-descent direction, with the gradient estimated from a second-order Taylor expansion with the current Hessian (Bofill update), followed by a modified Bulirsch–Stoer corrector on a DWI (distance-weighted interpolation) surface. A branch stops when the RMS gradient falls below 1 × 10⁻³ hartree/bohr after leaving the TS region, when the energy rises, when the energy changes by 1 × 10⁻⁶ hartree or less in one step, or at `--max-cycles` (default 125 steps per branch).
4. **Writing the path**: `irc` writes each branch, the whole path through the TS, and the end structures. For PDB/mmCIF input or with `--ref-pdb`, the trajectories and the two endpoint candidates are also converted to PDB.

---

## Judging the IRC

Even if the IRC does not converge, the result is usable when the endpoints, optimized with [`opt`](opt.md), reach the intended R and P.

| What to check | Where to look |
| --- | --- |
| The start was a TS | The console line `Transition vector is mode 0 with wavenumber … cm⁻¹.` shows a negative wavenumber |
| How each branch stopped | `forward_integration_converged` / `backward_integration_converged` in `result.json`: `true` when the RMS gradient fell below the threshold, `false` for an energy stop or the cycle limit; `forward_integration_stop_reason` / `backward_integration_stop_reason` give the reason |
| Bonds that change along the path | `bond_changes` in `result.json` (`formed` and `broken`, from `finished_first` to `finished_last`) |
| Which end is R and which is P | Optimize `forward_first.xyz` and `backward_last.xyz` with [`opt`](opt.md) and compare them with the intended R and P. The direction forward / backward does not decide it |

`irc` does not judge the endpoints; whether they are the intended R and P is for you to check.

Optimize both endpoints with `opt`; they are `.xyz` files, so pass the TS PDB with `--ref-pdb` for the atom order and the layers:

```bash
mlmm opt -i result_irc/forward_first.xyz --ref-pdb ts.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -q 0 -m 1 --out-dir ./result_opt_forward
mlmm opt -i result_irc/backward_last.xyz --ref-pdb ts.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -q 0 -m 1 --out-dir ./result_opt_backward
```

If the endpoints are not the intended R and P, see {ref}`When the TS search fails <ts-search-fails>`.

---

## Output files

When the run finishes, `--out-dir` (default `./result_irc/`) contains:

```text
result_irc/
├─ finished_irc_trj.xyz    # Whole IRC path through the TS
├─ finished_irc.pdb        # Same path as PDB
├─ finished_first.xyz      # First frame of the whole path (= forward_first.xyz when the forward branch runs)
├─ finished_last.xyz       # Last frame of the whole path (= backward_last.xyz when the backward branch runs)
├─ forward_irc_trj.xyz     # Forward branch, from the TS (when it runs)
├─ forward_irc.pdb         # Same branch as PDB
├─ forward_first.xyz       # End of the forward branch (endpoint candidate)
├─ forward_first.pdb       # Same structure as PDB
├─ backward_irc_trj.xyz    # Backward branch, from the TS (when it runs)
├─ backward_irc.pdb        # Same branch as PDB
├─ backward_last.xyz       # End of the backward branch (endpoint candidate)
├─ backward_last.pdb       # Same structure as PDB
└─ result.json             # Summary (--out-json)
```

The `.pdb` files are written for PDB/mmCIF input or with `--ref-pdb`. {ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Endpoint candidates**: `forward_first.xyz` and `backward_last.xyz` are the structures to optimize with [`opt`](opt.md). Each branch also writes its other end, next to the TS (`forward_last.xyz`, `backward_first.xyz`).
* **Path**: open `finished_irc_trj.xyz` or `finished_irc.pdb` in PyMOL or VMD to watch the reaction.
* **Summary**: with `--out-json`, [`result.json`](json-output.md) records the number of frames of each branch (`n_frames_forward`, `n_frames_backward`), how each branch stopped, `bond_changes`, the energies of the two ends and the TS (`energy_first_hartree`, `energy_ts_hartree`, `energy_last_hartree`), and the removed rigid motions and the starting Hessian under `rigid_projection`.
* **Console**: the step table of each branch and the elapsed time.

> **Note:** with YAML `irc.prefix: trial`, every file name other than `result.json` starts with `trial_`, as in `trial_finished_irc_trj.xyz`, and `files` in `result.json` records the prefixed names. A positive YAML `irc.dump_every` also writes the HDF5 checkpoint `irc_data.h5` during the run; it is off by default.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `irc`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | TS structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-l, --ligand-charge` | text | `None` | Total charge of unknown ligand residues (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `--max-cycles` | integer | `125` | Maximum number of IRC steps per branch |
| `--step-size` | float | `0.10` | Maximum step length in bohr (unweighted Cartesian coordinates) |
| `--root` | integer | `0` | Hessian eigenvector used as the reaction mode, counted from 0 in ascending order of eigenvalue |
| `--forward/--no-forward` | flag | `True` | Run the forward branch |
| `--backward/--no-backward` | flag | `True` | Run the backward branch |
| `--never-stop/--no-never-stop` | flag | `False` | Ignore the gradient and energy stop criteria and trace until `--max-cycles` |
| `-o, --out-dir` | path | `./result_irc/` | Output directory |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | How the ML backend computes the starting Hessian |
| `-b, --backend` | text | `uma` | ML-region backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--read-hess` | path | `None` | Start from the Hessian in a `.npy` file (for example from `freq` or `tsopt --dump-hess`) instead of computing it |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/irc.md) for every option.

> **Note:** in YAML (`--config`), every key of the `irc` block is listed under [`irc`](yaml-reference.md#irc-section) in the YAML Reference.

---

## Notes

* **A branch that stops at once**: when a branch ends after three frames or fewer before the cycle limit, the console warns `[irc] IRC stopped after only a few frames in …`. Retry with a smaller `--step-size`, for example `0.05`, before changing other settings; a large step can make EulerPC unstable.
* **`--never-stop` is off by default**: with it, the branches run past the physical end points to the cycle limit. Numerical failures and interruptions still stop the run. Inspect the trajectory and optimize the endpoints, and raise `--max-cycles` only when the extra path is useful.
* **`--root` counts from 0**: a successful TS optimization gives one imaginary mode along the reaction coordinate, so for a TS with n_imag = 1 keep `--root 0` (the only negative eigenvalue). Use `1`, `2`, … only when you know that spurious modes with lower (more negative) eigenvalues come before the reaction mode.
* **Cartesian coordinates**: `irc` always uses Cartesian coordinates, whatever YAML `geom.coord_type` says.
* **The `--read-hess` file** is the same `.npy` file as in [`freq`](freq.md), in Hartree/bohr², for all atoms or only the atoms in the Hessian calculation; pass a Hessian computed for the same geometry, charge, multiplicity, and calculator. It needs `irc.hessian_init: calc` (the default); when the file is used, `result.json["rigid_projection"]["hessian_source"]` is `"file"`.
* **Analytical Hessian and `--uma-workers`**: with UMA, `--hessian-calc-mode Analytical` cannot be combined with `--uma-workers` above 1 and stops with an error. Use `--uma-workers 1` for an [analytical Hessian](backends.md). Its speed and memory use depend on the backend and the system, so compare both modes on your system first.
* **Frozen atoms**: besides the frozen MM layer, `--freeze-atoms` freezes more atoms (1-based); how to choose them is described in {ref}`Frozen atoms and distance restraints <freeze-atoms-and-restraints>`.
* **Large systems**: `--hess-device cpu` keeps the starting Hessian and the IRC Hessian operations on the CPU, to stay within GPU memory.
* **At least one branch**: `--no-forward` together with `--no-backward` stops with an error.
* **One structure per run**: `-i` takes a single structure. Extract the frame you need from a trajectory to `.xyz` first and pass it with `--ref-pdb`.

---

## See also

* [tsopt](tsopt.md) — optimize the TS before running IRC
* [freq](freq.md) — check that the TS has one imaginary mode (n_imag = 1)
* [opt](opt.md) — optimize the IRC endpoints to R and P
* [all](all.md) — the full workflow, which runs IRC after `tsopt` and optimizes the endpoints
* [Troubleshooting](troubleshooting.md) — when a run fails
* [YAML Reference](yaml-reference.md) — every `irc` setting
* [Glossary](glossary.md) — IRC and other terms
* {ref}`Exit codes <exit-codes>` — what each exit status means
