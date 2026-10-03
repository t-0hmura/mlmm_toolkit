# `tsopt` (transition-state optimization)

## Overview

`tsopt` optimizes a transition-state (TS) candidate on a layered ML/MM enzyme model to a first-order saddle point, then computes the Hessian at the final geometry and counts its imaginary frequencies (n_imag). A successful TS optimization gives one imaginary mode along the reaction coordinate.

### What it is for

* **Refining a TS candidate**: turn the highest-energy image (HEI) of [`path-opt`](path-opt.md) / [`path-search`](path-search.md), or the top of a [`scan`](scan.md), into an optimized TS.
* **Checking a structure you built**: confirm that a hand-made candidate is a TS (n_imag = 1) and watch its reaction mode as an animation.
* **Rerunning the TS step of `all`**: optimize the TS of an [`all`](all.md) run again on its own, with different settings.

The ML region uses **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md), by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or **DFT**. The MM region uses the Amber parameters in `--parm7`.

If you do not have a candidate yet, make one first:

| What you have | Command that gives a candidate |
| --- | --- |
| The reactant **and** the product | [`path-opt`](path-opt.md) (two structures; `hei.xyz`) or [`path-search`](path-search.md) (two or more structures; one `hei_seg_NN.xyz` per segment with bond changes) |
| Only the reactant, or a bond you want to drive | [`scan`](scan.md) drives the reacting distance step by step and relaxes everything else |

Optimize the candidate with `tsopt`, then follow the TS with [`irc`](irc.md). The [`all`](all.md) command runs all of these steps in one go.

---

## Examples

In these examples, `ts_guess.pdb` is the full-system candidate that matches `real.parm7`, and `ml_region.pdb` defines the ML region (see [Building the ML region and layers](model-setup.md)).

### 1. Default run (RS-P-RFO)

Give the charge and the spin multiplicity of the ML region explicitly.

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_tsopt
```

### 2. Dimer method

Use this when computing the full Hessian again and again is too expensive, or to try a second method on a difficult candidate.

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode dimer --out-dir ./result_tsopt_dimer
```

### 3. Remove extra imaginary modes

Add `--flatten` when the candidate has more than one imaginary frequency.

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --flatten --out-dir ./result_tsopt_flatten
```

### 4. Start from a saved Hessian

Read a Hessian saved with `--dump-hess` at this same geometry (for example by `freq` on `ts_guess.pdb`), instead of computing it again.

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --read-hess ts_guess_hess.npy --out-dir ./result_tsopt
```

---

## How it works

1. **Reading the model and freezing the boundary**: the ML region comes from `--model-pdb` and the movable and frozen MM layers from the B-factors of the input PDB; the charge of the ML region comes from `-q` or `-l` (see {ref}`Charge specification <charge-specification>`). The frozen MM layer and the atoms in `--freeze-atoms` stay fixed, and the final Hessian covers the atoms chosen by `--active-dof-mode` (PHVA: partial Hessian vibrational analysis). An `.xyz` candidate needs the full-system PDB with `--ref-pdb` for the atom order and the layers.
2. **Choosing the optimizer** (`--opt-mode`): `hess` (default) runs **RS-P-RFO** (restricted-step partitioned rational function optimization), which uses the full Hessian; `rsirfo` and `trim` select RS-I-RFO (restricted-step image RFO) and TRIM (trust-region image minimization). `dimer` (or `grad`) runs the **Hessian-guided Dimer** method, which follows the lowest mode with gradients and refreshes its direction from an exact Hessian at intervals.
3. **Climbing along the reaction mode**: the optimizer goes uphill along the reaction mode and downhill along every other direction until the convergence criteria (`--thresh`, default `baker`) are met. `baker` needs all five at once (atomic units): max force below 3 × 10⁻⁴, RMS force below 2 × 10⁻⁴, max step below 3 × 10⁻⁴, RMS step below 2 × 10⁻⁴, and an energy change below 10⁻⁶ hartree, all tighter than Gaussian's default (`gau`). RS-P-RFO updates the Hessian with the Bofill formula and keeps each step within a trust radius of 0.1 bohr (`rsirfo.trust_max`). With microiteration (`--microiter`, on by default), each macro step moves the ML region and the boundary MM atoms bonded to it, and the movable MM atoms then relax with L-BFGS. With microiteration, the optimizer's own lines start with `[microiter]`; Dimer and `--no-microiter` print the `[tsopt]` form.
4. **Final check**: after convergence, `tsopt` computes the Hessian at the final geometry, counts n_imag, and writes each imaginary mode as an animation. Frozen atoms are handled as in [`freq`](freq.md#rigid-modes-with-frozen-boundaries).
5. **Removing extra imaginary modes (only with `--flatten`)**: if more than one imaginary mode remains, `tsopt` displaces the structure along the extra modes and optimizes again, until one mode is left or the round limit is reached. In Dimer mode, each round also computes the exact Hessian again, refreshes the dimer direction, and runs a short Dimer + L-BFGS segment.

---

## Reading the TS result

How the run ended decides what you get.

| How it ended | Console line | What `tsopt` leaves | What `all` does next |
| --- | --- | --- | --- |
| Converged | `[microiter] Converged!` or `[tsopt] Numerical optimization converged.`, then `[tsopt] Wrote N final imaginary mode(s).`; with n_imag ≥ 2, also `[tsopt] WARNING: Higher-order stationary point (n_imag=N). …`; with n_imag = 0, `[tsopt] No imaginary mode detected. …` | Final geometry, n_imag, and the imaginary modes | Continues to IRC if n_imag ≥ 1; stops before IRC if n_imag = 0 |
| Stopped on an energy plateau | `[microiter] Stalled (energy plateau; not converged)` or `[tsopt] Stalled (energy plateau; not converged)` | Final geometry and n_imag | Stops before IRC |
| Reached `--max-cycles` without converging | `[microiter] Reached max macro iterations (M).` or `[tsopt] Reached max cycles (N/M).` | Final geometry; no Hessian | Stops before IRC |
| Converged with `--skip-final-freq` | `[tsopt] WARNING: TS saddle-point order is not verified (--skip-final-freq).` | Final geometry; no Hessian | Stops before IRC, because the reaction mode is unchecked |
| The final Hessian failed | `[tsopt] ERROR: Terminal PHVA failed.` | Final geometry, with `hessian_status: failed` and the reason | Stops before IRC |

Read n_imag as follows:

| n_imag | Meaning |
| --- | --- |
| 1 | A first-order saddle point. Check that the mode moves the intended atoms, then run [`irc`](irc.md). |
| 0 | No imaginary mode: the structure has relaxed toward a minimum. |
| 2 or more | A higher-order saddle point: extra imaginary modes remain. `all` still runs IRC along the imaginary mode that the optimizer followed (the lowest imaginary mode when that one is not available), so the IRC endpoints show where that mode leads. |

A mode counts as imaginary when ν < −5.00 cm⁻¹. Values between −5.00 and 0 cm⁻¹ are treated as numerical noise.

(wrong-imaginary-mode-count)=
### Wrong imaginary-mode count after optimization

If n_imag is not 1, or the mode does not move the atoms of the intended reaction, try the following. The remedies can be combined.

| Result | What to try |
| --- | --- |
| n_imag = 0 | The candidate is not near a saddle. Get a better candidate from a path search or a scan. In `all`, `--refine-path` runs a recursive `path-search` that resolves the HEI more finely; it costs more, because each new step gets its own TS optimization and IRC. |
| n_imag ≥ 2 | Watch each mode. Re-optimize with `--flatten`, or compare `--precision fp32` / `fp64` on this candidate. |
| One mode, but the wrong motion | Check which atoms move, and start from a candidate closer to the intended reaction. |

For example, to retry in fp64 with flattening on:

```bash
mlmm tsopt -i ts_guess.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --precision fp64 --flatten -o result_tsopt
```

For more ideas, see {ref}`When the TS search fails <ts-search-fails>`; for other failures, see [Troubleshooting](troubleshooting.md).

---

## Output files

When the run finishes, `--out-dir` (default `./result_tsopt/`) contains:

```text
result_tsopt/
├─ final_geometry.xyz               # Final geometry (always written)
├─ final_geometry.pdb               # Same, as PDB (B-factors mark the layers)
├─ vib/
│  ├─ imag_01_-385.20cm-1_trj.xyz   # Animation of each imaginary mode
│  └─ imag_01_-385.20cm-1.pdb       # Same, as PDB
├─ optimization_all_trj.xyz         # Optimization trajectory (--dump)
├─ optimization_all.pdb             # Same, as PDB (--dump)
├─ .dimer_mode.dat                  # Current Dimer direction (Dimer only)
└─ result.json                      # Summary (--out-json)
```

mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers (see {ref}`mmCIF input <mmcif-input>`).

* **Final geometry**: `final_geometry.*` is the TS to pass to [`irc`](irc.md). In `final_geometry.pdb`, the B-factor is 0 for the ML region, 10 for movable MM atoms, and 20 for frozen MM atoms.
* **Reaction mode**: open `vib/imag_*_trj.xyz` in PyMOL or VMD and check that the atoms move along the bonds that form or break.
* **Summary**: with `--out-json`, `result.json` records how the run ended (`optimization_status`: `converged`, `stalled`, or `not_converged`), `hessian_status`, and n_imag (see [JSON Output Reference](json-output.md#tsopt)).

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `tsopt`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | One full-system structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`). For a trajectory, extract one frame to `.xyz` first (see {ref}`Extract one frame from a trajectory <trajectory-one-frame>`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Total charge of unknown ligands or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB input or `--ref-pdb`) |
| `-o, --out-dir` | path | `./result_tsopt/` | Output directory |
| `-b, --backend` | text | `uma` | Backend for the ML region (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--opt-mode` | `hess` / `dimer` / `rsirfo` / `trim` | `hess` | Optimizer: RS-P-RFO / Dimer / RS-I-RFO / TRIM (`rsprfo` = `hess`, `grad` = `dimer`). On `opt`, `grad` means L-BFGS (see {ref}`--opt-mode by command <opt-mode-semantics>`) |
| `--ref-mode` | path | `None` | Reference direction for the reaction mode (`.npz`, `.npy`, or text). `all` passes it from the MEP; you normally leave it unset. Not used by Dimer |
| `--microiter/--no-microiter` | flag | `True` | Alternate each macro TS step with an L-BFGS relaxation of the movable MM atoms. Not used by Dimer; with `--embedcharge`, the standard optimization runs instead |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | How the ML-region Hessian is computed |
| `--flatten/--no-flatten` | flag | `False` | Remove extra imaginary modes |
| `--freeze-atoms` | text | `None` | Extra atoms to freeze (1-based, comma-separated, for example `'1,3,5'`) |
| `--active-dof-mode` | `all` / `ml-only` / `partial` / `unfrozen` | `partial` | Atoms in the final frequency analysis: all atoms / ML region only / ML region and movable MM atoms / every atom outside the frozen layer |
| `--thresh` | preset | `baker` | Convergence criteria (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`) |
| `--max-cycles` | integer | `100000` | Maximum number of optimization cycles |
| `--stop-plateau/--no-stop-plateau` | flag | `False` | Stop when the energy stops changing (range below 1e-4 hartree over 50 cycles), then compute the Hessian |
| `--skip-final-freq/--no-skip-final-freq` | flag | `False` | Skip the final Hessian after convergence |
| `--read-hess` | path | `None` | Start from the Hessian in a `.npy` file instead of computing it (Cartesian, Hartree/bohr², for all atoms or only the movable ones) |
| `--dump-hess` | path | `None` | Save the final-geometry Hessian to a `.npy` file for `--read-hess` in `freq`, `tsopt`, or `irc`. Written only when the final Hessian is computed |
| `--precision` | `fp32` / `fp64` | per backend (`uma`: `fp32`; `orb`, `mace`: `fp64`) | Backend precision. `aimnet2` rejects `fp64` (see [MLIP Backends: Precision](backends.md#precision)) |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | Optimization coordinates: Cartesian / redundant internal / delocalized internal (DLC) / translation-rotation internal (TRIC) |
| `--config` | path | `None` | YAML file applied before the command-line options |
| `--dump/--no-dump` | flag | `False` | Write the optimization trajectory `optimization_all_trj.xyz` |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

For every option, run `mlmm tsopt --help-advanced` or see the [generated CLI reference](reference/commands/tsopt.md).

> **Note:** in YAML, Dimer reads the `hessian_dimer:` block, RS-P-RFO, RS-I-RFO, and TRIM share the `rsirfo:` block, and microiteration reads `microiter:`. Every key is listed under [`rsirfo`](yaml-reference.md#rsirfo), [`hessian_dimer`](yaml-reference.md#hessian_dimer), and [`microiter`](yaml-reference.md#microiter) in the YAML Reference.

> **Note:** if the reaction mode switches to another Hessian eigenvector (root) during the optimization, set `rsirfo.track_mode_by_overlap: true`.

> **Note:** if convergence is slow, lower `rsirfo.hessian_recalc` (default `500`) to 50–200 to recompute the exact Hessian more often, at the cost of more Hessian evaluations.

---

## Notes

(flatten-precedence-caveat)=
### When `--flatten` is on

`--flatten` is off by default. One YAML key, `hessian_dimer.flatten_max_iter`, sets the number of flatten rounds for every optimizer, Dimer and RS-P-RFO / RS-I-RFO / TRIM alike.

| Command line | Flatten rounds |
| --- | --- |
| Neither `--flatten` nor `--no-flatten` | `0` (off), unless YAML sets `hessian_dimer.flatten_max_iter` |
| `--flatten` | The YAML value if it is positive, otherwise `50` |
| `--no-flatten` | `0`, even if YAML sets a value |

`--flatten` removes extra imaginary modes. It cannot create a reaction mode that is missing; when n_imag = 0, get a better candidate instead.

### Other notes

* **Uphill steps are always allowed**: a saddle search has to go uphill along the reaction mode, so `tsopt` keeps `reject_uphill: false` even if YAML sets it. `--reject-uphill/--no-reject-uphill` belongs to `opt` and to the endpoint optimization in `all`.
* **Barrier from a product-side scan**: if the scan that made this candidate started from the product, read its barrier as described in [`scan` → Scan direction and barrier sign](scan.md#scan-direction-and-barrier-sign).
* **One root**: the optimizer climbs along one root (`0` is the lowest eigenvalue). Set it as a one-item list such as `rsirfo.roots: [0]`; Dimer uses `hessian_dimer.root`. `tsopt` has no `--root` flag.
* **Other RS-P-RFO settings**: `trust_norm: max_atom` limits the displacement of each atom instead of the whole step (Cartesian coordinates only), and `hessian_update: ts_bfgs` selects the TS-BFGS update instead of Bofill. Neither changes the trust radii.
* **Extra searches are opt-in**: after convergence, `tsopt` does not search further on its own, even when n_imag is not 1. Use `--flatten`, or set `rsirfo.saddle_recovery_max_cycles` above `0` (default `0`) to let RS-P-RFO / RS-I-RFO / TRIM step uphill when the exact Hessian shows no imaginary mode.
* **Flags that cannot be combined**: `--skip-final-freq` with `--dump-hess`; `--uma-workers` (parallel MLIP workers) above 1 with `--hessian-calc-mode Analytical` (see [Workers and Hessian mode](backends.md#workers-and-hessian-mode)).
* **`--skip-final-freq` and `--flatten`**: with RS-P-RFO / RS-I-RFO / TRIM, `--skip-final-freq` also skips `--flatten`, which needs the final Hessian.
* **`--read-hess` with RS-P-RFO / RS-I-RFO / TRIM**: the file replaces the first exact Hessian, so keep `rsirfo.hessian_init` at its default `calc`; other values stop with an error.
* **MM relaxation threshold in microiteration**: `microiter.micro_thresh` sets the convergence criteria of the MM relaxation; when it is unset, the MM relaxation uses the same criteria as the macro step.
* **Energy plateau and microiteration**: `--stop-plateau` watches the macro steps; it never stops the MM relaxation.
* **`--ref-mode` and frozen atoms**: `--ref-mode` only gives the reaction direction from the MEP; it does not change how frozen boundaries are treated.
* **At least one atom must move**: if every atom is frozen, `tsopt` stops with an error.
* **Errors**: invalid input or structure, or an error that keeps the optimizer from continuing, stops the run with an error message. Only the files written up to that point remain.
* **Option priority**: default < YAML < command line (see {ref}`Configuration precedence <configuration-precedence>`).

---

## See also

* [irc](irc.md) — follow the reaction path from the optimized TS
* [freq](freq.md) — full vibrational analysis and thermochemistry
* [path-opt](path-opt.md) / [path-search](path-search.md) / [scan](scan.md) — make a TS candidate
* [all](all.md) — model building, MEP, TS optimization, IRC, and frequencies in one run
* [Tips for studying reaction mechanisms](mechanism-tips.md) — what to try when the TS search fails
* [Troubleshooting](troubleshooting.md) — when a run fails
* [YAML Reference](yaml-reference.md) — every `rsirfo`, `hessian_dimer`, and `microiter` setting
* [Glossary](glossary.md) — TS, Dimer, Hessian, and other terms
* {ref}`Exit codes <exit-codes>` — what each exit status means
