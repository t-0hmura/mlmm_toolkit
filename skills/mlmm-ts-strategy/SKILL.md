---
name: mlmm-ts-strategy
description: >-
  Decision guidance for ML/MM enzyme reaction-barrier campaigns: backend-specific
  precision, TS-candidate routes, exact first-order-saddle validation, IRC early-stop
  handling, scan direction/staging, and controlled mutant comparisons with
  chemically corresponding regions. Use for barrier, imaginary-frequency,
  wrong-saddle, precision, MEP-vs-restraint, IRC connectivity, scan-direction, or
  mutant-comparison questions. Skip installation, pure structure-format editing,
  and MCP transport.
---

# mlmm ts-strategy

Cross-cutting decisions for getting a *correct* reaction barrier out of an ML/MM
ONIOM campaign. Every flag below is verified against `mlmm/cli/common_options.py`,
`mlmm/workflows/{scan,opt,path_search,all}.py`, and `mlmm/core/defaults.py`.

## 1. Precision: preserve backend defaults

| Backend | Unset default | Guidance |
|---|---|---|
| UMA | fp32 | Compare supported precisions on the target system when precision matters. |
| ORB | fp64 | Compare supported precisions on the target system when precision matters. |
| MACE | fp64 | Compare supported precisions on the target system when precision matters. |
| AIMNet2 | fp32 | No precision switch; explicit fp64 is rejected. |

- `--precision` is case-insensitive and unset resolves through the selected backend.
- Backend routing: `uma`→`uma_precision`; `orb`→`orb_precision` (`float32-high`|`float64`); `mace`→`mace_dtype`; `aimnet2`→fp32 is a no-op and **fp64 is rejected** (inputs cast to float32 upstream).
- Accepted on `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`/`scan2d`/`scan3d`, `path-opt`, `path-search`, `all`.
- fp64 changes numerical precision. `--deterministic` requests deterministic
  algorithms, but exact repeatability still requires verification of the
  installed backend/model/SDK and target stack (see `reproducibility.md`).

## 2. Two routes to a TS candidate

| Route | Subcommand | Mechanism | Use when |
|---|---|---|---|
| (a) MEP / path-search | `path-search` (or `path-opt` for one segment) | Recursive GSM/DMF segmentation; returns HEI and segment candidates for TS/IRC validation | You have R and P, optionally with ordered intermediates |
| (b) Distance-restrained build-up | `scan` (`scan2d`/`scan3d`) | Harmonic restraint `E = ½·k·(r_ij − target)²` (scan default `k=300` via `BIAS_KW`; the `10.0` in `restraints.py` `HarmonicBiasCalculator` is only an unused constructor fallback) drives the reacting distance(s) toward the barrier with L-BFGS relaxation | No usable second endpoint / TS guess — drive the reacting bond directly |

- There is **no `opt --restraint` flag**, but `opt` supports restrained optimization via `--dist-freeze` (with `--bias-k`, default k=300, the same `HarmonicBiasCalculator`); `scan` additionally drives staged target distances up to a TS candidate.
- `path-search` (`app.py`: "Search reaction pathways recursively.") auto-segments a multistep path; `path-opt` optimizes a single given segment.
- Feed a TS candidate from either route into `tsopt → irc` (or `all --tsopt`). Terminal PHVA checks saddle order; add `freq` for full modes or thermochemistry.

## 3. Wrong imaginary-frequency count at TS-opt

A certified first-order saddle has **exactly one** imaginary mode; inspect its
displacement and IRC connectivity. Two or more fail certification regardless
of magnitude.

| Symptom | Action |
|---|---|
| Extra imaginary modes | Inspect all mode displacements, the MEP guess, optimizer stop reason, and backend-specific numerical behavior; then retry an appropriate coordinate/flattening/precision setting and check the new terminal PHVA result. |
| Collapsed to `n_imag = 0` | Treat as failed, not as a TS. Improve the MEP/initial guess; `--flatten` only removes surplus modes and cannot create a missing reaction direction. |
| Poor MEP/HEI | In `all`, try opt-in `--refine-path` before TS optimization. It can split a poor path into several stages and increase cost, so it is off by default. |
| Still no clean saddle | Revise endpoints/scan coordinates and verify that the single imaginary mode moves the reacting atoms. |

- `--flatten`/`--no-flatten` runs the extra-imaginary-mode flattening loop. It is available on `opt`, `tsopt`, and `all`; the TS commands apply their TS-specific loop when `n_imag > 1`. `--flatten` uses the configured positive iteration cap and `--no-flatten` forces zero.
- `opt` and `tsopt` accept `--coord-type cart|redund|dlc|tric`; `all` accepts `cart|dlc`. The effective default is `cart`.
- `dlc` = delocalized internal coordinates. Its cost and convergence behavior
  are system-dependent; compare against `cart` on the same seed.
- `dlc` requires a **Hessian-based optimizer**: in `opt.py`, `--coord-type dlc` with L-BFGS (`--opt-mode grad`) is forced back to `cart` with a warning. Use it on `tsopt` (RS-P-RFO / RS-I-RFO / TRIM) or `opt --opt-mode hess`.
- `path-opt`/`path-search` have no `--coord-type` flag; they take the coordinate system from `--config` YAML (`geom.coord_type`), and pysisyphus ChainOfStates supports only `cart`/`dlc` there.
- `cart` is the default. Independently validate a change of coordinate system
  with frequency analysis and IRC connectivity.

`--ref-mode` is an advanced path-direction input, not a routine standalone
remedy. It accepts one or more atom-order-matched Cartesian 3N candidates from
`.npz`, `.npy`, or whitespace text and guides negative Hessian-root identity and
overlap; it does not replace the Hessian. `mlmm all` supplies CPU/file-cached
MEP candidates to Hessian TS optimizers by default. Dimer does not consume it.
With `all --no-tsopt-from-mep-tan`, cache creation/use is disabled and TSOPT
selects its initial root from the initial-structure Hessian modes.

Keep numerical convergence separate from saddle order. A converged
higher-order stationary point is not a first-order TS, although `all` may run
warning-labelled diagnostic IRC when a validated negative root exists.
Numerical non-convergence, no imaginary mode, failed/skipped PHVA, or no valid
negative root stops after preserving the TS artifacts and before IRC.

## 4. IRC stops too early

First reduce the step length, for example `mlmm irc ... --step-size 0.05`.
Use `--never-stop` (or `all --irc-never-stop`) when the intended operation is
unconditional tracing to the cycle cap. It ignores gradient and energy endpoint
criteria; numerical/integration failures still stop the run. Always inspect
both branches and bond connectivity. The mode is off by default.

## 5. Reading the barrier when the scan started from Product

If the scan/path **starts from P**, the raw reported barrier is the **reverse** direction.

| Quantity | Formula |
|---|---|
| Forward barrier | `E(TS) − E(reactant)` |
| Reverse barrier (raw P-start number) | `E(TS) − E(product)` |

- This is a *read-time interpretation*, **not a CLI flag**. Always confirm which endpoint is R vs P (read `segments/seg_NN/{reactant,product}.pdb` from the IRC, not the scan direction).

## 6. Staged vs concerted scan

`-s/--scan-lists` is `multiple=True` (`scan.py`). Help: "Multiple inline literals define sequential stages."

| Form | Invocation | Meaning | Needs mechanism up front? |
|---|---|---|---|
| Concerted | **single** `--scan-lists` literal with several `(i,j,target)` tuples | all coords driven together in one stage | No |
| Staged | one `--scan-lists` followed by several literals | each literal is one sequential restrained relaxation, written to `stage_NN/` | Yes — define the mechanism per stage |

```bash
# Concerted (one stage, two coords driven together):
mlmm scan -i r.pdb --parm e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# Staged (two sequential stages):
mlmm scan -i r.pdb --parm e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' '[(7,9,0.95)]' -o result_staged
```

- `path-search` does multistep auto-segmentation, so a **concerted** scan needs no mechanism breakdown.
- A **staged** scan needs the mechanism defined up front and exposes each prescribed coordinate change as a separate stage.
- A 4-tuple expands into 2 stages (bidirectional scan).

## 7. Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

Form each barrier within one chemical system (`TS - R` or `TS - P`) before
comparing mutant and wild type. Absolute energies of systems with different
compositions are not directly subtracted.

- Keep the backend/method, force field, convergence criteria, thermochemistry
  settings, and temperature matched.
- Use chemically corresponding ML and movable regions. A mutation may change
  the atom count; assign every new/deleted atom rather than requiring
  byte-identical layers.
- Determine charge and multiplicity independently for each model.
- Certify each TS independently with exactly one imaginary mode, inspect its
  displacement, and verify both IRC endpoint identities.

## See also

- `mlmm-cli/tsopt.md`, `scan.md`, `path-search.md`, `define-layer.md` — per-subcommand flags.
- `mlmm-cli/all-ts-only.md` — the full mutate→complete→transplant→run mutation recipe.
- `mlmm-workflows-output/SKILL.md` — IRC R/TS/P canonical paths and bond-change conventions.
- `mlmm-hpc/SKILL.md` — choosing the GPU class that determines §1.
- docs `reproducibility.md` — fp64 vs `--deterministic`.
