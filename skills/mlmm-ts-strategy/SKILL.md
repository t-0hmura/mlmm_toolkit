---
name: mlmm-ts-strategy
description: >-
  Decision guidance for ML/MM enzyme reaction-barrier campaigns: backend-specific
  precision, TS-candidate routes, exact first-order-saddle recovery, IRC early-stop
  handling, scan direction/staging, and controlled mutant comparisons with an
  identical atom/layer selection. Use for barrier, imaginary-frequency,
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
| UMA | fp32 | Try explicit fp64 for a numerically sensitive TS/Hessian; validate cost on the target GPU. |
| ORB | fp64 | Keep fp64 for TS/frequency work. Explicit fp32 is `float32-high`/TF32 and is screening-only. |
| MACE | fp64 | Keep the upstream float64 default for TS/frequency work. |
| AIMNet2 | fp32 | No precision switch; explicit fp64 is rejected. |

- `--precision` is case-insensitive and unset resolves through the selected backend.
- Backend routing: `uma`→`uma_precision`; `orb`→`orb_precision` (`float32-high`|`float64`); `mace`→`mace_dtype`; `aimnet2`→fp32 is a no-op and **fp64 is rejected** (inputs cast to float32 upstream).
- Accepted on `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`/`scan2d`/`scan3d`, `path-opt`, `path-search`, `all`.
- fp64 *reduces* GPU reduction-order drift but does **not** make a run bit-identical — only `--deterministic` does (see `reproducibility.md`).

## 2. Two routes to a TS candidate

| Route | Subcommand | Mechanism | Use when |
|---|---|---|---|
| (a) MEP / path-search | `path-search` (or `path-opt` for one segment) | Recursive GSM/DMF segmentation; brackets the TS between endpoints, auto-bridges gaps, one TS per segment | You have R (and optionally P / intermediates) and want the path discovered |
| (b) Distance-restrained build-up | `scan` (`scan2d`/`scan3d`) | Harmonic restraint `E = ½·k·(r_ij − target)²` (scan default `k=300` via `BIAS_KW`; the `10.0` in `restraints.py` `HarmonicBiasCalculator` is only an unused constructor fallback) drives the reacting distance(s) toward the barrier with L-BFGS relaxation | No usable second endpoint / TS guess — drive the reacting bond directly |

- There is **no `opt --restraint` flag**, but `opt` supports restrained optimization via `--dist-freeze` (with `--bias-k`, default k=300, the same `HarmonicBiasCalculator`); `scan` additionally drives staged target distances up to a TS candidate.
- `path-search` (`app.py`: "Search reaction pathways recursively.") auto-segments a multistep path; `path-opt` optimizes a single given segment.
- Feed a TS candidate from either route into `tsopt → irc → freq` (or `all --tsopt`).

## 3. Wrong imaginary-frequency count at TS-opt

A clean first-order saddle = **exactly one** dominant imaginary mode along the reaction coordinate.

| Symptom | Action |
|---|---|
| Extra imaginary modes | Keep ORB/MACE at fp64; consider UMA fp64, `--coord-type dlc`, and `--flatten`. Independently rerun `freq` and inspect the reaction mode. |
| Collapsed to `n_imag = 0` | Treat as failed, not as a TS. Improve the MEP/initial guess; `--flatten` only removes surplus modes and cannot create a missing reaction direction. |
| Poor MEP/HEI | In `all`, try opt-in `--refine-path` before TS optimization. It can split a poor path into several stages and increase cost, so it is off by default. |
| Still no clean saddle | Revise endpoints/scan coordinates and verify that the single imaginary mode moves the reacting atoms. |

- `--flatten`/`--no-flatten` (`tsopt.py` cli `default=None`): runs the extra-imaginary-mode flattening loop (`grad`: dimer loop; `hess`: post-RS-I-RFO step, triggered when `n_imag > 1`). `--flatten` uses `flatten_max_iter=50`; `--no-flatten` forces `flatten_max_iter=0`. Available on `tsopt`/`all`. Try `fp64` and/or `dlc` first, then add `--flatten` for a residual small mode — they are complementary.
- `--coord-type` = `click.Choice(['cart','redund','dlc','tric'])`, case-insensitive; `'dlc'` is present verbatim (`common_options.py` `add_coord_type_option`). Effective default `'cart'` (`defaults.py` `GEOM_KW_DEFAULT['coord_type']='cart'`).
- `dlc` = delocalized internal coordinates. Its cost and convergence behavior
  are system-dependent; compare against `cart` on the same seed.
- `dlc` requires a **Hessian-based optimizer**: in `opt.py`, `--coord-type dlc` with L-BFGS (`--opt-mode grad`) is forced back to `cart` with a warning. Use it on `tsopt` (RFO/RS-I-RFO) or `opt --opt-mode hess`.
- `path-opt`/`path-search` have no `--coord-type` flag; they take the coordinate system from `--config` YAML (`geom.coord_type`), and pysisyphus ChainOfStates supports only `cart`/`dlc` there.
- `cart` is the default. Independently validate a change of coordinate system
  with frequency analysis and IRC connectivity.

`--ref-mode` is an advanced path-direction input, not a routine standalone
remedy. `mlmm all` derives and supplies the normalized 3N Cartesian tangent
from its MEP so saddle recovery can reject a nearby minimum. Ordinary
standalone `tsopt` should omit it unless an externally derived, atom-order-
matched reaction mode is available.

## 4. IRC stops too early

First reduce the step length, for example `mlmm irc ... --step-size 0.05`.
This is the preferred response to an early plateau/energy-rise stop. If a small
shoulder must be crossed, opt in to `--never-stop` (or
`all --irc-never-stop`). It ignores energy-rise and plateau stop conditions;
integrator convergence, invalid values, and `--max-cycles` still stop the run.
Always inspect both branches and bond connectivity. The mode is off by default.

## 5. Reading the barrier when the scan started from Product

If the scan/path **starts from P**, the raw reported barrier is the **reverse** direction.

| Quantity | Formula |
|---|---|
| Forward barrier | `E(TS) − E(reactant)` |
| Reverse barrier (raw P-start number) | `E(TS) − E(product)` |

- This is a *read-time interpretation*, **not a CLI flag**. Always confirm which endpoint is R vs P (read `segments/seg_NN/{reactant,product}.pdb` from the IRC, not the scan direction).
- Concrete case (CM): the scan begins at P, so the campaign's "barrier" is reverse; the forward barrier = `E(TS) − E(R)`.

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
- A **staged** scan needs the mechanism defined up front, **but when the mechanism is known, staged gives cleaner per-step control and is generally preferred.**
- A 4-tuple expands into 2 stages (bidirectional scan).

## 7. Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

**Rule: every compared model MUST use the SAME atom set — identical atom count and residues — or it is not a controlled experiment.** A different atom set changes the energy reference and invalidates the ΔE‡ comparison; a geometrically re-derived ML/movable/frozen partition on the mutant also produces spurious soft modes (`tsopt.n_imaginary ≥ 2`, both tiny → IRC aborts).

mlmm recipe (transplant the WT ML/MM layer encoding onto the mutant; run with `--detect-layer`):

| step | action |
|---|---|
| 1 | Build the mutant structure; keep the same residue set as WT (only the mutated residue's identity differs). |
| 2 | Transplant WT's B-factor layer encoding onto the mutant by `(resid, atom-name)`: **ML=0.0, MovableMM=10.0, FrozenMM=20.0** (`common_options.py` `add_ml_layer_detection_options`). |
| 3 | Run with `--detect-layer` (default on) so the SAME layer assignment is reused. |

```bash
mlmm all -i mutant_layered.pdb -l 'LIG:Q' \
    --tsopt --thermo [--dft --dft-func-basis '...'] \
    -o result_mutant
```

| flag | do | why |
|---|---|---|
| `--detect-layer` | KEEP (default `True`) | reads transplanted B-factors → ML/movable/frozen byte-identical to WT |
| `-c/--center`, `-r/--radius` | **OMIT** | geometric extraction would re-derive a *different* pocket on the mutant (`all.py`: omitting `-c` skips extraction, uses the full structure as-is) |
| `-l/--ligand-charge` | KEEP | a non-standard ligand's charge isn't in the standard-AA table; `-l 'RES:Q'` auto-derives the total. Prefer over hardcoded `-q`. |
| `--movable-cutoff` | do NOT pass | it **disables** `--detect-layer` (`scan.py`) |

- Same principle applies to comparing two mechanisms on the same enzyme: identical atom set across both, only the reaction coordinate differs.

## See also

- `mlmm-cli/tsopt.md`, `scan.md`, `path-search.md`, `define-layer.md` — per-subcommand flags.
- `mlmm-cli/all-ts-only.md` — the full mutate→complete→transplant→run mutation recipe.
- `mlmm-workflows-output/SKILL.md` — IRC R/TS/P canonical paths and bond-change conventions.
- `mlmm-hpc/SKILL.md` — choosing the GPU class that determines §1.
- docs `reproducibility.md` — fp64 vs `--deterministic`.
