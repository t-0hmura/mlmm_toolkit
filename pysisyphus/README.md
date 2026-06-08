# `pysisyphus/` (bundled fork)

> **This is a repo-internal fork of [pysisyphus](https://github.com/eljost/pysisyphus), NOT the upstream PyPI package. Do not `pip install pysisyphus` alongside this package — it will silently overwrite the bundled copy and break long-running jobs.**

The fork is shipped inside this repository; treat it as part of `mlmm_toolkit` rather than a swappable upstream.

## Why a fork?

The upstream `pysisyphus` does not natively handle the constraints of full-protein ML/MM ONIOM jobs:

- **MLIP backends with autograd Hessians** evaluated on a GPU, while the ONIOM macro coordinates iterate on CPU
- **Macro / micro alternation** in the optimiser, where the ML region and the Movable MM region take alternating steps (chemistry-rule #3)
- **CPU-only bofill_update** for GPU OOM avoidance when the active-block partial Hessian is large
- **VRAM-aware stage handoff** — explicit `del` between IRC / tsopt / freq stages to free CUDA memory before the next stage loads its model
- **Initial-displacement memory hygiene** in IRC for full-protein systems (~10,000 atoms, 16 GB+ Hessians on the un-contracted ML-macro)
- **bofill_update advanced-indexing scatter** when only a subset of internal coordinates is updated (chemistry-rule #7)

The bundled fork patches these in five files and keeps everything else identical to upstream, so future upstream improvements remain easy to cherry-pick.

## Divergent files (do NOT replace with upstream)

| file | divergence | rule |
|------|------------|------|
| `pysisyphus/irc/IRC.py` | initial-displacement memory hygiene + contracted-ML-macro path + CPU stash of `forward_mw_hessian` to bound VRAM during backward integration; opt-in `require_pos_def_hessian` PSD convergence guard | freq-stage VRAM invariant, OOM bugfix |
| `pysisyphus/optimizers/HessianOptimizer.py` | opt-in `trust_band` rho-band trust update, `hessian_update_window` multistep TS-BFGS, `weighted_trust` per-coord-type L_inf trust; `get_xp`-based torch/numpy dispatch where shared API allows | TSopt step-control / trust radius |
| `pysisyphus/optimizers/hessian_updates.py` | `bofill_update` advanced-indexing scatter + CPU-only `bofill_update` for GPU OOM (CHEMISTRY-RULE:7); `multistep_ts_bfgs_update` helper; re-exports `_outer / _dot` from `_array` shim | scatter on subset of internals |
| `pysisyphus/optimizers/restrict_step.py` | `per_coord_type_weights` + `weighted_max_internal_step` helpers | weighted L∞ trust check |
| `pysisyphus/optimizers/gdiis.py` | `get_xp`-based torch/numpy dispatch (xp.linalg.norm / xp.sum) | torch/numpy backend share |
| `pysisyphus/tsoptimizers/TSHessianOptimizer.py` | RSIRFO kwargs and step-control for ONIOM macro/micro TS search | tsopt convergence + CHEMISTRY-RULE:3 |
| `pysisyphus/_array.py` | torch/numpy backend dispatch shim (`get_xp`, `_outer`, `_dot`, `_eigh`, `as_numpy`, `to_xp`) | used by `hessian_updates.py` + `HessianOptimizer.py` + `gdiis.py` |


## Scope vs upstream

This fork ships only the modules that mlmm_toolkit needs:
`Calculator` / `Dimer` (calculators), `cos.GrowingString`, `irc.EulerPC`,
`optimizers.{LBFGS, RFOptimizer, StringOptimizer, HessianOptimizer, gdiis,
restrict_step, hessian_updates,...}`, `tsoptimizers.{RSIRFOptimizer,
RSPRFOptimizer, TRIM, TSHessianOptimizer}`, `intcoords/`, `io/`, `cos/`,
`helpers`, `constants`, `elem_data`, `TablePrinter`, `xyzloader`,
`exceptions`, and the `_array` shim.

The upstream `pysisyphus run` CLI driver, the wavefunction-analysis subtree,
and the 30+ QM calculator backends are not part of this fork. Any external
caller that imported `from pysisyphus import run` must use mlmm's own
subcommand layer (`mlmm/cli/app.py`).

## release scope

During this release, **only annotation edits are allowed** on this directory:

- docstring additions / improvements
- type hints
- section banners (`# ===... ===`)
- per-file module docstring

**Forbidden** during polish:

- any change to numerical behaviour, control flow, or function signatures of the 4 divergent files above
- any new external dependency
- any rename of public symbols (would break `mlmm/tsopt.py`, `irc.py`, `path_opt.py`, `mlmm_calc.py` callers)

Logic edits to these files must be explicitly requested via a `[CHEMISTRY-RULE]` commit and that go through HEAVY benchmark verification before merge.

## Upstream compatibility

If you `pip install pysisyphus` into the same environment as `mlmm_toolkit`, Python's import machinery may resolve to the bundled copy or the upstream copy depending on `sys.path` order. The flat-top placement at `<repo-top>/pysisyphus/` is required for the bundled copy to take precedence in the editable-install path; do not move it under `mlmm/`.

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.3, §6 — repo-internal fork policy
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.2 — do-not-touch list (5 divergent files)
- `THIRD_PARTY_NOTICES.txt` — pysisyphus upstream attribution
