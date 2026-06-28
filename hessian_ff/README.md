# `hessian_ff/` (bundled module — no upstream package)

> **This is a repo-internal module with NO upstream PyPI package. `pip install hessian-ff` returns 404 — bundling is mandatory.**

`hessian_ff` provides analytical Hessian computation on the MM force field — specifically, the AMBER ff14SB-style harmonic bonds, angles, propers, impropers, and Lennard-Jones / Coulomb terms — without round-tripping through a finite-difference loop. It is consumed exclusively by `mlmm/backends/mlmm_calc.py` for the MM-region Hessian during the 3-layer 5-pass partial Hessian assembly (chemistry-rule #8) and the link-atom B-matrix projection (chemistry-rule #2).


## Why bundled?

There is no upstream `hessian_ff` package on PyPI or any public registry. The analytical-Hessian implementation is research code originating from the lab and is shipped as a sibling module to `mlmm/` because:

1. **Single consumer**: only `mlmm/backends/mlmm_calc.py` calls into it; `analytical_hessian.py` is the sole entry point.
2. **Tight coupling to ONIOM math**: the analytical Hessian must match the link-atom B-matrix projection convention (chemistry-rule #2) and the 3-layer 5-pass assembly order (chemistry-rule #8) used in `mlmm_calc.py`.
3. **PyTorch-based (CPU-only) + parmed**: numerical correctness is fingerprinted against the FD-Hessian path during smoke tests.

## File map

| file | role |
|------|------|
| `analytical_hessian.py` | **The single entry point** — `build_analytical_hessian(system, coords, active_atoms=...)` returns `(H, info)`: the dense MM Hessian for the active subset plus an info dict |
| `forcefield.py` | force-field term definitions (bond / angle / proper / improper / LJ / Coulomb) |
| `prmtop_parmed.py` | parmed-based parm7 reader (the atom-indexing helpers of chemistry-rule #9 live in `loaders.py`) |
| `loaders.py` | force-field parameter loading |
| `system.py` | atom / topology data classes |
| `constants.py` | unit conversion constants |
| `terms/` | per-term analytical derivative code (one file per term type) |
| `native/` | **required** C-accelerated kernels for `build_analytical_hessian` (compiled at install; the automatic torch fallback was removed) |
| `workflows.py` | dead code — production never reaches this, but `__all__` declared = retention policy says do not delete during polish |
| `tests/` | unit tests for individual force-field terms |

## release scope

During this release, **only annotation edits are allowed** on this directory:

- docstring additions / improvements
- type hints
- section banners (`# ===... ===`)
- per-file module docstring

**Forbidden** during polish:

- any change to numerical behaviour, control flow, or function signatures of `analytical_hessian.py`
- any change to the per-term derivative code in `terms/`
- any change to the parmed atom-indexing helpers in `prmtop_parmed.py` (chemistry-rule #9 lives here)
- any new external dependency
- any deletion of `workflows.py` (declared `__all__` symbols mean future external use is possible)

Logic edits must be explicitly requested via a `[CHEMISTRY-RULE:2]` or `[CHEMISTRY-RULE:8]` or `[CHEMISTRY-RULE:9]` commit and go through HEAVY benchmark verification before merge.

## `mlmm_calc.py` is the sole entry point

```python
from hessian_ff.analytical_hessian import build_analytical_hessian

# Called inside the 3-layer 5-pass partial Hessian assembly:
H_mm, info = build_analytical_hessian(system, coords, active_atoms=movable_mm_mask)
```

If you are adding a new MM term, edit `forcefield.py` + add a new per-term file under `terms/` + register it in `analytical_hessian.py` — that single entry point routes everything. Do not add a parallel entry point. (Existing terms such as `terms/cmap.py` are already wired in.)

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.3, §6 — repo-internal fork policy + chemistry-rule #2, #8, #9 file:line
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.2 — do-not-touch list
- `THIRD_PARTY_NOTICES.txt` — third-party attributions (parmed is `pip install parmed`; only the analytical-Hessian glue is bundled)
