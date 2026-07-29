# `hessian_ff/` (bundled module — no upstream package)

> **This is a repo-internal module with NO upstream PyPI package. `pip install hessian-ff` returns 404 — bundling is mandatory.**

`hessian_ff` provides analytical Hessian computation on the MM force field — specifically, the AMBER ff14SB-style harmonic bonds, angles, propers, impropers, and Lennard-Jones / Coulomb terms — without round-tripping through a finite-difference loop. It is consumed exclusively by `mlmm/backends/mlmm_calc.py` for the MM-region Hessian during the 3-layer 5-pass partial Hessian assembly (chemistry-rule #8) and the link-atom B-matrix projection (chemistry-rule #2).


## Why bundled?

There is no upstream `hessian_ff` package on PyPI or any public registry. The module is primarily lab-originated research code; its CMAP interpolation/coefficient material in `terms/cmap.py` and `native/bonded_ext.cpp` is adapted from OpenMM under the MIT license documented in `THIRD_PARTY_NOTICES.txt`. It is shipped as a sibling module to `mlmm/` because:

1. **Single consumer**: only `mlmm/backends/mlmm_calc.py` calls into it; `analytical_hessian.py` is the sole entry point.
2. **Tight coupling to ONIOM math**: the analytical Hessian must match the link-atom B-matrix projection convention (chemistry-rule #2) and the 3-layer 5-pass assembly order (chemistry-rule #8) used in `mlmm_calc.py`.
3. **PyTorch-based (CPU-only) + parmed**: numerical correctness is fingerprinted against the FD-Hessian path during smoke tests.

## File map

| file | role |
|------|------|
| `analytical_hessian.py` | **The single entry point** — `build_analytical_hessian(system, coords, active_atoms=...)` returns `(H, info)`: the dense MM Hessian for the active subset plus an info dict |
| `forcefield.py` | force-field term definitions (bond / angle / proper / improper / LJ / Coulomb) |
| `prmtop_parmed.py` | parmed-based parm7 reader; chemistry-rule #9 atom-index normalization lives in `mlmm/io/pdb_indexing.py` |
| `loaders.py` | force-field parameter loading |
| `system.py` | atom / topology data classes |
| `constants.py` | unit conversion constants |
| `terms/` | per-term analytical derivative code (one file per term type); `terms/cmap.py` includes OpenMM-derived MIT-licensed material |
| `native/` | **required** C-accelerated kernels for `build_analytical_hessian` (JIT-compiled on first use at runtime via `torch.utils.cpp_extension`; needs GCC ≥ 9 + ninja); `native/bonded_ext.cpp` includes OpenMM-derived MIT-licensed material |
| `workflows.py` | compatibility API declared through `__all__`; not imported by the production path |
| `tests/` | unit tests for individual force-field terms |

## Change policy

Numerical or control-flow changes require a demonstrated defect or feature
need, focused derivative/parity tests, and the relevant CPU/GPU validation.
Preserve the `analytical_hessian.py` public contract and the atom-index mapping
used by `prmtop_parmed.py`. `workflows.py` exposes compatibility symbols through
`__all__` even though the production path does not import it.

Validate logic changes with the relevant unit tests and scheduled numerical benchmark.

## `mlmm_calc.py` is the sole entry point

```python
from hessian_ff.analytical_hessian import build_analytical_hessian

# Called inside the 3-layer 5-pass partial Hessian assembly:
H_mm, info = build_analytical_hessian(system, coords, active_atoms=movable_mm_mask)
```

If you are adding a new MM term, edit `forcefield.py` + add a new per-term file under `terms/` + register it in `analytical_hessian.py` — that single entry point routes everything. Do not add a parallel entry point. (Existing terms such as `terms/cmap.py` are already wired in.)

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.1, §6 — chemistry-rule locations and repo-internal fork policy
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.3 — bundled-fork edit policy
- [`../THIRD_PARTY_NOTICES.txt`](../THIRD_PARTY_NOTICES.txt) — third-party attributions, including the exact OpenMM-derived files and license
