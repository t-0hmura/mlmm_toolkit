# MLIP Backends

mlmm-toolkit drives every ML/MM workflow stage (`opt`, `scan`, `tsopt`, `freq`,
`irc`, `path-search`,...) through a single `MLMMCore` ONIOM-coupling object.
`MLMMCore` dispatches the ML region to a per-backend adapter (`_UMABackend` /
`_OrbBackend` / `_MACEBackend` / `_AIMNet2Backend`) via the private
`_create_ml_backend` factory. This page documents how to select a backend, the
per-backend kwargs, and how to add a new backend.

## Public surface

```python
from mlmm.backends.mlmm_calc import MLMMCore, MLMMASECalculator, mlmm

# Primary entry: ML/MM ONIOM core. Takes parm7 + layered PDB + model-PDB and
# returns energy / forces / (partial) Hessian on the full system.
core = MLMMCore(
    input_pdb="layered.pdb",
    real_parm7="real.parm7",
    model_pdb="model.pdb",
    backend="uma",          # one of: "uma", "orb", "mace", "aimnet2"
    model_charge=0, model_mult=1,
    uma_model="uma-s-1p1",
    uma_precision="fp32",   # or "fp64" (full-precision base inference)
)

# ASE adapter (DMF and other ASE-based stages)
ase_calc = MLMMASECalculator(core)

# pysisyphus Calculator adapter (opt / tsopt / freq / irc / path-search stages)
pysis_calc = mlmm(
    input_pdb="layered.pdb",
    real_parm7="real.parm7",
    model_pdb="model.pdb",
    backend="uma",
    model_charge=0, model_mult=1,
)
```

Internally `MLMMCore.__init__` calls `_create_ml_backend(backend, ...)` (a
private factory in `mlmm/backends/mlmm_calc.py`) to instantiate the right
adapter. The factory raises `ValueError` for unknown backends; there is no
`'auto'` fallback in mlmm — workflow code passes the resolved backend
name from the CLI.

## File map

| file | role |
|------|------|
| `mlmm/backends/__init__.py` | `apply_precision_to_calc_cfg()` — routes the unified `--precision fp32\|fp64` CLI flag to each backend's native kwarg (`uma_precision` / `orb_precision` / `mace_dtype`) |
| `mlmm/backends/mlmm_calc.py` | `MLMMCore` (ML/MM ONIOM coupling) + `MLMMASECalculator` (ASE) + `mlmm` (pysisyphus Calculator) + per-backend adapters (`_UMABackend`, `_OrbBackend`, `_MACEBackend`, `_AIMNet2Backend`) + the private `_create_ml_backend` factory + FD-Hessian assembly + unit conversion |
| `mlmm/backends/xtb_embedcharge_correction.py` | xTB point-charge embedding correction (`--embedcharge` flag) for MM→ML environmental effects |

## Per-backend characteristics

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install fairchem-core` + HF auth | `uma-s-1p1` / `uma-m-1p1` | `uma_precision="fp32" \| "fp64"` |
| `orb` | `pip install orb-models` | `orb_v3_conservative_omol` | `orb_precision="float32-high" \| "float64"` (`"float32"` accepted as an alias) |
| `mace` | dedicated env: `pip uninstall -y fairchem-core && pip install mace-torch` (`e3nn` pin conflicts with UMA) | `MACE-OMOL-0` | `mace_dtype="float32" \| "float64"` |
| `aimnet2` | `pip install aimnet` | `aimnet2` | n/a |

### UMA fp64

Switching OMol-trained UMA from the default fp32 to fp64 can have a non-trivial
impact on TSopt and Hessian evaluation. Enable via:

```bash
mlmm tsopt -i ts.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
mlmm freq -i opt.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
mlmm irc -i ts.pdb --parm real.parm7 -q 0 -m 1 --precision fp64...
```

The unified `--precision` flag is routed to each backend's native kwarg
(`uma_precision` for UMA, `orb_precision` for ORB, `mace_dtype` for MACE)
by `apply_precision_to_calc_cfg` in `mlmm/backends/__init__.py`.

The unified `--backend-model NAME` flag likewise overrides the model variant
for the selected `--backend`, routed to the backend's model kwarg
(`uma_model` / `orb_model` / `mace_model` / `aimnet2_model`) by
`apply_backend_model_to_calc_cfg`. Unset keeps the backend's default model.

Or via YAML config (per-backend kwarg name):

```yaml
calc:
 uma_precision: fp64
```

Requires `fairchem-core ≥ 2.0` for the `InferenceSettings` API.

## Custom backend — bring your own ASE Calculator (`--calc-file`)

Beyond the built-in MLIP backends, the **ML region** can be driven by any
[ASE](https://wiki.fysik.dtu.dk/ase/) Calculator supplied at run time with
`--calc-file`, without modifying `mlmm_toolkit`. This couples the ML side of the
ML/MM ONIOM scheme to GFN-xTB (via `tblite` / `xtb-python`), DFTB+, ORCA, Psi4,
or any ASE-compatible engine — the boundary is the standard ASE Calculator
interface (energy in eV, forces in eV/Å), wrapped by the existing
`_ASEMLBackend` adapter.

Write a Python file exposing a `get_calculator` factory that returns an ASE
Calculator:

```python
# my_calc.py  (minimal illustrative example)
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

Swap `EMT()` for the engine you want — e.g. `tblite.ase.TBLite(...)` for
GFN-xTB, the DFTB+ ASE calculator, or `ase.calculators.orca.ORCA(...)`. Then
pass the file to any single-stage subcommand (it selects the `custom` ML
backend, overriding `--backend`):

    mlmm sp    -i complex.pdb --parm system.parm7 --calc-file my_calc.py -q 0 -m 1
    mlmm opt   -i complex.pdb --parm system.parm7 --calc-file my_calc.py
    mlmm freq  -i complex.pdb --parm system.parm7 --calc-file my_calc.py

Notes:

- The factory receives `charge`, `spin` (multiplicity; also offered as `mult` /
  `multiplicity`), and `device` when its signature accepts them, or
  unconditionally if it declares `**kwargs`, so engines that need the total
  charge (e.g. xTB) can be configured. Use a different factory name with
  `--calc-factory NAME`; a module-level Calculator instance is also accepted.
- The custom calculator drives the **ML region only**; the MM side keeps its
  usual `hessian_ff` / OpenMM backend and the ONIOM coupling is unchanged.
  Hessians use the finite-difference path, so `freq` and `tsopt --opt-mode hess`
  work with any engine. Frozen atoms are honored as usual.
- Available on the standalone subcommands (`sp`, `opt`, `tsopt`, `freq`, `irc`,
  `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`). For a permanent,
  installable backend with its own `--backend` name instead, see the recipe
  below.

## Add-a-backend recipe (5 steps)

To add a new backend `XYZModel` exposed as `--backend xyz`:

1. **Create backend adapter** in `mlmm/backends/mlmm_calc.py` (or a new file
 like `mlmm/backends/xyz.py` if it grows large): implement `_XYZBackend(_MLBackend)`
 inheriting from `_MLBackend` (the ABC) and a parallel `_XYZASEBackend` if an
 ASE path is needed. Both must accept the common kwargs `charge / spin / device /
 freeze_atoms / hessian_calc_mode / return_partial_hessian / hessian_double /
 print_timing / model` and any backend-specific kwargs (`precision`,
 `default_dtype`, etc.).
2. **Conform to `_MLBackend`**: implement the abstract methods
 `eval(atoms, need_grad=True) -> (E_eV, F_eV, opaque)` (energy in eV, forces in
 eV/Å, plus a backend-specific opaque object), `hessian_analytical(opaque, n_atoms,
 *, dtype) -> torch.Tensor` (returns the Hessian in eV/Å²), and the
 `supports_analytical_hessian` and `device` properties. Inherit from `_MLBackend`
 to pick up the generic finite-difference `hessian_fd(...)` for free (used when the
 backend has no analytical Hessian).
3. **Register in `_create_ml_backend`**: extend the factory in
 `mlmm/backends/mlmm_calc.py` to dispatch `backend == "xyz"` to `_XYZBackend(...)`.
 Add the new backend's kwargs to the `_create_ml_backend(...)` signature so
 `MLMMCore` can forward them.
4. **Wire the unified `--precision` flag** (optional): if the backend exposes
 a precision knob, add a `"xyz": (kw_name, kw_value)` entry under both `"fp32"`
 and `"fp64"` in `_PRECISION_DISPATCH` in `mlmm/backends/__init__.py` so the
 user-facing `--precision fp32|fp64` CLI flag routes correctly.
5. **Document + smoke**: add an entry to this page's file map / per-backend table,
 document model identifiers + install command, and add an `xyz` line to
 `tests/smoke/run.sh` so the new backend is exercised end-to-end.

## VRAM invariant (ML/MM-specific)

In an ML/MM ONIOM job the ML backend co-resides with PySCF (DFT correction),
parmed (parm7), and the MM force-field arrays on the same device. The per-direction
FD-Hessian loop in `mlmm/backends/mlmm_calc.py` is sized to fit this combined
footprint. **Do not refactor the per-direction loop into a batched tensor**
without re-running the full GPU smoke gate (`tests/smoke/run.sh`) and watching
peak VRAM — full-protein ML/MM `all` jobs OOM otherwise. The stage runners then `del calc`
between stages and the `all` workflow runs `gc.collect()` at stage boundaries;
this is part of the public contract and must not be removed when refactoring the
workflow.

## ONIOM coupling vs raw MLIP

The MLIP adapters in `mlmm/backends/mlmm_calc.py` evaluate the **ML region only**.
The subtractive ONIOM energy formula (`# CHEMISTRY-RULE:1`), the link-atom Hessian
B-matrix projection (`# CHEMISTRY-RULE:2`), and the 3-layer 5-pass partial Hessian
assembly (`# CHEMISTRY-RULE:8`) live in the same file. Backend authors adding a
new MLIP do not need to know the ONIOM coupling; they only need to expose a
calculator that returns ML-region energy / force / Hessian in the correct units.

## See Also

- [Python API](python-api.md) — `MLMMCore` / `MLMMASECalculator` / `mlmm` (pysisyphus Calculator) public surface.
- [Architecture](architecture.md) — 6-layer directory map + dependency direction.
- [CONTRIBUTING](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) — Recipe 3.2 "Add an MLIP backend" with full gate cycle references.
