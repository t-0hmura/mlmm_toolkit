# Orb backend (orb.md)

The Orb backend provides an energy-conservative MLIP option through
`orb-models`. Validate energies, forces, optimized structures, and
frequencies on the target system before selecting it for a workflow.

## Install

```bash
pip install 'mlmm-toolkit[orb]'         # pulls orb-models
```

The current ORB extra installs `orb-models`. If installation fails, inspect the
actual resolver error and `python -m pip check` instead of adding unrelated PyG
packages.

Or, if `mlmm-toolkit` is already installed:

```bash
pip install orb-models
```

Confirm:

```bash
python -c "import orb_models; print('orb backend OK:', orb_models.__version__)"
```

Orb model weights are downloaded on first use; no separate auth required.

## CLI usage

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b orb
```

Default model: `orb_v3_conservative_omol`. Check the installed
`orb-models` model card for checkpoint provenance, supported elements, and
runtime requirements; dataset coverage alone does not establish checkpoint
capability.

Inspect the default kwarg dict:

```bash
python -c "import mlmm.core.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

Orb accepts these `MLMM_CALC_KW` keys (the Orb-specific `orb_model` /
`orb_precision` are `_OrbBackend.__init__` parameters in `backends/mlmm_calc.py`;
the Hessian/calc keys below apply to every backend; defaults in
`core/defaults.py`):

| Key | Purpose |
|---|---|
| `model_charge`, `model_mult` | Total charge and spin multiplicity |
| `ml_device` | `'cuda'`, `'cpu'`, `'auto'` |
| `orb_model` | Override the default Orb checkpoint |
| `orb_precision` | `'float64'` (default) or `'float32-high'` for reduced-precision screening (key in `MLMM_CALC_KW`; `'float32'` normalizes to `'float32-high'`) |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `H_double` | Same as UMA |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| Conservative energy/force model with a reduced-precision option | Backend-specific TS and frequency behavior must be validated for the target system |
| Backend-specific precision selection | Backend-specific TS and frequency behavior must be validated for the target system |
| Easy installation through the extra | Check checkpoint element and state coverage before use |

Compare candidate geometries and frequencies against the backend selected for
production before mixing backends across a workflow.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| Extra imaginary modes | Inspect the modes and compare supported precisions on the target system before independently recomputing frequencies/IRC. |
| TS has more than one imaginary mode | The result is not a certified first-order saddle; tighten/restart and inspect all mode displacements. |

## See also

- `env-cuda.md` — torch / CUDA prerequisites.
- `uma.md` — UMA backend setup.
- `mace.md` — MACE backend setup (separate env).
- `mlmm-cli/tsopt.md` — diagnosing TS convergence problems.
