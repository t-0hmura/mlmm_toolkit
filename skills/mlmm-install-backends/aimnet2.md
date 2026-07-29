# AIMNet2 backend (aimnet2.md)

AIMNet2 is available through the `aimnet` package. Check the installed
checkpoint/model card for its element, charge, multiplicity, and system-size
domain, then validate energies, forces, structures, and frequencies on the
target system.

## Element coverage

Element support is checkpoint-specific. Do not infer support from the package
name or another model generation; inspect the installed checkpoint contract.

## Install

```bash
pip install 'mlmm-toolkit[aimnet]'         # pulls aimnet>=0.2.0
```

Or, if `mlmm-toolkit` is already installed:

```bash
pip install 'aimnet>=0.2.0'
```

Confirm:

```bash
python -c "import aimnet; print('aimnet:', aimnet.__version__)" && echo "mlmm + aimnet2 backend OK"
```

## CLI usage

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'GLU' -l 'GLU:-1' \
    --tsopt --thermo \
    -b aimnet2
```

Default model: `aimnet2`. Inspect:

```bash
python -c "import mlmm.core.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

AIMNet2 honors the keys below. (Note: `_AIMNet2Backend.__init__` takes only `aimnet2_model` / `model_charge` / `model_mult` / `ml_device`; the remaining keys are calculator/config-level `MLMM_CALC_KW`, defaults in `core/defaults.py`.)

| Key | Purpose |
|---|---|
| `model_charge`, `model_mult` | Total charge and spin multiplicity |
| `ml_device` | `'cuda'`, `'cpu'`, `'auto'` |
| `aimnet2_model` | Override the default checkpoint |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `H_double` | Standard cross-backend |

## When to use AIMNet2

| Use it when | Don't use it when |
|---|---|
| The installed checkpoint covers the target elements, charge, and multiplicity | The checkpoint contract excludes any target state |
| Target-system validation meets the required error tolerance | Energies, forces, or frequencies fail the target-system comparison |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `KeyError` on element during atom-type lookup | The installed checkpoint does not support the element; select a checkpoint/backend that does. |
| `RuntimeError: charge mismatch` | AIMNet2 charge is a per-atom-network output; supply `-q TOTAL` matching the cluster. |
| Unexpected behavior for a charge or multiplicity | Confirm that state is within the installed checkpoint's documented domain and compare against a reference method. |

## See also

- `env-cuda.md` — torch / CUDA prereq.
- `core.md` — `mlmm-toolkit` install.
- `uma.md` — UMA backend setup.
- `mlmm-structure-io/charge-multiplicity.md` — figuring out
  `-q` and `-m` for an unfamiliar substrate.
