# MACE backend (mace.md)

MACE-OMOL-0 is an MLIP trained on the OMol25 dataset for molecular
chemistry, including biomolecules and transition-metal complexes. Validate
the selected model on representative structures and stationary points for
your system.

## Critical: separate environment required

`mace-torch` pins a **different `e3nn` version** than `fairchem-core`
(UMA). The two cannot coexist. Plan: keep UMA in your default env and
put MACE in a separate env, e.g. `<your_mace_mlmm_env>`.

```bash
conda create -n <your_mace_mlmm_env> python=3.11
conda activate <your_mace_mlmm_env>

# torch matching your CUDA driver (see env-cuda.md)
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/<cu_index>

# Install mlmm first, then replace its incompatible UMA dependency with MACE.
pip install mlmm-toolkit
pip uninstall -y fairchem-core
pip install mace-torch            # resolves MACE's required e3nn==0.4.4 last
```

Keep this order: installing `mlmm-toolkit` after MACE would pull
`fairchem-core` back in and replace MACE's `e3nn==0.4.4` with an incompatible
`e3nn>=0.5` release.

If you accidentally install both UMA and MACE in one env, you'll see
errors like:

```
ImportError: e3nn 0.5.x requires ... but mace-torch installed e3nn 0.4.x
```

The fix is to remove the env and start over (`conda env remove -n <env>`).

## Confirm install

```bash
python -c "import mace; print('mace:', mace.__version__)" && echo "mlmm + mace backend OK"
```

## CLI usage

```bash
conda activate <your_mace_mlmm_env>
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b mace
```

Default model: `MACE-OMOL-0`. Inspect:

```bash
python -c "import mlmm.core.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

MACE accepts (the `_MACEBackend.__init__` parameters in `backends/mlmm_calc.py`; defaults in `core/defaults.py`):

| Key | Purpose |
|---|---|
| `model_charge`, `model_mult` | Total charge and spin multiplicity |
| `ml_device` | `'cuda'`, `'cpu'`, `'auto'` |
| `mace_model` | Override the default MACE checkpoint |
| `mace_dtype` | `'float64'` (default; matches mlmm's float64 Hessian assembly) or `'float32'` |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `H_double` | Standard cross-backend |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| Broad molecular and elemental coverage | Separate env needed |
| Analytical/native Hessian support in mlmm | Runtime and memory depend on system, device, and precision |
| Configurable model and precision | No multi-GPU sharding API |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` import error | UMA + MACE in the same env. Use a fresh env. |
| `RuntimeError: Expected all tensors to be on the same device` | Mixed `cpu`/`cuda` tensors after a `.to()` round-trip. Restart Python and ensure `device='cuda'` consistently. |
| Slow Hessian on `mace_dtype='float64'` | Float64 is usually more expensive than float32; benchmark both on the target system and retain the precision needed for stable curvature. |

## See also

- `env-cuda.md` — torch + CUDA prereq.
- `core.md` — `mlmm-toolkit` install (do this **inside** `<your_mace_mlmm_env>`).
- `uma.md` — alternate backend; keep it in a separate env.
- `mlmm-cli/tsopt.md` — TS solver choice (Dimer vs RS-I-RFO) interacts with backend.
