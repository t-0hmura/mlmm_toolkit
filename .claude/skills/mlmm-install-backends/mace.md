# MACE backend (mace.md)

MACE-OMOL-0 is a high-accuracy MLIP trained on the OMol25 dataset
(83 elements, organic + biomolecules + transition-metal complexes;
Levine et al. 2025, arXiv:2505.08762). It competes with UMA on
TS-region work and is one of the four "chemical accuracy"
(MAE ≲ 1 kcal/mol) backends in the Eastman et al. 2025 15-MLIP
benchmark (arXiv:2601.16331), alongside UMA-s-1.1, UMA-m-1.1, and
Orb-v3-omol.

## Critical: separate environment required

`mace-torch` pins a **different `e3nn` version** than `fairchem-core`
(UMA). The two cannot coexist. Plan: keep UMA in your default env and
put MACE in a separate env, e.g. `<your_mace_mlmm_env>`.

```bash
conda create -n <your_mace_mlmm_env> python=3.11
conda activate <your_mace_mlmm_env>

# torch matching your CUDA driver (see env-cuda.md)
pip install torch --index-url https://download.pytorch.org/whl/<cu_index>

# MACE + mlmm (without UMA-pulling extras)
pip install mace-torch
pip install mlmm          # WITHOUT [orb] / [aimnet] which would pull fairchem
```

If you accidentally install both UMA and MACE in one env, you'll see
errors like:

```
ImportError: e3nn 0.5.x requires ... but mace-torch installed e3nn 0.4.x
```

The fix is to nuke the env and start over (`conda env remove -n <env>`).

## Confirm install

```bash
python -c "import mace; print('mace:', mace.__version__)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='mace', charge=0, spin=1)"
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
python -c "import mlmm.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

MACE accepts (from `backends/__init__.py:_BACKEND_ACCEPTED_KEYS['mace']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default MACE checkpoint |
| `default_dtype` | `'float32'` (default) or `'float64'` |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double` | Standard cross-backend |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| Often the most accurate TS curvature on first-row organometallics | Separate env needed |
| Slower than Orb but ~2× faster than UMA-m | Larger first-load (~30 s) |
| Mature, widely benchmarked | No multi-GPU sharding API |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `e3nn` import error | UMA + MACE in the same env. Use a fresh env. |
| `RuntimeError: Expected all tensors to be on the same device` | Mixed `cpu`/`cuda` tensors after a `.to()` round-trip. Restart Python and ensure `device='cuda'` consistently. |
| Slow Hessian on `default_dtype='float64'` | Expected: float64 + 600-atom Hessian is ~4× slower than float32 with marginal accuracy gain. Use float64 only when you suspect a near-degenerate eigenvalue. |

## See also

- `env-cuda.md` — torch + CUDA prereq.
- `core.md` — `mlmm-toolkit` install (do this **inside** `<your_mace_mlmm_env>`).
- `uma.md` — primary backend; keep it in a separate env.
- `mlmm-cli/tsopt.md` — TS solver choice (Dimer vs RS-I-RFO) interacts with backend.