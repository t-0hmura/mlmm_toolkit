# Orb backend (orb.md)

Orb-v3 (Orbital Materials) is a fast MLIP useful for **screening** large
candidate sets where you can tolerate slightly lower TS-region accuracy
than UMA / MACE.

## Install

```bash
pip install 'mlmm[orb]'         # pulls orb-models
```

Or, if `mlmm-toolkit` is already installed:

```bash
pip install orb-models
```

Confirm:

```bash
python -c "import orb_models; print('orb_models:', orb_models.__version__)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='orb', charge=0, spin=1)"
```

Orb model weights are downloaded on first use; no separate auth required.

## CLI usage

```bash
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -b orb
```

Default model: `orb_v3_conservative_omol`. The `_conservative_` part
distinguishes this checkpoint (forces = ∇E, energy-conservative) from
the Orb paper's headline non-conservative variant (Neumann et al. 2024,
arXiv:2410.22570). The `_omol` suffix indicates training on the OMol25
dataset (Levine et al. 2025); element coverage matches OMol25's 83
elements (Eastman et al. 2025 benchmark, arXiv:2601.16331).

Inspect the default kwarg dict:

```bash
python -c "import mlmm.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

Orb accepts (canonical list in
`backends/__init__.py:_BACKEND_ACCEPTED_KEYS['orb']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default Orb checkpoint |
| `precision` | `'float32'` (default) or `'float64'` for tighter convergence |
| `compile_model` | `True` to torch-compile (faster after first call, slower start) |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double` | Same as UMA |

## Strengths and weaknesses

| Strength | Weakness |
|---|---|
| 5–10× faster per call than UMA-s on > 200-atom systems | TS curvature less accurate than UMA / MACE; many TS searches end up with multiple imaginary modes |
| Trained on broad organic chemistry | Not the right tool for fine-grained `wB97M-V` benchmarking |
| Easy install, no auth gate | Smaller user community than UMA |

In practice: use Orb-v3 to **filter** down a list of candidates, then
re-run survivors with UMA or MACE for the final TS / IRC.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `RuntimeError: ... mat1 and mat2 shapes ... ` during Hessian | `precision='float32'` precision insufficient on near-degenerate modes; try `precision='float64'`. |
| `compile_model=True` makes the first call 60+ s slow | Expected; subsequent calls are faster. Disable for short jobs. |
| TS converges with > 1 imaginary mode | Common with Orb on aromatic or metalloenzyme systems. Re-run that step with UMA/MACE. |

## See also

- `env-cuda.md` — torch / CUDA prerequisites.
- `uma.md` — recommended primary backend for production runs.
- `mace.md` — alternative high-accuracy backend (separate env).
- `mlmm-cli/tsopt.md` — diagnosing TS convergence problems.