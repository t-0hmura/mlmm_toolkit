# Orb backend (orb.md)

Orb-v3 (Orbital Materials) is a fast MLIP useful for **screening** large
candidate sets where you can tolerate slightly lower TS-region accuracy
than UMA / MACE.

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

Default model: `orb_v3_conservative_omol`. The `_conservative_` part
distinguishes this checkpoint (forces = ∇E, energy-conservative) from
the Orb paper's headline non-conservative variant (Neumann et al. 2024,
arXiv:2410.22570). The `_omol` suffix indicates training on the OMol25
dataset (Levine et al. 2025); element coverage matches OMol25's 83
elements (Eastman et al. 2026 benchmark, arXiv:2601.16331).

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
| 5–10× faster per call than UMA-s on > 200-atom systems | TS curvature less accurate than UMA / MACE; many TS searches end up with multiple imaginary modes |
| Trained on broad organic chemistry | Not the right tool for fine-grained `wB97M-V` benchmarking |
| Easy install, no auth gate | Smaller user community than UMA |

In practice: use Orb-v3 to **filter** down a list of candidates, then
re-run survivors with UMA or MACE for the final TS / IRC.

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| Extra small imaginary modes with explicit fp32 | Re-run at the fp64 default and independently recompute frequencies/IRC. |
| TS converges with > 1 imaginary mode | Common with Orb on aromatic or metalloenzyme systems. Re-run that step with UMA/MACE. |

## See also

- `env-cuda.md` — torch / CUDA prerequisites.
- `uma.md` — recommended primary backend for production runs.
- `mace.md` — alternative high-accuracy backend (separate env).
- `mlmm-cli/tsopt.md` — diagnosing TS convergence problems.
