# AIMNet2 backend (aimnet2.md)

AIMNet2 (IsayevLab) is a lightweight MLIP for **small organic
molecules**. It is the cheapest backend in `mlmm` but has
narrow element coverage; check whether your substrate is supported
before committing to it.

## Element coverage

AIMNet2 supports H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I.
**No first-row transition metals**, no Mg/Ca/Mn/Fe/Co/Ni/Cu/Zn —
metalloenzymes need UMA or MACE.

## Install

```bash
pip install 'mlmm[aimnet]'         # pulls aimnet>=0.1.0
```

Or, if `mlmm` is already installed:

```bash
pip install 'aimnet>=0.1.0'
```

Confirm:

```bash
python -c "import aimnet; print('aimnet:', aimnet.__version__)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='aimnet2', charge=0, spin=1)"
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
python -c "import mlmm.defaults as d; print(d.MLMM_CALC_KW)"
```

## Backend-specific flags

AIMNet2 accepts (from `backends/__init__.py:_BACKEND_ACCEPTED_KEYS['aimnet2']`):

| Key | Purpose |
|---|---|
| `charge`, `spin` | Total charge and spin multiplicity |
| `device` | `'cuda'`, `'cpu'`, `'auto'` |
| `model` | Override the default checkpoint |
| `freeze_atoms`, `hessian_calc_mode`, `return_partial_hessian`, `hessian_double` | Standard cross-backend |

## When to use AIMNet2

| Use it when | Don't use it when |
|---|---|
| Substrate is small organics (≤ 100 atoms, no metals) | Active site contains Zn, Mg, Mn, Fe, etc. |
| You want a cheap CPU baseline run | High accuracy on TS curvature is important |
| Pre-screening many transition states | Reproducing published `wB97M-V` benchmark numbers |

## Known gotchas

| Symptom | Cause / fix |
|---|---|
| `KeyError` on element during atom-type lookup | Unsupported element; switch to UMA or MACE. |
| `RuntimeError: charge mismatch` | AIMNet2 charge is per-atom-network output; supply `-q TOTAL` matching the cluster. |
| Convergence failures on radicals | AIMNet2 is closed-shell trained. Use `-m 1` only. |

## See also

- `env-cuda.md` — torch / CUDA prereq.
- `core.md` — `mlmm` install.
- `uma.md` — recommended for anything containing metals.
- `mlmm-structure-io/charge-multiplicity.md` — figuring out
  `-q` and `-m` for an unfamiliar substrate.