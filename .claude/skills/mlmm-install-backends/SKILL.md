---
name: mlmm-install-backends
description: How to install mlmm itself and each of its MLIP / DFT / xTB backends, and how to handle CUDA, PyTorch, and platform-specific quirks. Read this before trying to import mlmm or invoke any subcommand.
---

# mlmm — Install and Backends

## Overview

`mlmm-toolkit` is a Python package that depends on:

- a recent **PyTorch** wheel matching your CUDA driver,
- one or more **MLIP backends** (UMA / Orb / MACE / AIMNet2),
- optional **PySCF / GPU4PySCF** for DFT single points,
- optional **xtb** for ALPB solvent corrections,
- **AmberTools** for MM parameterization (`mlmm mm-parm`).

Bundled and installed automatically with the package: `pysisyphus` (a
GPU-tensor fork), `thermoanalysis`. They must not be installed separately.

This skill directory contains ten files; read them in this order:

| File | When to consult |
|---|---|
| `SKILL.md` (this file) | Orientation, decision tree, install order |
| `env-cuda.md` | After confirming you have GPU + driver, before installing torch |
| `core.md` | Installing `mlmm-toolkit` itself |
| `uma.md` | Installing UMA (default backend; HuggingFace auth required) |
| `orb.md` | Installing Orb-v3 |
| `mace.md` | Installing MACE — **separate environment** required |
| `aimnet2.md` | Installing AIMNet2 |
| `ambertools.md` | AmberTools (tleap / antechamber) for `mm-parm` |
| `dft.md` | PySCF + GPU4PySCF (handled separately per `[dft]` extra) |
| `xtb.md` | xTB ALPB solvent correction layer |

## Install order

1. **Check the env** — see `mlmm-env-detect/SKILL.md` to discover
   driver / scheduler / CUDA / conda before doing anything.
2. **CUDA + torch** — `env-cuda.md` decides which torch wheel to pull.
3. **`mlmm-toolkit` core** — `core.md`.
4. **At least one MLIP backend** — start with UMA (`uma.md`); add others
   only as you need them.
5. **DFT (optional)** — `dft.md`. Skip if you only need MLIP energies.
6. **xtb (optional)** — `xtb.md`. Skip unless you need implicit solvent.

## Decision tree: which backend?

```
Need TS + IRC validation on a known organic + 1st-row metal cluster?
    └── start with UMA-s-1.1 (uma.md)
        └── if accuracy is borderline, also install MACE-OMOL-0 (mace.md, separate env)

Need a fast screen across many candidates?
    └── Orb-v3 (orb.md)

Working on small organic molecules, no metals?
    └── AIMNet2 (aimnet2.md) — limited element coverage, but light

Need DFT//MLIP refinement?
    └── add dft.md regardless of MLIP choice
```

## Why two envs for MACE

`fairchem-core` (UMA) and `mace-torch` pin **different `e3nn` versions**
which cannot coexist. Solution: keep UMA in your default env (`<env_a>`)
and put MACE in a second env (`<env_b>`). Other backends (Orb, AIMNet2,
DFT, xTB) can sit in either.

`mlmm-toolkit` itself is the same code in both envs; only the calculator
plugin set differs.

`mlmm_toolkit` shares this `e3nn`/`pysisyphus` constraint and **cannot
coexist with mlmm in one env either**. Use `<mlmm_env>` if you
need both toolkits available on the same host.

## Conda env templates

Replace `<...>` with the values you discovered in `env-detect`. The
templates assume `python=3.11`; `mlmm-toolkit` requires Python ≥ 3.11.

`env_mlmm.yml` (UMA / Orb / AIMNet2 / DFT / xTB):

```yaml
name: <your_mlmm_env>
channels: [conda-forge, nvidia]
dependencies:
  - python=3.11
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch
      - mlmm[orb,aimnet,dft]   # extras: see core.md / per-backend md
      - xtb-python                     # if you need xtb.md features
```

`env_mlmm_mace.yml` (MACE only, separate env):

```yaml
name: <your_mace_mlmm_env>
channels: [conda-forge]
dependencies:
  - python=3.11
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch
      - mace-torch
      - mlmm       # without [orb,aimnet] to avoid fairchem deps
```

`<cu_index>` is one of `cpu`, `cu118`, `cu121`, `cu126`, `cu129` — see
`env-cuda.md` for the driver version → index mapping.

## Verify the install

```bash
mlmm --version
mlmm --help                     # subcommand list
python -c "import mlmm.defaults as d; print(sorted(n for n in dir(d) if not n.startswith('_')))"

# backend smoke checks (only those you installed)
python -c "from mlmm.backends import create_calculator; create_calculator(backend='uma',  charge=0, spin=1)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='orb',  charge=0, spin=1)"
python -c "from mlmm.backends import create_calculator; create_calculator(backend='mace', charge=0, spin=1)"   # only in MACE env
python -c "from mlmm.backends import create_calculator; create_calculator(backend='aimnet2', charge=0, spin=1)"
```

If any of these fails with `ImportError`, the corresponding backend's
`md` page lists missing dependencies. If a CUDA-related error appears,
go back to `env-cuda.md`.

## Common failure → fix

| Symptom | Likely cause | Where to look |
|---|---|---|
| `import torch` fails with `libcudart.so.12 not found` | Torch wheel CUDA index mismatches the driver | `env-cuda.md` (driver → cu index table) |
| `e3nn` version conflict on `pip install` | UMA + MACE in the same env | Use a separate env for MACE (this file, `mace.md`) |
| `gpu4pyscf` import fails on aarch64 | `gpu4pyscf-cuda12x` is x86_64 only | `dft.md` — fall back to CPU PySCF |
| `huggingface_hub.errors.GatedRepoError` on UMA load | UMA model is gated, not authenticated | `uma.md` — `huggingface-cli login` |
| `OSError: libcusolver.so.11 not found` | torch's bundled CUDA libs missing or shadowed | `env-cuda.md` — `LD_LIBRARY_PATH` order |
| `RuntimeError: CUDA out of memory` during freq | Hessian batch too large | reduce `hessian_calc_mode` batch in `mlmm.defaults.UMA_CALC_KW` |

## Live source

`pyproject.toml` lists the canonical extras and version pins. To inspect
without opening the file:

```bash
python -c "import importlib.metadata as m; print(m.metadata('mlmm').get_all('Provides-Extra'))"
python -c "import importlib.metadata as m; print(m.requires('mlmm'))"
```