# DFT backend — PySCF + GPU4PySCF (dft.md)

`mlmm dft` is a single-point DFT driver that re-evaluates
stationary-point energies (R / TS / IM / P) at a higher level of theory
than the MLIP used for the geometry. It runs through PySCF on CPU or
GPU4PySCF on CUDA-enabled x86_64.

The DFT backend is **optional** — it ships in the `[dft]` extras and is
not pulled by the default install.

## Install

```bash
pip install 'mlmm-toolkit[dft]'
```

This pulls (canonical pin in `pyproject.toml`):

| Package | Purpose | Platform |
|---|---|---|
| `pyscf>=2.13.0` | Reference SCF / DFT engine | All |
| `gpu4pyscf-cuda12x>=1.7.0` | CUDA acceleration of PySCF | **x86_64 only** |
| `cupy-cuda12x>=13.0,!=13.4.0` | Tensor backend for GPU4PySCF | x86_64 only |
| `basis-set-exchange>=0.11` | Programmatic basis-set lookup | All |

On `aarch64` (`uname -m`), `gpu4pyscf-cuda12x` is unavailable — the
extras install will succeed for `pyscf` and `basis-set-exchange` but
skip GPU4PySCF, leaving you on CPU PySCF.

Verify:

```bash
python -c "import pyscf; print('pyscf       :', pyscf.__version__)"
python -c "import gpu4pyscf; print('gpu4pyscf   :', gpu4pyscf.__version__)"   # only on x86_64
python -c "import cupy; print('cupy        :', cupy.__version__)"
```

## CPU vs GPU choice

| `--engine` | When to pick | Approximate cost |
|---|---|---|
| `gpu` (default) | x86_64 + a supported CUDA/GPU4PySCF stack. **Raises `ClickException` if GPU unavailable** — does **not** auto-fallback to CPU | Pilot the target system |
| `cpu` | aarch64, no supported GPU stack, or when you want to force CPU | Pilot the target system |

## CLI usage

```bash
mlmm dft -i ts.pdb --parm real.parm7 \
    -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-svp' \
    --engine gpu                  # default; use 'cpu' to force PySCF CPU
```

Common flag set:

| Flag | Purpose | Default |
|---|---|---|
| `-i, --input` | `.pdb`, or `.xyz` (with `--ref-pdb`) input | required |
| `--parm` | Amber parm7 topology for the full system | required |
| `-q, --charge` / `-l, --ligand-charge` | Total charge or per-residue mapping | `-q` or `-l` required for all inputs (charge cannot be auto-derived without one) |
| `-m, --multiplicity` | Spin multiplicity (2S+1) | 1 |
| `--func-basis` | `'FUNC/BASIS'` like `'wb97m-v/def2-tzvpd'` | `wb97m-v/def2-svp` |
| `--engine` | `gpu` / `cpu` | `gpu` |
| `-o, --out-dir` | Output directory | `./result_dft/` |

Inspect the live default kwargs:

```bash
python -c "import mlmm.core.defaults as d; print(d.GEOM_KW_DEFAULT, d.MLMM_CALC_KW, d.DFT_KW)"
```

## Failure diagnosis

| Symptom | Likely cause | Fix |
|---|---|---|
| `OSError: libcusolver.so.11 not found` | Missing CUDA libraries or a library-path conflict | Check the full error and installed CUDA packages; for path diagnostics, see `env-cuda.md`, Option 1 |
| `cupy.cuda.runtime.CUDARuntimeError: invalid device ordinal` | Requested index is outside the visible GPU set | Keep the scheduler's `CUDA_VISIBLE_DEVICES` and select a valid local index (usually 0 in a one-GPU job) |
| `RuntimeError: CUDA out of memory` mid-SCF | Calculation exceeds available VRAM | Use `--engine cpu` or a larger-memory GPU. Lowering `grid_level` or switching to `def2-svp` is also possible, but changes the calculation and requires validation |
| `gpu4pyscf` import succeeds but SCF stalls at start | Cause cannot be determined from this symptom alone | Inspect the full log. If cuTENSOR is needed, `pip install cutensor-cu12` adds it; check the [upstream CuPy/cuTENSOR compatibility guidance](https://github.com/pyscf/gpu4pyscf#installation) first |
| aarch64: `--engine gpu` requested but no `gpu4pyscf` | Architecture not supported | Raises `ClickException`; rerun with `--engine cpu` (or set `dft.engine: cpu` in YAML) |

## Resource sizing

Runtime and memory depend on elements, basis, functional, integration grid,
engine, and software stack. Run a representative pilot and size the production
job from measured peak memory and scheduler logs.

## See also

- `env-cuda.md` — `LD_LIBRARY_PATH` and torch CUDA pairing.
- `mlmm-cli/dft.md` — full subcommand flag reference.
- `mlmm-workflows-output/SKILL.md` — DFT//MLIP/MM single-point
  workflow (run `mlmm dft` after `mlmm all`).
