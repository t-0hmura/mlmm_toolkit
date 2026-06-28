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
| `gpu` (default) | x86_64 + CUDA, > 100 atoms, ωB97M-V or hybrid functional. **Raises `ClickException` if GPU unavailable** — does **not** auto-fallback to CPU | 1–10 h on consumer GPU per single point |
| `cpu` | aarch64, no GPU, small molecule, or when you want to force CPU | 10–100× slower than GPU |

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
| `--func-basis` | `'FUNC/BASIS'` like `'wb97m-v/def2-tzvpd'` | `wb97m-v/def2-tzvpd` |
| `--engine` | `gpu` / `cpu` | `gpu` |
| `-o, --out-dir` | Output directory | `./result_dft/` |

Inspect the live default kwargs:

```bash
python -c "import mlmm.core.defaults as d; print(d.GEOM_KW_DEFAULT, d.MLMM_CALC_KW, d.DFT_KW)"
```

## Failure diagnosis

| Symptom | Likely cause | Fix |
|---|---|---|
| `OSError: libcusolver.so.11 not found` | `LD_LIBRARY_PATH` shadowing torch's bundled CUDA libs | See `env-cuda.md`, Option 1 (`LD_LIBRARY_PATH` reorder) |
| `cupy.cuda.runtime.CUDARuntimeError: invalid device ordinal` | Torch and cupy disagree on CUDA visibility | `unset CUDA_VISIBLE_DEVICES` and let GPU4PySCF use device 0 |
| `RuntimeError: CUDA out of memory` mid-SCF | Functional / basis too heavy for VRAM | Lower `grid_level`, switch to `def2-svp`, or use `--engine cpu` |
| `gpu4pyscf` import succeeds but SCF stalls at start | cuTENSOR not installed | `pip install cutensor-cu12` (optional accelerator; no longer pulled by the `[dft]` extra) |
| aarch64: `--engine gpu` requested but no `gpu4pyscf` | Architecture not supported | Raises `ClickException`; rerun with `--engine cpu` (or set `dft.engine: cpu` in YAML); expect 10× slower |

## Memory rough-cuts

| Atoms | def2-SVP / wB97M-V | def2-TZVPD / wB97M-V |
|---|---|---|
| < 200 | < 8 GB VRAM | < 16 GB |
| 200–500 | 8–24 GB | 24–48 GB |
| 500–800 | 24–48 GB | use multi-GPU or CPU; check VRAM |
| > 800 | CPU only or smaller basis | not feasible on consumer GPU |

These are empirical — actual usage depends on grid_level and basis.

## See also

- `env-cuda.md` — `LD_LIBRARY_PATH` and torch CUDA pairing.
- `mlmm-cli/dft.md` — full subcommand flag reference.
- `mlmm-workflows-output/SKILL.md` — DFT//MLIP refinement
  workflow (run `mlmm dft` after `mlmm all`).