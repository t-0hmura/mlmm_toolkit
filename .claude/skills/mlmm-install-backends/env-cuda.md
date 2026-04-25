# CUDA + PyTorch setup (env-cuda.md)

This file picks up after `mlmm-env-detect/SKILL.md` — i.e. you
already know your driver version, your CPU architecture, and whether
CUDA is available via `module`, system install, or conda.

## Step 1. Pick the torch CUDA index from the driver

NVIDIA's userland driver is forward-compatible: a newer driver runs older
torch wheels, but not the other way around. Pick the **highest** torch
CUDA index whose minimum driver is ≤ your driver version.

| Driver version | Torch CUDA index (`<cu_index>`) | Notes |
|---|---|---|
| ≥ 450 (CUDA 11.8 baseline) | `cu118` | Maximum back-compat, slow on 5xx-series GPUs |
| ≥ 525 (CUDA 12.0+) | `cu121` | |
| ≥ 545 (CUDA 12.3+) | `cu124` (when available) | |
| ≥ 555 (CUDA 12.5+) | `cu126` | |
| ≥ 560 (CUDA 12.6+) | `cu126` | |
| ≥ 570 (CUDA 12.9+) | `cu129` | Recent / Blackwell-class GPUs |
| no GPU | `cpu` | CPU-only inference |

Check what wheels actually exist:

```bash
curl -s https://download.pytorch.org/whl/torch/ | grep -oE 'torch-[0-9.]+\+cu[0-9]+' | sort -u | tail
```

## Step 2. Make CUDA visible to the install

Three setups, depending on what `mlmm-env-detect/SKILL.md`
reported:

**Setup A — HPC modulefile** (`module avail cuda` had hits):

```bash
module load <CUDA_MODULE>           # exact name from `module avail cuda`
nvcc --version                      # confirm
echo "$CUDA_HOME"                   # often set by the module
```

Add `module load <CUDA_MODULE>` to **every** PBS / SLURM script that uses
the GPU (see `mlmm-hpc/SKILL.md`).

**Setup B — system install** (e.g. `/usr/local/cuda` on a workstation):

```bash
export CUDA_HOME=/usr/local/cuda
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$LD_LIBRARY_PATH"
nvcc --version
```

**Setup C — install inside a conda env** (no system CUDA, no modules):

```bash
conda activate <YOUR_ENV>
conda install -c nvidia cuda-toolkit=<MAJOR.MINOR>   # e.g. 12.6
```

## Step 3. Install torch matching `<cu_index>`

```bash
pip install torch --index-url https://download.pytorch.org/whl/<cu_index>
```

Verify:

```bash
python -c "
import torch
print('torch    :', torch.__version__)
print('cuda     :', torch.version.cuda)
print('available:', torch.cuda.is_available())
print('device 0 :', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'cpu')
"
```

If `torch.cuda.is_available()` is `False` despite a working `nvidia-smi`,
the torch wheel does not match your driver. Pick a **smaller** `<cu_index>`
and reinstall.

## Step 4. Avoid library-loading collisions

Torch ships its own CUDA libs in `<env>/lib/python<X.Y>/site-packages/nvidia/`.
A system `LD_LIBRARY_PATH` containing an older `libcusolver`,
`libcudnn`, or `libnvrtc` shadows torch's wheel. Two options:

```bash
# Option 1 — let torch's bundled libs win:
export LD_LIBRARY_PATH="$(python -c 'import torch, os; print(os.path.dirname(torch.__file__) + "/lib")'):$LD_LIBRARY_PATH"

# Option 2 — keep system libs, disable torch's preload:
export PYTORCH_NO_CUDA_PRELOAD=1
```

Symptoms that you have this problem:
`OSError: libcusolver.so.11: cannot open shared object file`,
`Could not load symbol cublasLtCreate`,
`undefined symbol: cusparseLoggerSetCallback`.

## Step 5. CPU-only fallback

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

`mlmm` runs MLIP backends on CPU but is **much slower**
(50–200×). DFT subcommand falls back to PySCF CPU automatically — see
`dft.md`.

## Architecture quirk: aarch64

If `uname -m` reports `aarch64` (e.g. ARM-based servers, Apple Silicon
under Linux containers, some HPC nodes):

- torch wheels exist for aarch64 + CUDA on recent versions; check
  `https://download.pytorch.org/whl/torch/`.
- **`gpu4pyscf-cuda12x` is x86_64 only.** DFT must use CPU PySCF on
  aarch64 — see `dft.md`.
- UMA / Orb / MACE / AIMNet2 wheels: check the backend's PyPI page.

## Loaded-state checks for an existing env

```bash
conda activate <YOUR_ENV>
python - <<'PY'
import torch, sys
print(f"python   : {sys.version.split()[0]}")
print(f"torch    : {torch.__version__}")
print(f"cuda     : {torch.version.cuda}")
print(f"cudnn    : {torch.backends.cudnn.version()}")
print(f"available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"device 0 : {torch.cuda.get_device_name(0)} ({torch.cuda.get_device_properties(0).total_memory // 1024**3} GB)")
PY
```

This is the canonical "is my CUDA + torch healthy?" probe used everywhere.

## See also

- `core.md` — install `mlmm-toolkit` itself (after torch is healthy).
- Backend mds (`uma.md`, `mace.md`, …) — extras that piggyback on the
  torch you just installed.
- `dft.md` — `gpu4pyscf-cuda12x` install + aarch64 fallback.
