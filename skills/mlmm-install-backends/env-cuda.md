# CUDA + PyTorch setup (env-cuda.md)

This file picks up after `mlmm-env-detect/SKILL.md` — i.e. you
already know your driver version, your CPU architecture, and whether
CUDA is available via `module`, system install, or conda.

## Step 1. Pick an official PyTorch 2.8 wheel

`mlmm-toolkit` pins `torch~=2.8.0`. PyTorch's official 2.8.0 matrix
publishes Linux/Windows wheels for `cu126`, `cu128`, `cu129`, and `cpu`.
It does not publish a 2.8.0 wheel on `cu118`, `cu121`, or `cu124`.

Pick the index supported by the site's driver **and the GPU architecture**:

- use the cluster administrator's tested module/wheel combination when one is
  supplied;
- `cu126` is the conservative starting point for pre-Blackwell hardware;
- use `cu128` or `cu129` when the GPU architecture or a dependency explicitly
  requires it;
- use `cpu` only when no NVIDIA GPU is assigned.

Do not convert the `nvidia-smi` "CUDA Version" banner directly into a wheel
index: it reports the newest CUDA version the driver advertises, not a locally
installed toolkit. The decisive check is the smoke test in Step 3.

## Step 2. Start from the NVIDIA driver, not a toolkit module

An official prebuilt PyTorch wheel carries its CUDA user-space libraries. It
needs a compatible NVIDIA driver and an allocated GPU; it does not require a
matching local CUDA toolkit, `nvcc`, `CUDA_HOME`, or a `cuda/<X.Y>` module.
Start with a clean runtime and verify the driver:

```bash
nvidia-smi
```

Load a CUDA toolkit and compiler only when building a C/CUDA extension from
source (for example a source-built GPU4PySCF or an architecture not covered by
available wheels). In that case, use the cluster administrator's compatible
module pair and keep it in both the build and job environments:

```bash
module load <CUDA_MODULE>           # exact site-provided name
module load <COMPILER_MODULE>       # only when required by that toolkit
nvcc --version
```

`mlmm-toolkit` parallelizes the MM side over CPU threads (default
`mm_threads=16`), so request `ppn`/`--cpus-per-task` ≥ `mm_threads`. No
cross-node MPI launcher is involved.

## Step 3. Install torch matching `<cu_index>`

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/<cu_index>
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
capture `python -m torch.utils.collect_env`, `python -m pip check`, and
`CUDA_VISIBLE_DEVICES` before changing wheel indexes. A CPU wheel, an
unassigned GPU, an unsupported architecture, or mixed module libraries can
produce the same symptom.

## Step 4. Avoid library-loading collisions

Torch wheels install CUDA libraries under `site-packages/nvidia/`. A system or
module `LD_LIBRARY_PATH` can select an incompatible `libcusolver`, `libcudnn`,
`libnvrtc`, or `libnvJitLink` first. Compare against a clean environment:

```bash
env -u LD_LIBRARY_PATH python -c "import torch; print(torch.cuda.is_available())"
python -m pip check
```

If the clean check works, remove the conflicting module/path entry from the
job. Do not use `PYTORCH_NO_CUDA_PRELOAD`; it is not a documented PyTorch
control variable.

Symptoms that you have this problem:
`OSError: libcusolver.so.11: cannot open shared object file`,
`Could not load symbol cublasLtCreate`,
`undefined symbol: cusparseLoggerSetCallback`.

## Step 5. CPU-only fallback

```bash
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
```

`mlmm-toolkit` runs MLIP backends on CPU but is usually much slower; benchmark
a representative structure. For DFT (`mlmm dft`), CPU PySCF is **not** an automatic
fallback — pass `--engine cpu` (or set `dft.engine: cpu` in YAML)
explicitly when the GPU backend is unavailable; with the default
`--engine gpu` the command raises a `ClickException` rather than
silently falling back. See `dft.md`.

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
