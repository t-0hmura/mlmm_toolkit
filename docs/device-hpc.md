# Device Configuration & HPC Setup

How to configure GPU/CPU devices for the ML/MM calculator and submit jobs on HPC clusters.

## Device Parameters

The ML/MM calculator (`mlmm_calc.mlmm`) uses separate device settings for the ML and MM backends:

| Parameter | Default | Description |
| --- | --- | --- |
| `ml_device` | `auto` | Device for MLIP inference. `auto` selects CUDA if available, otherwise CPU. |
| `ml_cuda_idx` | `0` | CUDA device index when `ml_device=cuda`. |
| `mm_backend` | `hessian_ff` | MM force field engine. `hessian_ff` (analytical, CPU-only) or `openmm` (supports CUDA). |
| `mm_device` | `cpu` | Device for MM backend. `cpu` for hessian_ff (required). `cuda` available for openmm. |
| `mm_cuda_idx` | `0` | CUDA device index when `mm_device=cuda` (openmm only). |
| `mm_threads` | `16` | Number of CPU threads for MM backend. |

### YAML configuration example

```yaml
calc:
  ml_device: cuda
  ml_cuda_idx: 0
  mm_backend: hessian_ff
  mm_device: cpu
  mm_threads: 16
```

### Using the OpenMM backend with CUDA

```yaml
calc:
  ml_device: cuda
  ml_cuda_idx: 0
  mm_backend: openmm
  mm_device: cuda
  mm_cuda_idx: 0
```

> **Note:** When both ML and MM use CUDA, they share GPU memory. For large systems, consider using `mm_device: cpu` to reduce VRAM consumption.

---

## VRAM Management

### Post-evaluation Hessian device (`--hess-device`)

The `freq` command supports `--hess-device` to control where the evaluated Hessian is placed and diagonalized. It does not change the device used by the calculator while evaluating the Hessian:

```bash
# Default: the resolved ML device
mlmm freq -i input.pdb --parm real.parm7 -q -1

# Move the evaluated Hessian to CPU for diagonalization
mlmm freq -i input.pdb --parm real.parm7 -q -1 --hess-device cpu
```

Use `--hess-device cpu` when:
- CPU diagonalization is preferable for the evaluated Hessian
- Retaining and diagonalizing the evaluated Hessian on GPU would add avoidable VRAM pressure

This option cannot prevent an out-of-memory failure that occurs inside the backend while the Hessian is being evaluated. Reduce the active region or select a lower-memory Hessian/backend configuration for that case.

### General VRAM tips

1. **Reduce the ML region size:** Use `mlmm extract` with a smaller `--radius`. Independently, tighten `define-layer --radius-freeze` to shrink the movable-MM shell and expand the frozen environment.
2. **Use hessian_ff (default):** The hessian_ff backend runs on CPU, avoiding an additional MM allocation on the GPU.
3. **Select the MM device deliberately:** When both ML and MM use CUDA, measure memory use on a representative pilot and use `mm_device: cpu` if needed.
4. **Monitor VRAM:** `print_vram` defaults to `True` (VRAM usage is printed during Hessian computation); set `print_vram: False` in YAML to suppress it.

---

## Backend precision defaults

`--precision` selects `fp32` or `fp64` (case-insensitive). When it is omitted,
the effective default is backend-specific:

| Backend | Default | Reason |
|---|---|---|
| UMA | fp32 | Upstream fairchem baseline. |
| ORB | fp64 | Backend default. |
| MACE | fp64 | Matches MACE's upstream `default_dtype="float64"`. |
| AIMNet2 | fp32 | No precision switch; explicit fp64 is rejected. |

Validate energies, forces, frequencies, runtime, and memory for both supported
precisions on the target backend, model, and system. Precision does not replace
an independent frequency and IRC check.

```bash
# Explicit fp64 UMA calculation
mlmm tsopt -i ts.pdb --parm enzyme.parm7 -q 0 -m 1 -b uma --precision fp64 -o result_ts

# Explicit fp32 ORB calculation
mlmm scan -i r.pdb --parm enzyme.parm7 -q 0 -b orb --precision fp32 --scan-lists '[(1,5,1.4)]' -o result_scan
```

`--precision` is accepted on every compute subcommand (`sp`, `opt`, `tsopt`, `freq`, `irc`, `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`, `all`) and is routed per backend (UMA precision, ORB precision, MACE `default_dtype`).

```{note}
For `-b aimnet2`, `fp32` is a no-op and `fp64` is *rejected* because model inputs are cast to float32 upstream. UMA, Orb, and MACE accept fp64. `--deterministic` requests deterministic algorithms but does not by itself guarantee end-to-end bit identity; verify the target backend/model/SDK and stack — see [Reproducibility](reproducibility.md).
```

---

## HPC Job Submission

### PBS example

```bash
#!/bin/bash
#PBS -N mlmm_opt
#PBS -q default
#PBS -l nodes=1:ppn=32:gpus=1,mem=120GB,walltime=72:00:00
#PBS -o ${PBS_JOBNAME}.o${PBS_JOBID}
#PBS -e ${PBS_JOBNAME}.e${PBS_JOBID}

set -euo pipefail
hostname
cd "${PBS_O_WORKDIR}"

# hessian_ff JIT-compiles C++ kernels on first use. If the system compiler is
# missing or cannot compile in C++20 mode, load the site's compiler module here:
# module load <COMPILER_MODULE>

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
command -v g++ >/dev/null || { echo "g++ is required for hessian_ff" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
  echo "hessian_ff requires a compiler supporting PyTorch's C++20 JIT flag" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "ninja is required for hessian_ff" >&2; exit 1; }

# Run optimization
mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### Slurm example

```bash
#!/bin/bash
#SBATCH --job-name=mlmm_opt
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail
hostname

# hessian_ff JIT-compiles C++ kernels on first use. If needed:
# module load <COMPILER_MODULE>
source ~/miniconda3/etc/profile.d/conda.sh
conda activate <your-env>
command -v g++ >/dev/null || { echo "g++ is required for hessian_ff" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
  echo "hessian_ff requires a compiler supporting PyTorch's C++20 JIT flag" >&2
  exit 1
fi
command -v ninja >/dev/null || { echo "ninja is required for hessian_ff" >&2; exit 1; }

mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### Key points

- **Single GPU for ML:** ML inference runs on one GPU. Request `gpus=1` (PBS) or `--gres=gpu:1` (Slurm); request a second GPU only if you place the OpenMM MM backend on a separate CUDA device (`mm_device: cuda`, `mm_cuda_idx: 1`).
- **CPU threads:** Request enough CPUs for the MM backend (`mm_threads`, default 16). Set `ppn=32` (PBS) or `--cpus-per-task=32` (Slurm) for a safety margin.
- **Memory:** size RAM from a representative pilot and scheduler peak-memory logs.
- **CUDA runtime:** Official PyTorch wheels carry CUDA user-space libraries; a compatible NVIDIA driver is normally sufficient. Load a site CUDA toolkit only for an extension that needs it.
- **C++ compiler:** The default `hessian_ff` MM backend JIT-compiles C++ kernels on first use, independently of CUDA. Every compute node needs a C++20-capable compiler and Ninja (GCC 13.3 was validated); load a compiler module when the system `g++` is absent or too old.

### Specifying a GPU index

If you are allocated multiple GPUs or want to target a specific GPU on a multi-GPU node:

```bash
# Option A: Environment variable (affects all CUDA programs)
export CUDA_VISIBLE_DEVICES=0

# Option B: YAML configuration (mlmm-specific)
# In config.yaml:
# calc:
#   ml_cuda_idx: 0
mlmm opt -i input.pdb --parm real.parm7 -q -1 --config config.yaml
```

---

## Limitations

- **No ML multi-GPU parallelism:** ML inference runs on a single GPU. The OpenMM MM backend may use a separate CUDA device (`mm_device: cuda`, `mm_cuda_idx`); the default hessian_ff MM backend is CPU-only.
- **No distributed computing:** workflows run on one node. Configurations with
  `workers > 1` may spawn local worker processes but do not distribute across
  nodes.
- **hessian_ff is CPU-only:** the default MM backend runs on CPU; `mm_device` must be `cpu`/`auto` — `mm_device: cuda` raises a `ValueError` rather than silently falling back.

---

## See Also

- [Getting Started](getting-started.md) — Installation and CUDA setup
- [ML/MM Calculator](mlmm-calc.md) — Calculator architecture and parameters
- [YAML Reference](yaml-reference.md) — Full configuration reference
- [freq](freq.md) — `--hess-device` option details
- [Troubleshooting](troubleshooting.md) — Common error fixes
