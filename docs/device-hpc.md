# Device Configuration & HPC Setup

This page sets where the ML and MM parts of the ML/MM calculator run (GPU or CPU) and gives job scripts for PBS and Slurm. By default, ML inference runs on the GPU when CUDA is available, and the MM force field (`hessian_ff`) runs on the CPU.

## Device Parameters

The devices are set in the `calc` section of the YAML file (`--config`).

| Parameter | Default | Description |
| --- | --- | --- |
| `ml_device` | `auto` | Device for ML inference: `auto`, `cuda`, or `cpu`. `auto` selects CUDA when it is available, otherwise CPU. With `--backend dft`, `--dft-engine` (`calc.dft.engine`) sets the device instead. |
| `ml_cuda_idx` | `0` | CUDA device index for ML inference on CUDA. |
| `mm_backend` | `hessian_ff` | MM engine: `hessian_ff` (CPU only) or `openmm` (CPU or CUDA). |
| `mm_device` | `cpu` | Device for the MM engine. `hessian_ff` takes `cpu` or `auto` and runs on the CPU. `openmm` also takes `cuda`, and its `auto` selects CUDA when OpenMM has a CUDA platform. |
| `mm_cuda_idx` | `0` | CUDA device index when OpenMM runs on CUDA. |
| `mm_threads` | `16` | Number of CPU threads for the MM engine. |

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

---

## VRAM Management

### Hessian device (`--hess-device`)

`freq` and `irc` take `--hess-device`: `cuda`, `cpu`, or `auto` (default), which follows `ml_device`. In `freq`, it sets where the evaluated Hessian is kept and diagonalized. In `irc`, it sets where the initial Hessian is stored and where the IRC operations run.

```bash
# Default: the ML device
mlmm freq -i r_complex_layered.pdb --parm7 real.parm7 -q -1

# Move the evaluated Hessian to CPU for diagonalization
mlmm freq -i r_complex_layered.pdb --parm7 real.parm7 -q -1 --hess-device cpu
```

Use `--hess-device cpu` when:
- keeping and diagonalizing the Hessian on the GPU would use VRAM that the calculation needs

### General VRAM tips

1. **Reduce the ML region size:** Use `mlmm extract` with a smaller `--radius`. See [Make the model smaller](model-setup.md#make-the-model-smaller).
2. **Use hessian_ff (default):** The hessian_ff backend runs on CPU, avoiding an additional MM allocation on the GPU.
3. **Monitor VRAM:** `print_vram` defaults to `true` and prints the peak VRAM usage during Hessian computation.

---

## Precision in scheduled jobs

Choose precision by backend and purpose, then measure its cost on the allocated GPU; see [MLIP Backends › Precision](backends.md#precision).

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
  --parm7 real.parm7 \
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
  --parm7 real.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### Key points

- **GPUs:** Request one GPU for ML inference (`gpus=1` for PBS, `--gres=gpu:1` for Slurm); request a second GPU only if you place the OpenMM MM backend on a separate CUDA device (`mm_device: cuda`, `mm_cuda_idx: 1`).
- **CPU threads:** Request enough CPUs for the MM backend (`mm_threads`, default 16). The examples request 32 (`ppn=32`, `--cpus-per-task=32`) as a margin.
- **Memory:** Size RAM from a representative test run and the scheduler's peak-memory log.
- **CUDA runtime:** Official PyTorch wheels carry CUDA user-space libraries; a matching NVIDIA driver is normally sufficient. Load a site CUDA toolkit only for an extension that needs it.

### Specifying a GPU index

If you are allocated multiple GPUs or want to target a specific GPU on a multi-GPU node:

```bash
# Option A: Environment variable (affects all CUDA programs).
# Under a scheduler, keep the value it sets; set it yourself only outside one.
export CUDA_VISIBLE_DEVICES=0

# Option B: YAML configuration (mlmm-specific)
# In config.yaml:
# calc:
#   ml_cuda_idx: 0
mlmm opt -i r_complex_layered.pdb --parm7 real.parm7 -q -1 --config config.yaml
```

---

## Notes

* **hessian_ff runs only on the CPU**: with the default `mm_backend: hessian_ff`, `mm_device` takes `cpu` or `auto`, and `mm_device: cuda` stops the run with an error. Use `mm_backend: openmm` to run MM on CUDA.
* **One GPU for ML inference**: with the default `--uma-workers 1`, ML inference runs on the one GPU set by `ml_cuda_idx`. For `--uma-workers` above 1, see [MLIP Backends › Workers and Hessian mode](backends.md#workers-and-hessian-mode).
* **ML and MM on one GPU share its memory**: when both use CUDA on the same device, measure the peak memory on a representative test run, and use `mm_device: cpu` for a large system.
* **`--hess-device cpu` does not prevent every out-of-memory error**: an out-of-memory error inside the backend while the Hessian is being evaluated happens before the Hessian is moved. Make the model smaller, or choose a Hessian or backend setting that uses less memory.
* **C++ compiler on every compute node**: the default `hessian_ff` MM backend JIT-compiles C++ kernels on first use, independently of CUDA. Every compute node needs a C++20-capable compiler and Ninja (GCC 13.3 was validated); load a compiler module when the system `g++` is absent or too old. The job scripts above check both.

---

## See Also

- [Installation](installation.md) — installation, CUDA, and the C++ compiler
- [ML/MM Calculator](mlmm-calc.md) — calculator architecture and parameters
- [MLIP Backends](backends.md) — precision, workers, and Hessian mode
- [YAML Reference](yaml-reference.md) — full configuration reference
- [freq](freq.md) · [irc](irc.md) — `--hess-device`
- [Troubleshooting](troubleshooting.md) — common error fixes
