# Device Configuration & HPC Setup

## Overview

> **Summary:** How to configure GPU/CPU devices for the ML/MM calculator and submit jobs on HPC clusters.

---

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

### Hessian device (`--hess-device`)

The `freq` command supports `--hess-device` to control where Hessian assembly and diagonalization run:

```bash
# Default: same device as ml_device (typically CUDA)
mlmm freq -i input.pdb --parm real.parm7 -q -1

# Force CPU for Hessian assembly (saves VRAM for large systems)
mlmm freq -i input.pdb --parm real.parm7 -q -1 --hess-device cpu
```

Use `--hess-device cpu` when:
- The active region is large (> ~500 unfrozen atoms)
- You encounter CUDA out-of-memory errors during frequency calculations
- VRAM is limited (< 16 GB)

### General VRAM tips

1. **Reduce the ML region size:** Use `mlmm extract` with a smaller `--radius` or `mlmm define-layer` with a tighter `--radius-freeze`.
2. **Use hessian_ff (default):** The hessian_ff backend is CPU-only, leaving all VRAM for MLIP inference.
3. **Avoid OpenMM CUDA for large systems:** If both ML and MM use CUDA, VRAM pressure doubles.
4. **Monitor VRAM:** Set `print_vram: True` in YAML to print VRAM usage after each ML evaluation.

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

# Load environment modules
. /home/apps/Modules/init/profile.sh
module load cuda/12.9

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mlmm

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

module load cuda/12.9

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mlmm

mlmm opt \
  -i r_complex_layered.pdb \
  --parm p_complex.parm7 \
  -q -1 -m 1 \
  --opt-mode grad \
  --out-dir opt_result
```

### Key points

- **Single GPU:** mlmm-toolkit uses one GPU per job. Request `gpus=1` (PBS) or `--gres=gpu:1` (Slurm).
- **CPU threads:** Request enough CPUs for the MM backend (`mm_threads`, default 16). Set `ppn=32` (PBS) or `--cpus-per-task=32` (Slurm) for safety margin.
- **Memory:** 120 GB is typically sufficient for enzyme active-site models. Increase for very large systems.
- **CUDA module:** Load CUDA **before** activating conda to ensure PyTorch finds the correct CUDA runtime.

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

- **No multi-GPU parallelism:** Unlike pdb2reaction (which supports multi-worker inference via Ray), mlmm-toolkit runs on a single GPU. The hessian_ff backend does not support parallel predictors.
- **No distributed computing:** All calculations run within a single process on a single node.
- **hessian_ff is CPU-only:** The default MM backend always runs on CPU regardless of `mm_device` setting.

---

## See Also

- [Getting Started](getting-started.md) -- Installation and CUDA setup
- [ML/MM Calculator](mlmm-calc.md) -- Calculator architecture and parameters
- [YAML Reference](yaml-reference.md) -- Full configuration reference
- [freq](freq.md) -- `--hess-device` option details
- [Troubleshooting](troubleshooting.md) -- Common error fixes
