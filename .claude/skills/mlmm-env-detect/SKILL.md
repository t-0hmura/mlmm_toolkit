---
name: mlmm-env-detect
description: Fallback skill for detecting the current compute environment (local / PBS / SLURM, CPU architecture, GPU, CUDA, conda env). Use only when the env is unknown — other skills assume placeholders are filled by the user or context.
---

# mlmm Environment Detection

## When to use this skill

Most `mlmm-toolkit` skills assume you already know your environment and
fill in placeholders like `<YOUR_QUEUE>`, `<NCPU>`, `<NGPU>`, `<CUDA_MODULE>`,
`<YOUR_ENV>` directly. **Use this skill only when those are unknown** — for
example, the first time you run on a new host, or when an automated agent
has no prior context about the cluster.

After running this discovery sequence, you should be able to populate
every placeholder used by `mlmm-hpc/SKILL.md`,
`mlmm-install-backends/{env-cuda,core}.md`, and the per-backend
install pages.

## Discovery sequence

Run these one by one. Each step only takes a second; the goal is to
populate the placeholder table at the bottom.

### 1. Scheduler

```bash
if   command -v qsub   >/dev/null; then SCHED=pbs    # Torque or PBSPro
elif command -v sbatch >/dev/null; then SCHED=slurm
else SCHED=local
fi
echo "scheduler: $SCHED"
```

### 2. Architecture and OS

```bash
uname -m                          # x86_64 / aarch64
uname -s                          # Linux / Darwin
lscpu | grep -E "^(Architecture|Model name|CPU\(s\)):"
```

`aarch64` (ARM64) means **`gpu4pyscf-cuda12x` is not available** — DFT must
fall back to CPU PySCF. Other backends (UMA / MACE / Orb / AIMNet2) work
on aarch64 if the wheels exist for your driver.

### 3. GPU

```bash
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>/dev/null \
  || echo "no NVIDIA GPU detected"
```

If no GPU: stay on CPU (`--engine cpu`, no `gpus=N` request in PBS).
If a GPU is present, note the **driver version** and **VRAM** — both
constrain which torch CUDA index and which model size you can use.

### 4. CUDA toolkit

Three places it might live (skill content does not transcribe paths):

```bash
# 4a. HPC modulefile (most common on managed clusters)
command -v module >/dev/null && module avail cuda 2>&1 | head -20

# 4b. System-wide install (apt/yum, /usr/local/cuda)
command -v nvcc && nvcc --version

# 4c. Inside a conda env
ls "$(conda info --base 2>/dev/null)/envs"/*/bin/nvcc 2>/dev/null
```

If none of the three returns anything, CUDA is not installed locally —
either install it or run CPU-only.

### 5. PBS scheduler details (when `SCHED=pbs`)

```bash
qstat -Q                          # available queues + max walltime
pbsnodes -a 2>/dev/null | grep -E "^[a-z0-9]|^ *(np|properties|gpus)" | head
qstat -u "$USER"                  # your running / queued jobs
```

The output of `pbsnodes -a` tells you per-node `np` (CPU count) and
`gpus` (GPU count) — these become `<NCPU>` and `<NGPU>` in PBS preambles.

### 6. SLURM scheduler details (when `SCHED=slurm`)

```bash
sinfo -o "%P %l %N %G"            # partition, max time, nodes, gres
scontrol show partition           # full partition table
squeue -u "$USER"
```

### 7. Conda envs

```bash
conda env list                    # all envs
# is mlmm installed somewhere?
for env in $(conda env list | awk '/^[a-zA-Z]/{print $1}'); do
  conda run -n "$env" python -c 'import mlmm; print("'"$env"': ", mlmm.__version__)' 2>/dev/null
done
```

The env that successfully imports `mlmm-toolkit` is your `<YOUR_ENV>`.
If none does, see `mlmm-install-backends/SKILL.md`.

### 8. Loaded modules

```bash
command -v module >/dev/null && module list 2>&1
```

## Mapping discovery output to placeholders

| Placeholder | How to fill it from the output above |
|---|---|
| `<YOUR_QUEUE>` | A queue from `qstat -Q` (PBS) whose `max_walltime` covers your job |
| `<YOUR_PARTITION>` | A partition from `sinfo -o "%P %l %N %G"` (SLURM) whose `TIMELIMIT` covers your job |
| `<NCPU>` | `np` from `pbsnodes -a` (PBS) or `--cpus-per-task` budget (SLURM) |
| `<NGPU>` | `gpus = N` from `pbsnodes -a` (PBS) or `--gres=gpu:N` (SLURM) |
| `<MEM>` | A safe fraction of the per-node memory: `pbsnodes -a \| grep totalmem` (PBS) or `sinfo -o "%m"` (SLURM) |
| `<CUDA_MODULE>` | A line from `module avail 2>&1 \| grep -i cuda` (naming varies: `cuda`, `cudatoolkit`, `nvhpc`, …) |
| `<YOUR_ENV>` | The conda env that imported `mlmm-toolkit` in step 7 |
| `<HH:MM:SS>` | Your estimated walltime, capped by the queue's `resources_max.walltime` |

## Recipe: full one-shot probe

Paste this into the host once. The output is enough to populate every
placeholder used by other `mlmm-*` skills.

```bash
{
  echo "=== Scheduler ==="
  command -v qsub   >/dev/null && echo "PBS"
  command -v sbatch >/dev/null && echo "SLURM"
  command -v qsub >/dev/null || command -v sbatch >/dev/null || echo "local only"

  echo
  echo "=== Architecture ==="
  uname -mrs

  echo
  echo "=== GPU ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>&1 || echo "no GPU"

  echo
  echo "=== CUDA modulefiles ==="
  command -v module >/dev/null && module avail cuda 2>&1 | head -20

  echo
  echo "=== PBS queues (if PBS) ==="
  command -v qstat >/dev/null && qstat -Q 2>/dev/null

  echo
  echo "=== SLURM partitions (if SLURM) ==="
  command -v sinfo >/dev/null && sinfo -o "%P %l %N %G" 2>/dev/null

  echo
  echo "=== Conda envs with mlmm ==="
  for env in $(conda env list 2>/dev/null | awk '/^[a-zA-Z]/{print $1}'); do
    conda run -n "$env" python -c \
      'import mlmm; print("'"$env"':", mlmm.__version__)' 2>/dev/null
  done
} 2>&1 | tee env_probe.txt
```

The resulting `env_probe.txt` is a self-contained summary you can refer
back to for the rest of your session.

## Cross-references

- `mlmm-hpc/SKILL.md` — uses `<YOUR_QUEUE>`, `<NCPU>`, `<NGPU>`,
  `<MEM>`, `<HH:MM:SS>`, `<CUDA_MODULE>`, `<YOUR_ENV>` placeholders.
- `mlmm-install-backends/env-cuda.md` — uses driver version
  and `<CUDA_MODULE>` to pick the right torch CUDA wheel.
- `mlmm-install-backends/dft.md` — uses `uname -m` to decide
  between `gpu4pyscf-cuda12x` (x86_64) and CPU PySCF (aarch64).