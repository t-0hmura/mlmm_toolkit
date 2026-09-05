---
name: mlmm-hpc
description: >-
  PBS and SLURM submission guidance for mlmm-toolkit, including generic
  placeholder-based job templates, resource budgeting, monitoring, UMA
  predictor workers, and dynamic dispatch for many independent systems. Use
  for qsub, sbatch, walltime, GPU/CPU resources, workers, pbsdsh, flock, or
  batch-campaign questions. Skip local runs, installation, and output parsing.
---

# mlmm HPC

## Purpose

`mlmm-toolkit` is a CPU+GPU Python program; on HPC clusters you typically
submit it as a PBS or SLURM job that requests one node with one GPU by default.
This skill provides **generic templates** with placeholders — fill in
your queue / module / env names from `mlmm-env-detect/SKILL.md`.

## When the env is unknown

If you don't know the cluster's queue / GPU / module configuration,
read `mlmm-env-detect/SKILL.md` first. It walks through the
discovery commands (`qstat -Q`, `pbsnodes -a`, `nvidia-smi`,
`module avail cuda`, `conda env list`) and tells you how to fill the
placeholders this skill uses.

## PBS preamble template (Torque / PBSPro)

```bash
#!/usr/bin/env bash
#PBS -N <jobname>
#PBS -q <YOUR_QUEUE>
#PBS -l nodes=1:ppn=<NCPU>:gpus=<NGPU>,mem=<MEM>GB,walltime=<HH:MM:SS>
#PBS -o <jobname>.out
#PBS -e <jobname>.err
set -euo pipefail
cd "${PBS_O_WORKDIR}"

# Preflight: fail fast if the env or the CUDA driver is missing.
command -v conda >/dev/null || { echo "conda not on PATH"; exit 1; }
nvidia-smi -L >/dev/null     || { echo "no GPU visible"; exit 1; }

# Prebuilt PyTorch/backend wheels need a compatible NVIDIA driver, not a local
# CUDA toolkit module. hessian_ff still JIT-compiles C++ kernels on first use;
# if the system g++ cannot compile in C++20 mode, load <COMPILER_MODULE> here.
# Load <CUDA_MODULE> separately only for an extension that needs that toolkit.
# The default workers=1 run needs no MPI launcher.

# Conda env (env-detect outputs <YOUR_ENV>)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>
command -v g++ >/dev/null || { echo "g++ is required for hessian_ff" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
    echo "hessian_ff requires a compiler supporting PyTorch's C++20 JIT flag" >&2
    exit 1
fi
command -v ninja >/dev/null || { echo "ninja is required for hessian_ff" >&2; exit 1; }

# Optional: torch CUDA tuning
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all > mlmm.log 2>&1
```

PBSPro syntax differs slightly (`#PBS -l select=1:ncpus=<NCPU>:ngpus=<NGPU>:mem=<MEM>gb`).
Both are accepted by most modern Torque + PBSPro installations; check
`man qsub` on your cluster.

## SLURM preamble template

```bash
#!/usr/bin/env bash
#SBATCH --job-name=<jobname>
#SBATCH --partition=<YOUR_PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=<NCPU>
#SBATCH --gres=gpu:<NGPU>
#SBATCH --mem=<MEM>G
#SBATCH --time=<HH:MM:SS>
#SBATCH --output=%x.%j.out
set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
# Preflight: confirm conda + GPU before launching
command -v conda >/dev/null || { echo "ERROR: conda not on PATH"; exit 1; }
command -v nvidia-smi >/dev/null && nvidia-smi -L || echo "WARN: nvidia-smi not found; continuing"
# Prebuilt wheels need no CUDA toolkit module. hessian_ff needs a C++20-capable compiler;
# load <COMPILER_MODULE> here if the system g++ is missing or too old.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>
command -v g++ >/dev/null || { echo "g++ is required for hessian_ff" >&2; exit 1; }
if ! g++ -std=c++20 -x c++ -fsyntax-only /dev/null; then
    echo "hessian_ff requires a compiler supporting PyTorch's C++20 JIT flag" >&2
    exit 1
fi
command -v ninja >/dev/null || { echo "ninja is required for hessian_ff" >&2; exit 1; }
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all
```

## Walltime budgeting

Pilot one representative segment on the target backend and node before
requesting production walltime. Path cost scales with `--max-nodes` and
optimizer cycles; TS/frequency cost depends on Hessian mode and active degrees
of freedom; DFT cost depends strongly on elements, basis, functional, grid, and
engine. Add margin for retries and first-use compilation.

## CPU vs GPU choice

| Workload | CPU | GPU |
|---|---|---|
| MLIP inference | Supported; benchmark the selected backend | Backend/model support varies; benchmark the target system |
| `mlmm dft` | Supported | Supported with a compatible GPU4PySCF stack |
| Analytical MLIP Hessian | Supported by selected backends | Runtime and memory are backend/model/system dependent; compare with finite difference on a pilot |

Check `mlmm-install-backends/dft.md` for `--engine gpu` / `cpu`
specifics, including the aarch64 caveat (CPU PySCF only).

## Monitoring and control

PBS:

```bash
qstat -u "$USER"                 # state: Q (queued), R (running), C (complete)
qstat -f <jobid>                 # full job info
qdel <jobid>                     # cancel a single job by id
```

SLURM:

```bash
squeue -u "$USER"
scontrol show job <jobid>
scancel <jobid>
```

Before cancellation, inspect the owner, name, and state, then cancel the
specific job ID:

```bash
qstat -f <jobid> && qdel <jobid>
scontrol show job <jobid> && scancel <jobid>
```

Do not derive cancellation IDs from an unreviewed bulk pipeline; a broad
filter can cancel an unrelated job in the same account.

## Failed jobs / restart

`mlmm all` doesn't auto-resume by default; re-running creates
a fresh `result_all/`. Several stages support manual continuation:

- `tsopt`, `freq`, `irc`, `dft` — re-run on the previous output.
- `path-search` — needs **≥2** input structures (reactant/product endpoints); pass them as repeated `-i` (a lone `mep.pdb` is rejected).

For walltime-truncated jobs, write the per-stage outputs to a
persistent location and resume from the last completed stage.

## UMA predictor workers

The default `--workers 1` uses one in-process UMA predictor. `--workers N`
with `N > 1` selects fairchem's `ParallelMLIPPredictUnit`; install
`fairchem-core[extras]`, request enough GPU/process resources, and set
`--workers-per-node` to match the allocation. The exact multi-node launcher is
site/fairchem specific and is intentionally absent from these generic templates.

The parallel predictor exposes no autograd model. Therefore an explicit
`--hessian-calc-mode Analytical` combined with `--workers > 1` is a hard error,
not a finite-difference fallback. Use `--workers 1` or explicitly select
`FiniteDifference`. ORB, MACE, AIMNet2, and custom calculators do not use the
UMA worker pool.

## Parallel job submission patterns

### Fan-out (one job per task)

```bash
for ts in seg_*.pdb; do
    jobid=$(qsub -v TS="$ts" generic_dft.sh)
    echo "submitted $ts as $jobid"
done
```

Each `qsub` produces an independent PBS job; the scheduler load-balances
them.

### Dynamic dispatch (one job, N nodes pull tasks)

When you have many short tasks and want to amortize the queue wait,
use the flock + pbsdsh pattern documented in `dynamic-dispatch.md`. One
qsub grabs N nodes, each node runs a worker that pulls tasks from a
shared list with file-lock-protected counter increment.

## Useful environment variables

| Variable | Purpose |
|---|---|
| `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` | Reduce torch memory fragmentation |
| `CUDA_VISIBLE_DEVICES=0` | Restrict to a single GPU per worker |
| `OMP_NUM_THREADS=<NCPU>` | Limit OpenMP threads (avoid oversubscription) |
| `MKL_NUM_THREADS=<NCPU>` | Intel MKL thread cap |
| `LD_LIBRARY_PATH=<torch lib>:...` | Override system CUDA libs (see env-cuda.md) |

## ssh-based remote submission

Generally avoided in shared distribution skills (depends on
per-user ssh config). If your cluster requires `ssh <login> qsub`,
add that as a wrapper around the PBS / SLURM command above; do **not**
embed it inside the skill template.

## See also

- `dynamic-dispatch.md` — flock + pbsdsh template for many short tasks.
- `mlmm-env-detect/SKILL.md` — discover queue / module / env
  values for the placeholders above.
- `mlmm-install-backends/env-cuda.md` — driver / torch CUDA
  pairing.
- `mlmm-cli/all.md` — the typical workload submitted to HPC.
