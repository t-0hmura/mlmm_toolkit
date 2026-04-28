---
name: mlmm-hpc
description: How to run mlmm on PBS (Torque/PBSPro) and SLURM clusters — generic preamble templates with placeholders, walltime guidance, CPU vs GPU choice, monitoring, and the dynamic-dispatch (flock + pbsdsh) recipe in dynamic-dispatch.md.
---

# mlmm HPC

## Overview

`mlmm-toolkit` is a CPU+GPU Python program; on HPC clusters you typically
submit it as a PBS or SLURM job that requests one node with one GPU.
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
#PBS -o /dev/null
#PBS -e /dev/null
cd "${PBS_O_WORKDIR}"

# CUDA + toolchain: HPC modulefiles (env-detect outputs <CUDA_MODULE>)
# - gcc: load when the system default is too old for the CUDA toolkit or
#   when pip will compile a C/CUDA extension from source.
# - <OPENMPI_MODULE>: only for multi-node Ray (`workers > 1`); omit on
#   single-node jobs.
command -v module >/dev/null && module load <CUDA_MODULE> gcc <OPENMPI_MODULE>

# Conda env (env-detect outputs <YOUR_ENV>)
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>

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

cd "${SLURM_SUBMIT_DIR}"
# CUDA + toolchain (see PBS template above for when gcc / OpenMPI are needed)
command -v module >/dev/null && module load <CUDA_MODULE> gcc <OPENMPI_MODULE>
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate <YOUR_ENV>
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    --out-dir result_all
```

## Walltime budgeting

Empirical rough-cuts for clusters of ~200–700 atoms with UMA-s-1.1
on a single mid-range GPU. Adjust generously.

| Stage | Per-segment time | Notes |
|---|---|---|
| `extract` | < 1 min | Pure Python, CPU |
| `path-search` (GSM) | 5–30 min | Scales with `--max-nodes` |
| `path-search` (DMF) | 10–60 min | Slower than GSM but more robust |
| `tsopt` (RS-I-RFO) | 5–60 min | Hessian rebuilds dominate |
| `tsopt` (Dimer) | 1–10 min | Hessian-free; cheaper |
| `irc` | 5–30 min | Forward + backward; default 125 cycles each |
| `freq` | 5–30 min | Hessian once + diagonalization |
| `dft` (ωB97M-V/def2-svp, GPU) | 30 min – 6 h | Heavy; ~1–10 h on TZVPD |
| `dft` (CPU) | 10× GPU time | Use only for small clusters |

For an `all` run with 2 segments + DFT: budget **6–24 h walltime**.
For pure MLIP `all` (no DFT): **2–6 h** is usually enough.

## CPU vs GPU choice

| Workload | CPU | GPU |
|---|---|---|
| MLIP inference (any backend) | Possible but ~50–200× slower | **Required for production** |
| `mlmm dft` with ωB97M-V on > 200 atoms | Slow (10–100 h) | Recommended |
| `mlmm dft` with cheap functional / small molecule | Fine | Marginal speedup |
| Hessian (analytical, UMA) | OK if VRAM-limited | Faster |

Check `mlmm-install-backends/dft.md` for `--engine gpu` / `cpu`
specifics, including the aarch64 caveat (CPU PySCF only).

## Monitoring and control

PBS:

```bash
qstat -u "$USER"                 # state: Q (queued), R (running), C (complete)
qstat -f <jobid>                 # full job info
qdel <jobid>                     # cancel YOUR job
```

SLURM:

```bash
squeue -u "$USER"
scontrol show job <jobid>
scancel <jobid>
```

**`qdel` / `scancel` only kill jobs you own.** Don't kill jobs from
other users on shared clusters.

## Failed jobs / restart

`mlmm all` doesn't auto-resume by default; re-running creates
a fresh `result_all/`. Several stages support manual continuation:

- `tsopt`, `freq`, `irc`, `dft` — re-run on the previous output.
- `path-search` — pass the partial `mep.pdb` as `-i`.

For walltime-truncated jobs, write the per-stage outputs to a
persistent location and resume from the last completed stage.

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
