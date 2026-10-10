---
name: mlmm-hpc
description: "PBS and SLURM submission for mlmm-toolkit: placeholder-based job templates, resource budgeting, job monitoring, UMA predictor workers, and the dynamic-dispatch recipe for many independent systems in `dynamic-dispatch.md`. TRIGGER on `qsub` / `sbatch` / walltime / GPU or CPU resources / workers / `pbsdsh` / `flock` / batch-campaign questions. SKIP for local runs, installation, or output parsing."
---

# mlmm HPC

## Purpose

`mlmm-toolkit` is a CPU+GPU Python program; on HPC clusters you typically
submit it as a PBS or SLURM job that requests one node with one GPU by default.
This skill provides **generic templates** with placeholders — fill in
your queue / module / env names from [`mlmm-install/backends.md`](../mlmm-install/backends.md#probe-the-compute-environment).

## When the env is unknown

If you don't know the cluster's queue / GPU / module configuration,
read [`mlmm-install/backends.md`](../mlmm-install/backends.md#probe-the-compute-environment) first. It walks through the
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

# Conda env (<YOUR_ENV> from backends.md, Probe the compute environment)
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

GPU DFT (`-b dft` and `mlmm dft`) grows steeply with the ML-region size: for
`wb97m-v/def2-svp` on a 16 GB consumer GPU, a 63-atom single point took about
8 min and an 87-atom one about 18 min, roughly the 2.4–2.6 power of the atom
count. ML regions of several hundred atoms need a GPU with strong FP64
throughput and 24 GB or more; time one structure on the production GPU and
extrapolate before a batch. The first SCF of each run converges on a coarse
grid first (`--scf-stepwise-grid`, on by default), which shortened it 1.4–1.9
times from about 60 atoms up; small systems can be slightly slower, so
`--no-scf-stepwise-grid` turns it off.

## CPU vs GPU choice

| Workload | CPU | GPU |
|---|---|---|
| MLIP inference | Supported; benchmark the selected backend | Backend/model support varies; benchmark the target system |
| `mlmm dft` | Supported | Supported with a compatible GPU4PySCF stack |
| Analytical MLIP Hessian | Supported by selected backends | Runtime and memory are backend/model/system dependent; compare with finite difference on a pilot |

Check [`mlmm-install/backends.md`](../mlmm-install/backends.md#dft-pyscf-gpu4pyscf) for `--dft-engine gpu` / `cpu`
specifics, including the GPU4PySCF source build on aarch64.

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

## Before and after submitting

- Run one real job of the batch first and read its log; submit the rest only after it passes.
- Check the plan without computing: `mlmm all ... --dry-run` runs the preparation and the charge and electron-parity checks, prints the plan, and skips the calculations; `mlmm sp ... --show-config` prints the merged configuration and exits.
- Before `qsub` / `sbatch`, check that no job with the same name is queued; afterwards, confirm that exactly one was created.
- Judge success from what the job wrote, not from the job leaving the queue. Write the exit code to a file from the job script (`trap 'echo "rc=$?" > "$PBS_O_WORKDIR/$PBS_JOBID.exit"' EXIT`) and set no second EXIT trap after it. A walltime kill skips the trap, so with no exit file, read the scheduler history (`qstat -x -f <jobid>` on PBSPro, `sacct -j <jobid>` on SLURM).
- Keep heavy I/O and per-job environments on node-local scratch (`$TMPDIR`, or `/var/tmp/$PBS_JOBID`). Stop when the job ID is empty, and remove only that job's directory at the end.
- Throttle large copies to a shared file system (`rsync --bwlimit=...`); many concurrent writes can fail with I/O errors on some NFS servers.

## Failed jobs / restart

For a completed MEP with failed segment post-processing, repeat the original
`all` command against the persistent `--out-dir` and add
`--resume-segment N`. Keep the original inputs, topology, layers, extraction,
path, and calculator settings; post-processing settings may change. If the MEP
itself did not complete, restart the path calculation from its endpoint
structures.

## UMA predictor workers

The default `--uma-workers 1` uses one in-process UMA predictor. `--uma-workers N`
with `N > 1` selects fairchem's `ParallelMLIPPredictUnit`; install
`fairchem-core[extras]`, request enough GPU/process resources, and set
`--uma-workers-per-node` to match the allocation. The exact multi-node launcher is
site/fairchem specific and is intentionally absent from these generic templates.

The parallel predictor exposes no autograd model. Therefore an explicit
`--hessian-calc-mode Analytical` combined with `--uma-workers > 1` is a hard error,
not a finite-difference fallback. Use `--uma-workers 1` or explicitly select
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
| `CUDA_VISIBLE_DEVICES` | Normally leave the scheduler-provided mapping unchanged. Set it manually only outside scheduler isolation or as part of a tested worker-launch scheme; device indices inside a job are local to that mapping. |
| `OMP_NUM_THREADS=<NCPU>` | Limit OpenMP threads (avoid oversubscription) |
| `MKL_NUM_THREADS=<NCPU>` | Intel MKL thread cap |
| `CUPY_CACHE_DIR`, `CUDA_CACHE_PATH` | Put the CuPy and CUDA kernel caches in a work directory when the home directory has a file-count quota |
| `LD_LIBRARY_PATH=<torch lib>:...` | Override system CUDA libs (see backends.md, CUDA and PyTorch) |

## ssh-based remote submission

Generally avoided in shared distribution skills (depends on
per-user ssh config). If your cluster requires `ssh <login> qsub`,
add that as a wrapper around the PBS / SLURM command above; do **not**
embed it inside the skill template.

## See also

- `dynamic-dispatch.md` — flock + pbsdsh template for many short tasks.
- [`mlmm-install/backends.md`](../mlmm-install/backends.md#probe-the-compute-environment) — discover queue / module / env
  values for the placeholders above.
- [`mlmm-install/backends.md`](../mlmm-install/backends.md#cuda-and-pytorch) — driver / torch CUDA
  pairing.
- `mlmm-cli/all.md` — the typical workload submitted to HPC.
