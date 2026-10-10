# Backends and environment

Per-backend install notes, the CUDA and PyTorch pairing, and probes for an
unknown host. The install order, verification, and failure table are in
[SKILL.md](SKILL.md).

- [UMA](#uma)
- [ORB](#orb)
- [MACE (separate environment)](#mace-separate-environment)
- [AIMNet2](#aimnet2)
- [DFT (PySCF, GPU4PySCF)](#dft-pyscf-gpu4pyscf)
- [xTB point-charge correction](#xtb-point-charge-correction)
- [CUDA and PyTorch](#cuda-and-pytorch)
- [Probe the compute environment](#probe-the-compute-environment)

## UMA

UMA (Universal Model for Atoms, Meta FAIR) is the default backend and covers
the broadest element and chemistry range of the four. It comes through
`fairchem-core`, a core dependency, so no extra is needed:

```bash
pip install mlmm-toolkit                      # fairchem-core comes along
python -c "import fairchem; print('fairchem :', fairchem.__version__)"
```

The weights are gated on Hugging Face. Accept the FAIR Chemistry License at
<https://huggingface.co/facebook/UMA>, then log in:

```bash
hf auth login               # paste a Read token from huggingface.co/settings/tokens
hf auth whoami
```

The token is cached in `~/.cache/huggingface/`; later runs and batch jobs pick
it up. Model variants are selected by name, not by separate repos.

`-b uma` is the default. Pick another variant with `--backend-model` or
`calc.uma_model` in YAML:

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -b uma
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -b uma --backend-model uma-m-1p1
```

| Model (`--backend-model`) | Notes |
|---|---|
| `uma-s-1p2` (default) | Small model |
| `uma-m-1p1` | Larger model; benchmark accuracy and cost on the target system |

Pitfalls:

- `GatedRepoError` or `401 Client Error: Unauthorized`: the token is missing
  or lacks access to `facebook/UMA`. Accept the license and re-run
  `hf auth login`.
- `e3nn` install conflict: `fairchem-core` clashes with `mace-torch`; put MACE
  in its own env ([MACE](#mace-separate-environment)).
- `--uma-workers` above 1 needs `fairchem-core[extras]` and finite-difference
  Hessians; an analytical Hessian needs one worker.
- A frequency calculation runs out of VRAM: compare Hessian modes and model
  sizes on a representative pilot, or move Hessian assembly to the CPU.
- The first call is slower than later ones: the model is downloaded once
  (cache in `~/.cache/huggingface/hub/`) and the kernels are compiled.

## ORB

ORB (`-b orb`) is an energy-conserving MLIP from `orb-models`. Use Python 3.12
(recommended; installs ORB 0.7) or 3.11 (installs ORB 0.5.x).

```bash
pip install 'mlmm-toolkit[orb]'   # pulls orb-models
pip install orb-models            # when mlmm-toolkit is already installed
python -c "import orb_models; print('orb backend OK:', orb_models.__version__)"
```

If installation fails, read the resolver error and `python -m pip check`
instead of adding unrelated PyG packages. The weights download on first use
without authentication.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -b orb
```

The default model is `orb_v3_conservative_omol` in fp64. Check the installed
`orb-models` model card for checkpoint provenance, supported elements, and
runtime requirements; dataset coverage alone does not establish what a
checkpoint can do. ORB is a conservative energy/force model with a
reduced-precision option and an easy install, but its TS and frequency
behavior must be validated on the target system. Compare candidate geometries
and frequencies against the production backend before mixing backends in one
workflow.

Pitfalls:

- `--precision fp32` selects ORB's reduced `float32-high` mode; use it for
  screening only.
- Extra imaginary modes: inspect the modes and compare the supported
  precisions on the target system before recomputing frequencies or IRC.
- A TS candidate with n_imag above 1 is not a first-order saddle point;
  tighten or restart the optimization and inspect every mode's displacement.

## MACE (separate environment)

MACE-OMOL-0 is an MLIP trained on the OMol25 dataset for molecular chemistry,
including biomolecules and transition-metal complexes. Validate it on
representative structures and stationary points for your system.

`mace-torch` pins `e3nn==0.4.4`, while `fairchem-core` (UMA) needs
`e3nn>=0.5`, so the two cannot share an env. Keep UMA in your default env and
put MACE in a second env. ORB goes in the default env (0.7 on 3.12, 0.5.x on 3.11); AIMNet2, DFT,
and xTB can sit in either. mlmm-toolkit is the same code in both envs; only
the backend set differs.

```bash
conda create -n <YOUR_MACE_ENV> python=3.11
conda activate <YOUR_MACE_ENV>

# torch matching your CUDA driver (see CUDA and PyTorch)
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/<cu_index>

# Install mlmm first, then replace its incompatible UMA dependency with MACE.
pip install mlmm-toolkit
pip uninstall -y fairchem-core
pip install mace-torch            # resolves MACE's required e3nn==0.4.4 last

python -c "import mace; print('mace:', mace.__version__)" && echo "mlmm + mace backend OK"
```

Keep this order. `fairchem-core` is a core dependency of mlmm-toolkit, so it
must be removed before MACE is installed; installing mlmm-toolkit after MACE
would pull `fairchem-core` back in and replace `e3nn==0.4.4` with an
incompatible `e3nn>=0.5`. The same env can come from the
[conda template](SKILL.md#conda-env-template) with `python=3.11` and plain
`mlmm-toolkit`; run the two swap commands above after activating it.

```bash
conda activate <YOUR_MACE_ENV>
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -b mace
```

The default model is `MACE-OMOL-0` in fp64. MACE has broad molecular and
elemental coverage, analytical Hessians in mlmm, and a configurable model and
precision. It needs the separate env, its runtime and memory depend on the
system, device, and precision, and it has no multi-GPU sharding.

Pitfalls:

- An error such as
  `ImportError: e3nn 0.5.x requires ... but mace-torch installed e3nn 0.4.x`
  means UMA and MACE share one env. Remove the env
  (`conda env remove -n <env>`) and start over.
- `RuntimeError: Expected all tensors to be on the same device`: mixed CPU and
  CUDA tensors after a `.to()` round trip. Restart Python and keep the device
  `cuda` throughout.
- Slow Hessians in fp64: fp64 usually costs more than fp32. Benchmark both on
  the target system and keep the precision that gives stable curvature.

## AIMNet2

AIMNet2 (`-b aimnet2`) comes from the `aimnet` package. Element support is
checkpoint-specific: check the installed checkpoint's card for its elements,
charge, multiplicity, and system size, and do not infer support from the
package name or another model generation. Validate energies, forces,
structures, and frequencies on the target system.

```bash
pip install 'mlmm-toolkit[aimnet]'         # pulls aimnet>=0.2.0
pip install 'aimnet>=0.2.0'                # when mlmm-toolkit is already installed
python -c "import aimnet; print('aimnet:', aimnet.__version__)" && echo "mlmm + aimnet2 backend OK"
mlmm all -i R.pdb P.pdb -c 'A:LIG:301' -l 'LIG:-1' --tsopt --thermo -b aimnet2
```

The default model is `aimnet2`, fp32 only. Use it when the installed
checkpoint covers the target elements, charge, and multiplicity and the
target-system comparison meets the required error tolerance; do not use it
when the checkpoint excludes any target state or the comparison fails.

Pitfalls:

- `--precision fp64` stops with an error.
- `KeyError` on an element during atom-type lookup: the checkpoint does not
  support that element; choose a checkpoint or backend that does.
- Give the ML-region total charge `-q` and multiplicity `-m`; they are model
  inputs, not per-atom charges.
- Unexpected behavior for a charge or multiplicity: confirm the state is in
  the checkpoint's documented domain and compare against a reference method.

## DFT (PySCF, GPU4PySCF)

`mlmm dft` and `-b dft` run PySCF on the CPU or GPU4PySCF on CUDA x86_64. The
DFT pieces are optional and are not part of the default install.

```bash
pip install 'mlmm-toolkit[dft]'
# CUDA 12 systems (cu126 wheel): pip install 'mlmm-toolkit[dft-cuda12]'
```

`[dft]` pulls `pyscf>=2.13.0` and `basis-set-exchange>=0.11` on every platform,
and `gpu4pyscf-cuda13x>=1.8.1,<2` with `cupy-cuda13x>=13.6,<15` on Linux x86_64
only; `[dft-cuda12]` pulls the CUDA 12 builds instead. On aarch64
(`uname -m`), the extra installs PySCF and basis-set-exchange but no
GPU4PySCF; build GPU4PySCF from source (<https://github.com/pyscf/gpu4pyscf>)
or run with `--dft-engine cpu`.

A source build that works with Python 3.12 and a CUDA 12 toolkit module:

```bash
pip install 'pyscf>=2.13.0' pyscf-dispersion cupy-cuda12x
git clone --depth 1 --branch v1.8.1 https://github.com/pyscf/gpu4pyscf.git
cd gpu4pyscf
cmake -S gpu4pyscf/lib -B build/temp.gpu4pyscf -DCUDA_ARCHITECTURES=90-real -DBUILD_LIBXC=ON
cmake --build build/temp.gpu4pyscf -j "$(nproc)"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
```

Set `CUDA_ARCHITECTURES` to the GPU's compute capability (`90-real` for
Hopper). Build in place and use `PYTHONPATH`: the package's `setup.py`
expects the separate libxc wheel, which has no aarch64 build. The libxc step
downloads its sources, so the build node needs network access. Before
production, run one small GPU SCF in the same environment, for example water
with `wb97m-v/def2-svp`.

```bash
python -c "import pyscf; print('pyscf       :', pyscf.__version__)"
python -c "import gpu4pyscf; print('gpu4pyscf   :', gpu4pyscf.__version__)"   # only on x86_64
python -c "import cupy; print('cupy        :', cupy.__version__)"
```

Usage and the CPU or GPU choice are in [dft](../mlmm-cli/dft.md).

## xTB point-charge correction

For MLIP/MM commands, `--embedcharge` adds this correction:

```text
ΔE = E_xTB(ML + MM point charges) - E_xTB(ML)
```

Forces and Hessians use the corresponding difference. It is off by default
and expensive, because each correction runs xTB with and without the MM point
charges. Keep the ML region to roughly 200–300 atoms or fewer as a practical
guideline, then benchmark the actual system, point-charge count, and hardware;
this is not a hard atom limit. In `mlmm dft`, the same flag instead adds the
Amber MM point charges directly to the PySCF Hamiltonian, without xTB.

The MLIP/MM correction calls the standalone `xtb` executable:

```bash
conda install -c conda-forge xtb
xtb --version
```

The executable must stay on `PATH` in batch jobs. Configure the correction
under `calc` as listed in `docs/yaml-reference.md`; use `--help-advanced` for
the corresponding CLI options.

## CUDA and PyTorch

FAIR-Chem selects PyTorch 2.13. PyTorch's official 2.13.0
matrix publishes Linux/Windows wheels for `cu126`, `cu130`, `cu132`, and `cpu`.

### Pick a wheel

`nvidia-smi` shows `CUDA Version` at its top right, the newest CUDA the driver
supports; choose a CUDA wheel at or below it. `cu130` is the recommended
choice, and newer GPU architectures may need a newer wheel. Use the cluster
administrator's tested module and wheel pair when one is supplied, and `cpu`
only when no NVIDIA GPU is assigned. The decisive check is the smoke test
below.

### Driver, not toolkit

An official prebuilt PyTorch wheel carries its CUDA user-space libraries. It
needs a compatible NVIDIA driver and an allocated GPU; it does not require a
matching local CUDA toolkit, `nvcc`, `CUDA_HOME`, or a `cuda/<X.Y>` module.
Start with a clean runtime and check the driver with `nvidia-smi`.

Load a CUDA toolkit and compiler only when building a C/CUDA extension from
source, for example a source-built GPU4PySCF or an architecture not covered by
the wheels. Then use the administrator's compatible module pair and keep it in
both the build and the job environment:

```bash
module load <CUDA_MODULE>           # exact site-provided name
module load <COMPILER_MODULE>       # only when required by that toolkit
nvcc --version
```

The MM side runs on CPU threads (`calc.mm_threads`, default 16), so request
`ppn` or `--cpus-per-task` of at least that many. No cross-node MPI launcher
is involved.

### Install and check torch

```bash
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/<cu_index>
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
`CUDA_VISIBLE_DEVICES` before changing the wheel index. A CPU wheel, an
unassigned GPU, an unsupported architecture, or mixed module libraries give
the same symptom.

### Library-loading collisions

Torch wheels install CUDA libraries under `site-packages/nvidia/`. A system or
module `LD_LIBRARY_PATH` can load an incompatible `libcusolver`, `libcudnn`,
`libnvrtc`, or `libnvJitLink` first. The symptoms are
`OSError: libcusolver.so.11: cannot open shared object file`,
`Could not load symbol cublasLtCreate`, and
`undefined symbol: cusparseLoggerSetCallback`. Compare against a clean
environment:

```bash
env -u LD_LIBRARY_PATH python -c "import torch; print(torch.cuda.is_available())"
python -m pip check
```

If the clean check works, remove the conflicting module or path entry from the
job. Do not use `PYTORCH_NO_CUDA_PRELOAD`; it is not a documented PyTorch
control variable.

### CPU only

```bash
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cpu
```

MLIP backends run on the CPU but usually much slower; benchmark a
representative structure. `mlmm dft` does not switch to the CPU by itself:
with the default `--dft-engine gpu` and no working GPU stack it stops with an
error, so pass `--dft-engine cpu` (or `dft.engine: cpu` in YAML).

### aarch64

When `uname -m` reports `aarch64` (ARM servers, Apple Silicon under Linux
containers, some HPC nodes):

- Torch wheels exist for aarch64 with CUDA in recent versions; check
  <https://download.pytorch.org/whl/torch/>.
- `gpu4pyscf-cuda13x` is x86_64 only; see [DFT](#dft-pyscf-gpu4pyscf).
- For UMA, ORB, MACE, and AIMNet2 wheels, check each backend's PyPI page.
- The aarch64 wheel of `warp-lang` 1.18.0, which UMA loads through `fairchem-core`, needs glibc 2.35 and stops the import with `GLIBC_2.35' not found` on older systems (`ldd --version`); there, install `warp-lang==1.17.0` (`nvalchemi-toolkit-ops` needs 1.13.0 or newer).

### Check an existing env

This is the standard check of a CUDA and torch env:

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

## Probe the compute environment

Use these probes only when the host is unknown. The report below fills every
placeholder used by the other mlmm skills.

```bash
{
  echo "=== Scheduler ==="
  command -v qsub   >/dev/null && echo "PBS"      # Torque or PBSPro
  command -v sbatch >/dev/null && echo "SLURM"
  command -v qsub >/dev/null || command -v sbatch >/dev/null || echo "local only"

  echo; echo "=== Architecture ==="
  uname -mrs                        # x86_64 / aarch64, Linux / Darwin
  lscpu | grep -E "^(Architecture|Model name|CPU\(s\)):"

  echo; echo "=== GPU ==="
  nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv 2>&1 || echo "no GPU"

  echo; echo "=== CUDA toolkit ==="
  command -v module >/dev/null && module avail cuda 2>&1 | head -20    # HPC modulefile
  command -v nvcc && nvcc --version                                    # system install
  ls "$(conda info --base 2>/dev/null)/envs"/*/bin/nvcc 2>/dev/null    # inside a conda env

  echo; echo "=== PBS queues and nodes (if PBS) ==="
  command -v qstat >/dev/null && qstat -Qf 2>/dev/null
  pbsnodes -a 2>/dev/null | grep -E "^[a-z0-9]|^ *(np|properties|gpus)" | head

  echo; echo "=== SLURM partitions (if SLURM) ==="
  command -v sinfo >/dev/null && sinfo -o "%P %l %N %G" 2>/dev/null

  echo; echo "=== Conda envs with mlmm ==="
  for env in $(conda env list 2>/dev/null | awk '/^[a-zA-Z]/{print $1}'); do
    conda run -n "$env" python -c \
      'import mlmm; print("'"$env"':", mlmm.__version__)' 2>/dev/null
  done

  echo; echo "=== Loaded modules ==="
  command -v module >/dev/null && module list 2>&1
} 2>&1
```

Follow-up per scheduler:

```bash
qstat -u "$USER"                  # PBS: your running and queued jobs
scontrol show partition           # SLURM: full partition table
squeue -u "$USER"                 # SLURM: your jobs
```

Reading the report:

- `aarch64`: the `gpu4pyscf-cuda13x` wheel is missing, so DFT runs on CPU
  PySCF unless GPU4PySCF is built from source. The MLIP backends work where
  wheels exist for the driver.
- No GPU: the MLIP backends run on the CPU (`ml_device: auto`), `dft` needs
  `--dft-engine cpu`, and any `gpus=N` request is dropped from the PBS
  preamble.
- With a GPU, note the driver version and VRAM; they bound the torch wheel and
  the model size.
- No CUDA toolkit found does not block a prebuilt CUDA wheel, which carries its
  runtime libraries and needs only the driver. Load a toolkit only to compile a
  CUDA extension from source ([CUDA and PyTorch](#cuda-and-pytorch)).
- `pbsnodes -a` gives each node's `np` (CPU count) and `gpus` (GPU count).
- The env where `import mlmm` succeeds is `<YOUR_ENV>`. If none does, follow
  [Install order](SKILL.md#install-order).

| Placeholder | How to fill it |
|---|---|
| `<YOUR_QUEUE>` | A queue from `qstat -Q` (PBS) whose `resources_max.walltime` covers the job |
| `<YOUR_PARTITION>` | A partition from `sinfo -o "%P %l %N %G"` (SLURM) whose `TIMELIMIT` covers the job |
| `<NCPU>` | `np` from `pbsnodes -a` (PBS) or the `--cpus-per-task` budget (SLURM) |
| `<NGPU>` | `gpus = N` from `pbsnodes -a` (PBS) or `--gres=gpu:N` (SLURM) |
| `<MEM>` | A safe fraction of the node memory: `pbsnodes -a \| grep totalmem` (PBS) or `sinfo -o "%m"` (SLURM) |
| `<CUDA_MODULE>` | A line from `module avail 2>&1 \| grep -i cuda` (names vary: `cuda`, `cudatoolkit`, `nvhpc`); only for source builds |
| `<YOUR_ENV>` | The conda env that imported mlmm-toolkit |
| `<HH:MM:SS>` | The estimated walltime, capped by the queue's `resources_max.walltime` |

Do not save the raw report in a project or repository: it can contain private
host names, paths, scheduler policy, and env names. If you need a file, write
it under `${TMPDIR:-/tmp}` with mode `0600`, redact it, and delete it after
copying the placeholder values.

## See also

- [SKILL.md](SKILL.md): install order, verification, failure table.
- [ambertools.md](ambertools.md): AmberTools for `mm-parm`.
- [mlmm-hpc](../mlmm-hpc/SKILL.md): job scripts that use the placeholders above.
- [tsopt](../mlmm-cli/tsopt.md) and [freq](../mlmm-cli/freq.md): Hessian mode and the TS checks per backend.
- [dft](../mlmm-cli/dft.md): the `dft` command and DFT//MLIP/MM single points.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md#charge-and-multiplicity): choosing `-q` and `-m`.
- Docs: [Installation](../../docs/installation.md), [MLIP Backends](../../docs/backends.md), [Refine an MLIP TS with DFT](../../docs/dft-backend.md).
