---
name: mlmm-install
description: "Install recipes for mlmm-toolkit core, AmberTools, and the ML backends (UMA, Orb, MACE, AIMNet2), plus the optional DFT (PySCF/GPU4PySCF) and xTB pieces, the CUDA/PyTorch wheel choice, and a probe for an unknown machine (scheduler, GPU, CUDA, conda env). `backends.md` holds per-backend notes and the probe; `ambertools.md` covers tleap and antechamber. TRIGGER on install / setup / `pip install` / `conda env` / `ImportError` / CUDA mismatch / 'GPU not detected' / `huggingface` auth / e3nn conflict / `tleap` not found / a `Missing system prerequisite:` message (compiler, Chrome libraries), or when the compute environment is unknown. SKIP when mlmm imports cleanly and the user is running subcommands; the CLI skill covers usage."
---

# Install mlmm-toolkit

Create one conda-forge env with AmberTools and PyTorch for your driver, then `pip install mlmm-toolkit` (UMA comes with the core); MACE needs its own env.

mlmm-toolkit needs a PyTorch build that matches the NVIDIA driver, at least one
MLIP backend, AmberTools for `mlmm mm-parm` ([ambertools.md](ambertools.md)),
and a C++20 compiler for the bundled `hessian_ff` kernels. PySCF/GPU4PySCF (DFT)
and the `xtb` executable are optional. The bundled `pysisyphus` fork,
`thermoanalysis`, and `hessian_ff` install with the package; do not install
them separately. The install is done when `mlmm --version` prints the version
and the checks in [Verify the install](#verify-the-install) pass.

## Install order

1. On a new or unknown host, run the probes in
   [backends.md](backends.md#probe-the-compute-environment) first.
2. Create a conda env from conda-forge with Python 3.12, AmberTools,
   PDBFixer, a C++ compiler, cyipopt (for DMF; `pydmf` comes with
   mlmm-toolkit), and PyTorch in one command.
3. Choose the PyTorch build in that command. `nvidia-smi` shows `CUDA Version`
   at its top right, the newest CUDA the driver supports; use `cuda130`, or
   `cuda129` when it is below 13.0. Do not add a pip `torch` to this env.
   Details are in [CUDA and PyTorch](backends.md#cuda-and-pytorch).
4. Install mlmm-toolkit and headless Chrome for Plotly PNG export.
5. Accept the UMA license on Hugging Face and log in ([UMA](backends.md#uma)).
6. Add other backends only when needed, after checking the model
   domain and a pilot on the target system. MACE needs its own env
   ([MACE](backends.md#mace-separate-environment)). DFT is optional
   ([DFT](backends.md#dft-pyscf-gpu4pyscf)); skip it if you only need MLIP
   energies. xTB, in the command below, is used only by the MLIP/MM
   `--embedcharge` correction or by a custom calculator that runs xTB
   ([xTB](backends.md#xtb-point-charge-correction)).

```bash
conda create -n <YOUR_ENV> -c conda-forge --override-channels python=3.12 ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler cyipopt xtb "pytorch=2.13.0=cuda130*" cuda-nvcc -y
conda activate <YOUR_ENV>
pip install mlmm-toolkit                 # UMA via fairchem-core
plotly_get_chrome -y                     # headless Chrome; needs network
hf auth login                            # after accepting the license on facebook/UMA
```

## System packages

pip and conda do not install the libraries of Chrome, and an environment built
without `cxx-compiler` has no compiler; a bare container or cluster node often
lacks both. UMA builds code the first time it runs: a GPU run needs a C compiler
(`gcc`), a CPU run a C++ compiler (`g++`), and `hessian_ff` a C++20 compiler.
PNG export needs the operating-system libraries of the Chrome that
`plotly_get_chrome -y` downloads. The conda command above installs the
compilers; check with `gcc --version` (`g++ --version` on a CPU), then install
what is missing:

```bash
sudo apt-get update && sudo apt-get install -y gcc g++    # Ubuntu/Debian
sudo dnf install -y gcc gcc-c++                            # Rocky/RHEL
```

Without sudo, use `module load gcc` or `conda install -c conda-forge c-compiler cxx-compiler`;
ask the administrator for Chrome's libraries (`all` and `path-search` then skip
only the PNG files). The Chrome libraries are, for Ubuntu/Debian (Ubuntu 22.04 and
Debian 12: `libasound2` instead of `libasound2t64`):

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates fonts-liberation libasound2t64 libatk-bridge2.0-0 \
  libatk1.0-0 libatspi2.0-0 libcairo2 libcups2 libcurl4 libdbus-1-3 libexpat1 libgbm1 \
  libglib2.0-0 libgtk-3-0 libnspr4 libnss3 libpango-1.0-0 libudev1 libvulkan1 libx11-6 libxcb1 \
  libxcomposite1 libxdamage1 libxext6 libxfixes3 libxkbcommon0 libxrandr2 wget xdg-utils
```

and for Rocky/RHEL:

```bash
sudo dnf install -y alsa-lib at-spi2-atk at-spi2-core atk ca-certificates cairo cups-libs \
  dbus-libs expat glib2 gtk3 libX11 libXcomposite libXdamage libXext libXfixes libXrandr \
  libgcc liberation-fonts libxcb libxkbcommon mesa-libgbm nspr nss nss-util pango systemd-libs \
  vulkan-loader wget xdg-utils
```

A missing one stops the command with `Missing system prerequisite:` and prints the same
commands. Docs: [System prerequisites](../../docs/installation.md#system-prerequisites).

## Install the core

`pip install mlmm-toolkit` needs no C/C++ build. The bundled `hessian_ff`
compiles its C++ kernels on first use through `torch.utils.cpp_extension`;
Ninja comes as a dependency, and a working C++20 compiler must be on `PATH`.
The bundled packages install as separate top-level packages next to `mlmm`.

```bash
conda activate <YOUR_ENV>
pip install mlmm-toolkit                                           # core only (UMA)
pip install 'mlmm-toolkit[aimnet,dft]'       # extras as needed
```

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; `-b uma` works |
| `[orb]` | `orb-models>=0.7.0` | `-b orb`, in a separate env |
| `[aimnet]` | `aimnet>=0.2.0` | `-b aimnet2` |
| `[dft]` | PySCF, CUDA 13 GPU4PySCF and CuPy on Linux x86_64 | DFT with the `cuda130` PyTorch build |
| `[dft-cuda12]` | PySCF, CUDA 12 GPU4PySCF and CuPy on Linux x86_64 | DFT with the `cuda129` PyTorch build |
| `[dev]` | `pytest` family | Contributing |

There is no `[mace]` extra, because `mace-torch` and `fairchem-core` need
different `e3nn` versions; MACE goes in a separate env
([MACE](backends.md#mace-separate-environment)). ORB also goes in its own env
([ORB](backends.md#orb-separate-environment)). `pyproject.toml` holds the
full list of extras and pins.

Contributors install from source; `pip install -e` picks up edits without a
reinstall: `git clone https://github.com/t-0hmura/mlmm_toolkit.git mlmm && cd mlmm && pip install -e '.[aimnet,dft]'`.

## Conda env template

Replace `<...>` with the values from the probe. The template uses
`python=3.12`; mlmm-toolkit supports Python 3.12 and 3.13.
`mm-parm` also needs AmberTools ([ambertools.md](ambertools.md)).

`env_mlmm.yml` (UMA, AIMNet2, DFT, xTB):

```yaml
name: <YOUR_ENV>
channels: [conda-forge, nodefaults]
dependencies:
  - python=3.12
  - ambertools=24.8
  - "numpy>=2,<2.5"
  - pdbfixer
  - cxx-compiler
  - cyipopt                            # only for DMF
  - pytorch=2.13.0=cuda130*            # cuda129* below CUDA 13.0
  - cuda-nvcc
  - xtb                                # only for the MLIP/MM correction
  - pip
  - pip:
      - mlmm-toolkit[aimnet,dft]
```

With `cuda129*`, use `[dft-cuda12]` in place of `[dft]`; see
[CUDA and PyTorch](backends.md#cuda-and-pytorch). The MACE env is built
separately ([MACE](backends.md#mace-separate-environment)).

## Choose a backend

| `-b` | Default model | Install | Notes |
|---|---|---|---|
| `uma` (default) | `uma-s-1p2p1`; `uma-s-1p2p1-ft-d4rev`, `uma-s-1p2`, `uma-m-1p1` via `--backend-model` | Core dependency; gated weights need a Hugging Face login | fp32 by default |
| `orb` | `orbmol-v2`; `orb_v3_conservative_omol` via `--backend-model` | `[orb]` extra in a separate env | fp64 by default; `--precision fp32` is for screening |
| `mace` | `MACE-OMOL-0`; `MACE-POLAR-1-S`/`-M`/`-L` (mace-torch >= 0.3.16) via `--backend-model` | Separate env (`e3nn` conflict) | fp64 by default |
| `aimnet2` | `aimnet2`; `aimnet2-2025`, `aimnet2-nse` (open-shell), `aimnet2-pd` via `--backend-model` | `[aimnet]` extra | fp32 only; `--precision fp64` stops with an error |

Select the backend with `-b` on any ML/MM calculation command. Check each
candidate model's card for supported elements, charge, multiplicity, and
training domain, then compare energies, forces, frequencies, runtime, and
memory on a representative system. Add the `[dft]` extra when DFT//MLIP/MM
single points are needed. For a non-MLIP engine, use a custom calculator.

## Custom backend (--calc-file)

To drive the ML region with an engine that is not a built-in MLIP (GFN-xTB,
DFTB+, ORCA, Psi4, and so on), write a Python file with a `get_calculator()`
factory that returns an [ASE](https://wiki.fysik.dtu.dk/ase/) Calculator and
pass it with `--calc-file`. It selects the `custom` ML backend and overrides
`-b`.

```bash
# my_calc.py:
#   from ase.calculators.emt import EMT
#   def get_calculator(charge=0, spin=1, device="auto", **kwargs):
#       return EMT()              # swap for tblite.ase.TBLite(...) etc.
mlmm sp -i complex.pdb --parm7 system.parm7 --calc-file my_calc.py -q 0 -m 1
```

The custom calculator drives the ML region only; the MM side keeps its
`hessian_ff` or OpenMM engine and the ONIOM coupling is unchanged. It works on
every ML/MM calculation command and on `all`, which passes it to every stage.
Rename the factory with `--calc-file-func-name NAME`. Hessians use finite
differences. The full guide is in [MLIP Backends](../../docs/backends.md).

## Verify the install

```bash
mlmm --version
mlmm --help                     # subcommand list
hf auth whoami                  # Hugging Face account used for UMA downloads
python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

# backend checks (only those you installed); a subcommand --help loads no backend
python -c "import fairchem.core" && echo "uma backend OK"
python -c "import orb_models" && echo "orb backend OK"
python -c "import mace"       && echo "mace backend OK"     # MACE env only
python -c "import aimnet"     && echo "aimnet2 backend OK"

# installed extras and requirements
python -c "import importlib.metadata as m; print(m.metadata('mlmm-toolkit').get_all('Provides-Extra'))"
python -c "import importlib.metadata as m; print(m.requires('mlmm-toolkit'))"
```

Success is a version string, a list of 22 subcommands (`all`, `mm-parm`,
`extract`, `path-search`, `path-opt`, `opt`, `sp`, `tsopt`, `freq`, `irc`,
`dft`, `scan`, `scan2d`, `scan3d`, `oniom-export`, `oniom-import`,
`define-layer`, `trj2fig`, `energy-diagram`, `add-elem-info`, `fix-altloc`,
`bond-summary`), your Hugging Face user name, `CUDA: True` with the GPU name,
and an `OK` line for each installed backend. An `ImportError` points to that
backend's section in [backends.md](backends.md); a CUDA error points to
[CUDA and PyTorch](backends.md#cuda-and-pytorch).

## Common failure → fix

| Symptom | Likely cause | Fix |
|---|---|---|
| `import torch` fails with `libcudart.so.12 not found` | Incomplete or mixed wheel install, or a library-path collision; not by itself evidence against the driver | [CUDA and PyTorch](backends.md#cuda-and-pytorch) |
| `e3nn` version conflict on `pip install` | UMA and MACE in the same env | Separate env for MACE ([MACE](backends.md#mace-separate-environment)) |
| `gpu4pyscf` import fails on aarch64 | `gpu4pyscf-cuda13x` is x86_64 only | Build GPU4PySCF from source or run with `--dft-engine cpu` ([DFT](backends.md#dft-pyscf-gpu4pyscf)) |
| `huggingface_hub.errors.GatedRepoError` on UMA load | License not accepted or not logged in | Accept the license on `facebook/UMA`, run `hf auth login`, confirm with `hf auth whoami` ([UMA](backends.md#uma)) |
| `OSError: libcusolver.so.11 not found` | torch's bundled CUDA libraries missing or shadowed by `LD_LIBRARY_PATH` | [CUDA and PyTorch](backends.md#cuda-and-pytorch) |
| `RuntimeError: CUDA out of memory` during `freq` | Active Hessian exceeds GPU memory | Reduce the active region, keep `return_partial_hessian: true`, and compare `Analytical` with `FiniteDifference` on a pilot ([freq](../mlmm-cli/freq.md)) |
| `tleap: command not found` | AmberTools missing or not sourced | [ambertools.md](ambertools.md) |
| `Missing system prerequisite: C compiler (gcc).` or `C++ compiler (g++).`; a PyTorch error `Failed to find C compiler` or `No working C++ compiler found` | No compiler on `PATH`; UMA builds code the first time it runs | Install `gcc` (`g++` on a CPU): [System packages](#system-packages); on a cluster, `module load gcc` in the job script |
| `Missing system prerequisite: Chrome for PNG export.`, or `path-search` / `all` finish without the diagram PNG | No headless Chrome for Plotly, or an operating-system library it needs is missing | `plotly_get_chrome -y`, then the libraries in [System packages](#system-packages) |

## Upgrade and clean rebuild

```bash
pip install --upgrade mlmm-toolkit
mlmm --version                    # confirm the new version
```

Across minor versions, also re-check `mlmm <subcommand> --help` (the flag set
may change) and the `summary.json` keys
([outputs](../mlmm-overview/outputs.md)). To remove the package, run
`pip uninstall mlmm-toolkit`, or drop the whole env with
`conda env remove -n <YOUR_ENV>`.

## When the environment is unknown

Other mlmm skills assume placeholders such as `<YOUR_QUEUE>`, `<NCPU>`,
`<NGPU>`, `<CUDA_MODULE>`, and `<YOUR_ENV>` are already known. When they are
not, for example on the first run on a new host or for an agent without prior
context, run the probes in
[Probe the compute environment](backends.md#probe-the-compute-environment).
The output fills every placeholder used by [mlmm-hpc](../mlmm-hpc/SKILL.md)
and [backends.md](backends.md).

## Next step

- Run commands: [mlmm-cli](../mlmm-cli/SKILL.md).
- Prepare structures, layers, and charges: [mlmm-model-setup](../mlmm-model-setup/SKILL.md).
- Write PBS or SLURM job scripts: [mlmm-hpc](../mlmm-hpc/SKILL.md).
- Per-backend steps, CUDA diagnostics, and environment probes: [backends.md](backends.md); AmberTools: [ambertools.md](ambertools.md).
- Docs: [Installation](../../docs/installation.md) and [MLIP Backends](../../docs/backends.md).
