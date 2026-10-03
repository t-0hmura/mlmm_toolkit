---
name: mlmm-install-backends
description: "Install recipes for mlmm-toolkit core, AmberTools, and the ML backends (UMA, Orb, MACE, AIMNet2), plus the optional DFT (PySCF/GPU4PySCF) and xTB pieces, the CUDA/PyTorch wheel choice, and a probe for an unknown machine (scheduler, GPU, CUDA, conda env). `backends.md` holds per-backend notes and the probe; `ambertools.md` covers tleap and antechamber. TRIGGER on install / setup / `pip install` / `conda env` / `ImportError` / CUDA mismatch / 'GPU not detected' / `huggingface` auth / e3nn conflict / `tleap` not found, or when the compute environment is unknown. SKIP when mlmm imports cleanly and the user is running subcommands; the CLI skill covers usage."
---

# Install mlmm-toolkit

Install torch for your driver, then `pip install mlmm-toolkit`, then at least one ML backend; MACE needs its own env, and `mm-parm` needs AmberTools.

mlmm-toolkit needs a PyTorch wheel that matches the NVIDIA driver, at least one
MLIP backend, AmberTools for `mlmm mm-parm` ([ambertools.md](ambertools.md)),
and a C++20 compiler for the bundled `hessian_ff` kernels. PySCF/GPU4PySCF (DFT)
and the `xtb` executable are optional. The bundled `pysisyphus` fork,
`thermoanalysis`, and `hessian_ff` install with the package; do not install
them separately. The install is done when `mlmm --version` prints the version
and the checks in [Verify the install](#verify-the-install) pass.

## Install order

1. On a new or unknown host, run the probes in
   [backends.md](backends.md#probe-the-compute-environment) first.
2. Create a conda env with Python 3.12 (3.11 or newer is required; ORB needs
   3.11 or 3.12) and install AmberTools, PDBFixer, and a C++ compiler.
3. Install PyTorch. `nvidia-smi` shows `CUDA Version` at its top right, the
   newest CUDA the driver supports; choose a wheel at or below it (`cu126`,
   `cu130`, or `cu132`). `cu130` is the recommended choice. Details are in
   [CUDA and PyTorch](backends.md#cuda-and-pytorch).
4. Install mlmm-toolkit and headless Chrome for Plotly PNG export.
5. Accept the UMA license on Hugging Face and log in ([UMA](backends.md#uma)).
6. Add at least one more backend only when needed, after checking the model
   domain and a pilot on the target system. MACE needs its own env
   ([MACE](backends.md#mace-separate-environment)). DFT is optional
   ([DFT](backends.md#dft-pyscf-gpu4pyscf)); skip it if you only need MLIP
   energies. Install xTB only for the MLIP/MM `--embedcharge` correction or xTB
   through a custom calculator ([xTB](backends.md#xtb-point-charge-correction)).

```bash
conda create -n <YOUR_ENV> python=3.12 -y
conda activate <YOUR_ENV>
conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
pip install mlmm-toolkit                 # UMA via fairchem-core
plotly_get_chrome -y                     # headless Chrome; needs network
hf auth login                            # after accepting the license on facebook/UMA
```

## Install the core

`pip install mlmm-toolkit` needs no C/C++ build. The bundled `hessian_ff`
compiles its C++ kernels on first use through `torch.utils.cpp_extension`;
Ninja comes as a dependency, and a working C++20 compiler must be on `PATH`.
The bundled packages install as separate top-level packages next to `mlmm`.

```bash
conda activate <YOUR_ENV>
pip install mlmm-toolkit                                           # core only (UMA)
pip install --only-binary=dm-tree 'mlmm-toolkit[orb,aimnet,dft]'   # extras as needed
```

| Extra | Pulls in | When you need it |
|---|---|---|
| (none) | UMA via `fairchem-core`, base deps | Default; `-b uma` works |
| `[orb]` | `orb-models` (0.7 or newer on Python 3.12, 0.5.x on 3.11) | `-b orb` |
| `[aimnet]` | `aimnet>=0.2.0` | `-b aimnet2` |
| `[dft]` | PySCF, CUDA 13 GPU4PySCF and CuPy on Linux x86_64 | DFT with a `cu130` or `cu132` wheel |
| `[dft-cuda12]` | PySCF, CUDA 12 GPU4PySCF and CuPy on Linux x86_64 | DFT with a `cu126` wheel |
| `[mcp]` | `mcp[cli]>=1.29,<2` | Running the MCP server |
| `[dev]` | `pytest` family | Contributing |

There is no `[mace]` extra, because `mace-torch` and `fairchem-core` need
different `e3nn` versions; MACE goes in a separate env
([MACE](backends.md#mace-separate-environment)). `pyproject.toml` holds the
full list of extras and pins.

Contributors install from source; `pip install -e` picks up edits without a
reinstall: `git clone https://github.com/t-0hmura/mlmm_toolkit.git mlmm && cd mlmm && pip install --only-binary=dm-tree -e '.[orb,aimnet,dft]'`.

## Conda env template

Replace `<...>` with the values from the probe. The template uses
`python=3.12` for ORB; mlmm-toolkit itself requires Python 3.11 or newer.
`mm-parm` also needs AmberTools ([ambertools.md](ambertools.md)).

`env_mlmm.yml` (UMA, ORB, AIMNet2, DFT, xTB):

```yaml
name: <YOUR_ENV>
channels: [conda-forge, nvidia]
dependencies:
  - python=3.12
  - xtb                                # only for the MLIP/MM correction
  - pip
  - pip:
      - --extra-index-url https://download.pytorch.org/whl/<cu_index>
      - torch==2.13.0
      - mlmm-toolkit[orb,aimnet,dft]
```

`<cu_index>` is one of `cpu`, `cu126`, `cu130`, `cu132`, the indexes in
PyTorch's official 2.13.0 matrix; see
[CUDA and PyTorch](backends.md#cuda-and-pytorch). The MACE env is built
separately ([MACE](backends.md#mace-separate-environment)).

## Choose a backend

| `-b` | Default model | Install | Notes |
|---|---|---|---|
| `uma` (default) | `uma-s-1p2`; `uma-m-1p1` via `--backend-model` | Core dependency; gated weights need a Hugging Face login | fp32 by default |
| `orb` | `orb_v3_conservative_omol` | `[orb]` extra | fp64 by default; `--precision fp32` is for screening |
| `mace` | `MACE-OMOL-0` | Separate env (`e3nn` conflict) | fp64 by default |
| `aimnet2` | `aimnet2` | `[aimnet]` extra | fp32 only; `--precision fp64` stops with an error |

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
python -c "import fairchem"   && echo "uma backend OK"
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
| `import torch` fails with `libcudart.so.12 not found` | Wheel CUDA index and driver do not match, or mixed CUDA libraries | [CUDA and PyTorch](backends.md#cuda-and-pytorch) |
| `e3nn` version conflict on `pip install` | UMA and MACE in the same env | Separate env for MACE ([MACE](backends.md#mace-separate-environment)) |
| `gpu4pyscf` import fails on aarch64 | `gpu4pyscf-cuda13x` is x86_64 only | Build GPU4PySCF from source or run with `--dft-engine cpu` ([DFT](backends.md#dft-pyscf-gpu4pyscf)) |
| `huggingface_hub.errors.GatedRepoError` on UMA load | License not accepted or not logged in | Accept the license on `facebook/UMA`, run `hf auth login`, confirm with `hf auth whoami` ([UMA](backends.md#uma)) |
| `OSError: libcusolver.so.11 not found` | torch's bundled CUDA libraries missing or shadowed by `LD_LIBRARY_PATH` | [CUDA and PyTorch](backends.md#cuda-and-pytorch) |
| `RuntimeError: CUDA out of memory` during `freq` | Active Hessian exceeds GPU memory | Reduce the active region, keep `return_partial_hessian: true`, and compare `Analytical` with `FiniteDifference` on a pilot ([freq](../mlmm-cli/freq.md)) |
| `tleap: command not found` | AmberTools missing or not sourced | [ambertools.md](ambertools.md) |

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
- Prepare structures, layers, and charges: [mlmm-structure-io](../mlmm-structure-io/SKILL.md).
- Write PBS or SLURM job scripts: [mlmm-hpc](../mlmm-hpc/SKILL.md).
- Per-backend steps, CUDA diagnostics, and environment probes: [backends.md](backends.md); AmberTools: [ambertools.md](ambertools.md).
- Docs: [Installation](../../docs/installation.md) and [MLIP Backends](../../docs/backends.md).
