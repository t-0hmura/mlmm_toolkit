# Installation

`mlmm-toolkit` is intended for Linux environments (local workstations or HPC clusters), and production runs normally use a CUDA-capable GPU. The MM part needs **AmberTools** (`tleap`) to build the topology and a **C++20 compiler** for the `hessian_ff` kernels. The conda command below installs both.

## Quick start

For PyTorch, `nvidia-smi` shows `CUDA Version` at its top right, the newest CUDA the driver supports. Choose a wheel at or below it (`cu126`, `cu130`, or `cu132`). The commands below use the recommended `cu130`.

### Required

```bash
# 1) Create a conda environment with AmberTools, PDBFixer, and a C++ compiler
# 2) Install a CUDA-enabled PyTorch build
# 3) Install mlmm-toolkit
# 4) Install headless Chrome for Plotly static image export (PNG)
#    Downloads a Chromium binary; requires internet access.

conda create -n mlmm-toolkit python=3.12 -y
conda activate mlmm-toolkit
conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
TORCH_INDEX=cu130  # recommended; or cu126 / cu132
pip install 'torch==2.13.0' --index-url "https://download.pytorch.org/whl/${TORCH_INDEX}"
pip install mlmm-toolkit
plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. It needs a free HF account with read-only token. Accept the FAIR Chemistry License v1 at <https://huggingface.co/facebook/UMA> first:

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

You only need to do this once per machine / environment. Then check the installation with `mlmm --version`.

### Optional

For DMF, also install cyipopt and pydmf ([step 3 below](#step-by-step-installation)). ORB, AIMNet2, MACE, and DFT are installed in [step 7](#step-by-step-installation).

(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load a CUDA toolkit only when the site/build requires one**

    A prebuilt PyTorch wheel does not require `nvcc`. If a dependency must be
    built from source, use `module avail cuda` and load the compiler/toolkit
    combination documented by the cluster:

    ```bash
    module load cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
    ```

2. **Create a conda environment with AmberTools**

    On a cluster with a system AmberTools module, run `module unload amber` first so that it does not conflict with the conda AmberTools.

    ```bash
    conda create -n <your-env> python=3.12 -y
    conda activate <your-env>
    conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
    ```

3. **Install cyipopt and pydmf**
    Required if you want to use the DMF method (`--mep-mode dmf`) in MEP search; neither is installed with `mlmm-toolkit`. You can skip this step if you only use GSM. If `--mep-mode dmf` still stops with an import error, see {ref}`Installation / environment <installation--environment>`.

    ```bash
    conda install -c conda-forge cyipopt -y
    pip install 'pydmf[torch]>=1.2'   # for --dmf-backend cpu only: pip install 'pydmf>=1.2'
    ```

4. **Install PyTorch with the right CUDA build**

    Recommended example (`cu130`):

    ```bash
    pip install 'torch==2.13.0' --index-url https://download.pytorch.org/whl/cu130
    ```

    The official 2.13.0 matrix also provides `cu126`, `cu132`, and `cpu`.
    Choose the wheel with the `nvidia-smi` rule in the Quick start above, then check GPU access in step 8. See [PyTorch's version matrix](https://pytorch.org/get-started/previous-versions/).

5. **Install `mlmm-toolkit` itself and Chrome for visualization**

    ```bash
    pip install mlmm-toolkit
    plotly_get_chrome -y
    ```

    The `hessian_ff` kernels are built automatically on first use. If the build fails, see {ref}`hessian_ff build / import <hessian_ff-build--import>` for a manual rebuild.

6. **Log in to Hugging Face Hub (UMA model)**

    ```bash
    hf auth login
    ```

    For license requirements and non-interactive login, see the Required section above.

    Refer to the upstream projects for additional details:

    - fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
    - Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

7. **(Optional) Install additional MLIP backends and extras**

    mlmm-toolkit uses UMA by default. For another backend, install its extra and select it with `-b/--backend` (for example, `-b orb`):

    **ORB** (requires Python 3.11 or 3.12; 3.12 recommended):

    ```bash
    pip install "mlmm-toolkit[orb]"
    ```

    **AIMNet2**:

    ```bash
    pip install "mlmm-toolkit[aimnet]"
    ```

    **MACE**: `mace-torch` pins `e3nn==0.4.4`, which conflicts with UMA's `fairchem-core`, so install it in a separate environment built with steps 2-6.

    ```bash
    pip uninstall -y fairchem-core
    pip install mace-torch
    ```

    **DFT** (`-b dft`, `--dft`, `mlmm dft`): `[dft]` installs the CUDA 13 GPU4PySCF build on Linux x86_64, for the cu130 / cu132 PyTorch wheels of step 4; with the cu126 wheel, install `[dft-cuda12]` instead. On aarch64, build [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) from source.

    ```bash
    pip install "mlmm-toolkit[dft]"
    ```

    For when to use DFT/MM and how to check an MLIP/MM TS with it, see [Refine an MLIP TS with DFT](dft-backend.md). Three more extras are available: `[openmm]` (OpenMM as the MM backend), `[mcp]` (the `mlmm-mcp` server for agent clients), and `[pdbfixer]` (PDBFixer through pip instead of conda).

8. **Verify installation**

    ```bash
    mlmm --version
    mlmm -h
    hf auth whoami
    ```

    The first line should display the installed version, the second the list of subcommands, and the third your Hugging Face user name. To verify GPU access:

    ```bash
    python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"
    ```

    If `CUDA: False`, inspect the installed wheel, scheduler GPU visibility,
    driver, and environment libraries before changing versions:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## System requirements

**OS.** Linux is recommended. Native Windows is not supported because AmberTools (`tleap`) is not available there.

**Python.** 3.12 is recommended (3.11 at minimum); the ORB backend needs 3.11 or 3.12.

**GPU / CUDA.** An NVIDIA GPU whose driver supports the chosen wheel (see Quick start); newer GPU architectures may need a newer wheel. CPU-only execution works but is usually much slower.

**AmberTools and compiler.** AmberTools (`tleap`) builds the topology in `mm-parm` and `all`; a matching `--parm7` from an earlier run skips that step. The default `hessian_ff` MM backend needs a C++20 compiler (validated with GCC 13.3). PDBFixer is needed only for `mm-parm --add-h`.

**VRAM, RAM, and disk.** Memory grows with the backend, the atom count, and the Hessian mode, and the disk holds the environment, the model weights, the topology, and the trajectories and Hessians; run one representative calculation on the target node and watch the peak use.

## Next steps

- [Getting Started](getting-started.md) — the shortest run, and which page to read next
- [Quickstart: `mlmm all`](quickstart-all.md) — build an MEP from R and P
- [Quickstart: `mlmm all --scan-lists`](quickstart-scan.md) — build a path from one structure
- [Quickstart: TS-only mode](quickstart-tsopt.md) — optimize and check a TS candidate
- [Refine an MLIP TS with DFT](dft-backend.md) — refine and check the TS with DFT/MM
- [Common options and selectors](cli-conventions.md) — shared options, and how to give residues and atoms
- [Device Configuration & HPC Setup](device-hpc.md) — GPU settings and job scripts for clusters
- [Troubleshooting](troubleshooting.md) — common errors and what to try
