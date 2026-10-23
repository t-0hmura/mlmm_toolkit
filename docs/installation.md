# Installation

`mlmm-toolkit` is intended for Linux environments (local workstations or HPC clusters), and production runs normally use a CUDA-capable GPU. The MM part needs **AmberTools** (`tleap`) to build the topology and a **C++20 compiler** for the `hessian_ff` kernels. The conda command below installs both.

## Quick install

`nvidia-smi` shows `CUDA Version` at its top right: the newest CUDA the driver supports. The commands below install the conda-forge PyTorch for CUDA 13.0 (`cuda130`); if this number is below 13.0, replace `cuda130` with `cuda129`.

### Required

```bash
# 1) Create a conda environment with AmberTools, PDBFixer, a C++ compiler, cyipopt, and CUDA-enabled PyTorch
# 2) Install mlmm-toolkit
# 3) Install headless Chrome for Plotly static image export (PNG)
#    Downloads a Chromium binary; requires internet access.

conda create -n mlmm-toolkit -c conda-forge --override-channels python=3.12 ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler cyipopt "pytorch=2.13.0=cuda130*" cuda-nvcc -y
conda activate mlmm-toolkit

pip install mlmm-toolkit

plotly_get_chrome -y
```

Finally, log in to **Hugging Face Hub** so that UMA models can be downloaded. It needs a free HF account with a read-only token. Accept the FAIR Chemistry License v1 at <https://huggingface.co/facebook/UMA> first:

```bash
hf auth login
# or, with an access token in scripts:
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

You only need to do this once per machine / environment.

Check the installation with `mlmm --version`, which prints the installed version.

### Optional

ORB, AIMNet2, MACE, and DFT are installed in step 7 of {ref}`Step-by-step installation <step-by-step-installation>`.

(step-by-step-installation)=
## Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load a CUDA toolkit only when the site/build requires one**

    The conda-forge PyTorch of step 2 does not need a CUDA module. If a dependency must be
    built from source, use `module avail cuda` and load the compiler/toolkit
    combination documented by the cluster:

    ```bash
    module load cuda/<your-version>   # e.g. cuda/12.6 or cuda/12.9
    ```

2. **Create a conda environment with AmberTools and PyTorch**

    On a cluster with a system AmberTools module, run `module unload amber` first so that it does not conflict with the conda AmberTools. PyTorch comes from conda-forge in the same command (choose its CUDA build as in step 4); do not install a pip PyTorch into this environment.

    ```bash
    conda create -n <your-env> -c conda-forge --override-channels python=3.12 ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler cyipopt "pytorch=2.13.0=cuda130*" cuda-nvcc -y
    conda activate <your-env>
    ```

3. **cyipopt and pydmf (DMF)**
    The DMF method (`--mep-mode dmf`) in MEP search uses cyipopt, installed in step 2, and pydmf, which comes with `mlmm-toolkit` in step 5. If `--mep-mode dmf` still stops with an import error, see {ref}`Installation / environment <installation-environment-problems>`.

4. **Choose the PyTorch CUDA build**

    Step 2 installs PyTorch built for CUDA 13.0 (`cuda130`). If the `CUDA Version` of `nvidia-smi` is below 13.0, write `"pytorch=2.13.0=cuda129*"` in step 2 instead. Check GPU access in step 8.

5. **Install `mlmm-toolkit` itself and Chrome for visualization**

    ```bash
    pip install mlmm-toolkit
    plotly_get_chrome -y
    ```

    The `hessian_ff` kernels are built automatically on first use. If the build fails, see {ref}`hessian_ff build / import <hessian_ff-build--import>` for a manual rebuild.

    On aarch64, the `warp-lang` 1.18.0 wheel, which UMA loads through `fairchem-core`, needs glibc 2.35 and stops the import with `GLIBC_2.35' not found` on older systems (`ldd --version`); there, install `warp-lang==1.17.0` (`nvalchemi-toolkit-ops` needs 1.13.0 or newer).

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

    **ORB**:

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
    pip install 'mace-torch>=0.3.16'
    ```

    **DFT** (`-b dft`, `--dft`, `mlmm dft`): `[dft]` installs the CUDA 13 GPU4PySCF build on Linux x86_64; with the `cuda129` PyTorch build (step 4), install `[dft-cuda12]` instead. On aarch64, build [GPU4PySCF](https://github.com/pyscf/gpu4pyscf) from source. Without a GPU, `[dft]` still installs PySCF; run DFT with `--dft-engine cpu`.

    ```bash
    pip install "mlmm-toolkit[dft]"
    ```

    For when to use DFT/MM and how to check an MLIP/MM TS with it, see [Refine an MLIP TS with DFT](dft-backend.md). Two more extras are available: `[openmm]` (OpenMM as the MM backend) and `[pdbfixer]` (PDBFixer through pip instead of conda).

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

    If `CUDA: False`, inspect the installed PyTorch build, scheduler GPU visibility,
    driver, and environment libraries before changing versions:

    ```bash
    python -m torch.utils.collect_env
    python -m pip check
    ```

## System requirements

**OS.** Linux is recommended. Native Windows is not supported because AmberTools (`tleap`) is not available there.

**Python.** 3.12 is required.

**GPU / CUDA.** An NVIDIA GPU whose driver supports the chosen PyTorch build (see Quick install); newer GPU architectures may need a newer build. CPU-only execution works but is usually much slower.

**AmberTools and compiler.** AmberTools (`tleap`) builds the topology in `mm-parm` and `all`; a matching `--parm7` from an earlier run skips that step. The default `hessian_ff` MM backend needs a C++20 compiler (validated with GCC 13.3). PDBFixer is needed only for `mm-parm --add-h`.

**VRAM, RAM, and disk.** Memory grows with the backend, the atom count, and the Hessian mode, and the disk holds the environment, the model weights, the topology, and the trajectories and Hessians; for a very large system, run one representative calculation on the target node and watch the peak use.

## Next steps

- [Getting Started](getting-started.md) — the shortest run, and which page to read next
- [Quickstart: `mlmm all`](quickstart-all.md) — build an MEP from R and P
- [Quickstart: `mlmm all --scan-lists`](quickstart-scan.md) — build a path from one structure
- [Quickstart: TS-only mode](quickstart-tsopt.md) — optimize and check a TS candidate
- [Refine an MLIP TS with DFT](dft-backend.md) — refine and check the TS with DFT/MM
- [Common options and selectors](cli-conventions.md) — shared options, and how to give residues and atoms
- [Device Configuration & HPC Setup](device-hpc.md) — GPU settings and job scripts for clusters
- [Troubleshooting](troubleshooting.md) — common errors and what to try
