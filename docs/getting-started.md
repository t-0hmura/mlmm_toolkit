# Getting Started

## Overview

`mlmm-toolkit` is a Python CLI toolkit for computing **enzymatic reaction pathways** using an **ML/MM** (Machine Learning / Molecular Mechanics) approach. It couples an MLIP (Machine Learning Interatomic Potential) backend for the reactive (ML) region with a bundled MM force field engine (`hessian_ff`) for the surrounding protein environment, using an ONIOM-like energy decomposition. The default MLIP backend is **UMA** (Meta's FAIR-Chem); alternative backends (`orb`, `mace`, `aimnet2`) can be selected via `--backend`.

In many workflows, a **single command** is enough to generate a useful **first-pass** reaction path:
```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

---
You can also run **ML/MM model setup, MEP search, TS optimization, IRC, thermochemistry, and single-point DFT** in a single run by adding `--tsopt --thermo --dft`:
```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft
```
---

Given **(i) two or more full protein-ligand PDB files** (R,..., P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `mlmm-toolkit` automatically:

- defines the **ML region** around user-defined substrates,
- generates **Amber parm7/rst7** topology files for the MM region ([`mm-parm`](mm-parm.md)),
- assigns a **3-layer ML/MM partitioning** ([`define-layer`](define-layer.md)),
- explores **minimum-energy paths (MEPs)** with path optimization methods such as the Growing String Method (GSM) and Direct Max Flux (DMF),
- _optionally_ optimizes **transition states**, runs **vibrational analysis**, **IRC calculations**, and **single-point DFT calculations**.

At the ML stage, the reactive region uses a machine learning interatomic potential (MLIP). The default backend is UMA (Meta's FAIR-Chem); alternative backends include ORB, MACE, and AIMNet2 (selected via `-b/--backend`). The MM region uses `hessian_ff`, a bundled C++ native extension that computes Amber force field energies, forces, and Hessians. The total energy follows an ONIOM-like decomposition:

```
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

where REAL is the full system, MODEL is the ML region, "high" is the selected MLIP backend (default: UMA), and "low" is `hessian_ff`.

The CLI is designed to generate **multi-step enzymatic reaction mechanisms** with minimal manual intervention. The same workflow also works for small-molecule systems. When you skip pocket extraction (omit `--center/-c` and `--ligand-charge`), you can also use `.xyz` inputs.

```{important}
- Input PDB files must already contain **hydrogen atoms**.
- When you provide multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ); otherwise an error is raised.
- Most subcommands require `--parm` (Amber topology) and `--model-pdb` (ML region definition). The [`all`](all.md) command generates these automatically.
```

```{tip}
If you are new to the project, read [Concepts & Workflow](concepts.md) first.
For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md).
If you encounter an error during setup or runtime, refer to [Troubleshooting](troubleshooting.md).
```

### CLI conventions

| Convention | Example |
|------------|---------|
| **Residue selectors** | `'SAM,GPP'` or `'A:123,B:456'` |
| **Charge mapping** | `-l 'SAM:1,GPP:-3'` |
| **Atom selectors** | `'TYR,285,CA'` or `'TYR 285 CA'` |

For full details, see [CLI Conventions](cli-conventions.md).

`path-search` naming note: the CLI subcommand is `path-search`, while the documentation filename is [`path-search.md`](path-search.md).


### Recommended tools for hydrogen addition

If your PDB lacks hydrogen atoms, use one of the following tools before running mlmm-toolkit:

| Tool | Example Command | Notes |
|------|-----------------|-------|
| **reduce** (Richardson Lab) | `reduce input.pdb > output.pdb` | Fast, widely used for crystallographic structures |
| **pdb2pqr** | `pdb2pqr --ff=AMBER input.pdb output.pqr` | Adds hydrogens and assigns partial charges |
| **Open Babel** | `obabel input.pdb -O output.pdb -h` | General-purpose cheminformatics toolkit |
| **mm-parm --add-h** | `mlmm mm-parm -i input.pdb --add-h` | Hydrogen addition via PDBFixer (through AmberTools) |

To ensure identical atom ordering across multiple PDB inputs, apply the same hydrogen-addition tool with consistent settings to all structures.

```{warning}
This software is still under development. Please use it at your own risk.
```

---

## Installation

`mlmm-toolkit` is intended for Linux environments (local workstations or HPC clusters) with a CUDA-capable GPU. Several dependencies -- notably **PyTorch**, **fairchem-core (UMA)**, **gpu4pyscf-cuda12x**, and **hessian_ff** -- expect a working CUDA installation. Alternative MLIP backends (ORB, MACE, AIMNet2) have their own optional dependencies; see the install extras below.

### Prerequisites

mlmm-toolkit uses the following components:

- **MLIP backends**: Energy, force, and Hessian calculations for the ML region. The default is UMA (fairchem-core). ORB (`pip install "mlmm-toolkit[orb]"`) and AIMNet2 (`pip install "mlmm-toolkit[aimnet2]"`) are also available. MACE (`pip install "mlmm-toolkit[mace]"`) is also available.
- **hessian_ff**: Amber force field calculations for the MM region (requires building the C++ extension).
- **AmberTools**: Automatic parm7/rst7 generation via the `mm-parm` subcommand (tleap, antechamber, parmchk2).

Refer to the upstream projects for additional details:

- fairchem / UMA: <https://github.com/facebookresearch/fairchem>, <https://huggingface.co/facebook/UMA>
- Hugging Face token & security: <https://huggingface.co/docs/hub/security-tokens>

### Quick start

Below is a minimal setup example that works on many CUDA 12.9 clusters. Adjust module names and versions to match your system. This example assumes the default GSM MEP mode (no DMF). For DMF, install cyipopt via conda first.

```bash
# 1) Install a CUDA-enabled PyTorch build
# 2) Install mlmm-toolkit
# 3) Build the hessian_ff C++ native extension
# 4) Install a headless Chrome for Plotly figure export

pip install torch --index-url https://download.pytorch.org/whl/cu129
pip install mlmm-toolkit

# Previous stable release (v0.1.1):
# pip install git+https://github.com/t-0hmura/mlmm_toolkit.git@v0.1.1

# Optional: install alternative MLIP backends
pip install "mlmm-toolkit[orb]"       # ORB backend
pip install "mlmm-toolkit[aimnet2]"   # AIMNet2 backend
pip install "mlmm-toolkit[mace]"      # MACE backend

cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
plotly_get_chrome -y
```

> **Note:** If you switch to a different runtime environment (node/container/Python/PyTorch), rebuild `hessian_ff` in that environment.  
> On many clusters, install Ninja first:
>
> ```bash
> conda install -c conda-forge ninja -y
> cd hessian_ff/native && make clean && make
> ```

Next, log in to **Hugging Face Hub** so that UMA models can be downloaded (only required when using the default `uma` backend). Either:

```bash
# Hugging Face CLI
hf auth login --token '<YOUR_ACCESS_TOKEN>' --add-to-git-credential
```

or

```bash
# Classic CLI
huggingface-cli login
```

You only need to do this once per machine / environment.

- If you want to use the Direct Max Flux (DMF) method for MEP search, create a conda environment and install cyipopt before installing mlmm-toolkit.
 ```bash
 # Create and activate a dedicated conda environment
 conda create -n mlmm python=3.11 -y
 conda activate mlmm

 # Install cyipopt (required for the DMF method in MEP search)
 conda install -c conda-forge cyipopt -y
 ```

- If you are on an HPC cluster that uses *environment modules*, load CUDA **before** installing PyTorch, like this:
 ```bash
 module load cuda/12.9
 ```

- **AmberTools** is required for the `mlmm mm-parm` subcommand (Amber topology generation). Install it separately:
 ```bash
 conda install -c conda-forge ambertools -y
 ```
 Even without AmberTools, other subcommands work if you provide `--parm` manually.

### AmberTools installation

The `mm-parm` subcommand (automatic parm7/rst7 generation) requires AmberTools. The easiest way to install it is via conda:

```bash
conda install -c conda-forge ambertools -y
```

Even without AmberTools, other subcommands work if you provide `--parm` manually.

### Step-by-step installation

If you prefer to build the environment piece by piece:

1. **Load CUDA (if you use environment modules on an HPC cluster)**

 ```bash
 module load cuda/12.9
 ```

2. **Create and activate a conda environment**

 ```bash
 conda create -n mlmm python=3.11 -y
 conda activate mlmm
 ```

3. **Install cyipopt**
 Required if you want to use the DMF method in MEP search.

 ```bash
 conda install -c conda-forge cyipopt -y
 ```

4. **Install AmberTools**
 Required for the `mlmm mm-parm` subcommand (Amber topology generation with tleap/antechamber).

 ```bash
 conda install -c conda-forge ambertools -y
 ```

5. **Install PyTorch with the right CUDA build**

 For CUDA 12.9:

 ```bash
 pip install torch --index-url https://download.pytorch.org/whl/cu129
 ```

 (You may use another compatible version if your cluster recommends it.)

6. **Install `mlmm-toolkit`**

 ```bash
 pip install mlmm-toolkit
 ```

 To install with optional MLIP backends:

 ```bash
 pip install "mlmm-toolkit[orb]"       # ORB backend
 pip install "mlmm-toolkit[aimnet2]"   # AIMNet2 backend
 pip install "mlmm-toolkit[mace]"      # MACE backend
 ```

7. **Build the `hessian_ff` C++ native extension**

 ```bash
 cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make
 ```

 This compiles the C++ code that provides fast Amber force field energy, force, and Hessian calculations.

 > **Note:** If you move to a different runtime environment, install Ninja and rebuild in that environment:
 >
 > ```bash
 > conda install -c conda-forge ninja -y
 > cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
 > ```

8. **Install Chrome for visualization**

 ```bash
 plotly_get_chrome -y
 ```

9. **Log in to Hugging Face Hub (required for UMA backend only)**

 ```bash
 huggingface-cli login
 ```

 See also:

 - <https://github.com/facebookresearch/fairchem>
 - <https://huggingface.co/facebook/UMA>
 - <https://huggingface.co/docs/hub/security-tokens>

10. **Verify installation**

 ```bash
 mlmm --version
 ```

 This should display the installed version (e.g., `0.x.y`; the exact output depends on the git tag).

---

## Multi-backend examples

The default MLIP backend is UMA. Use `-b/--backend` to switch to an alternative backend, and `--embedcharge` to enable xTB point-charge embedding:

```bash
# Use ORB backend
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b orb

# Use MACE backend
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b mace

# Enable xTB point-charge embedding
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 --embedcharge
```

---

## Quickstart routes (recommended)

- [Quickstart: run `mlmm all`](quickstart-all.md)
- [Quickstart: run `mlmm scan`](quickstart-scan-spec.md)
- [Quickstart: validate TS with `mlmm tsopt`](quickstart-tsopt-freq.md)

---

## Typical workflow

The `mlmm all` command orchestrates a multi-step pipeline. When run individually, the typical workflow is:

```text
1. extract - Define ML region from full protein-ligand PDB
2. mm-parm - Generate Amber parm7/rst7 topology (requires AmberTools)
3. define-layer - Assign 3-layer ML/MM partitioning (B-factor encoding)
4. path-search - Recursive MEP search (Growing String Method)
5. tsopt - Transition state optimization
6. freq - Vibrational analysis and thermochemistry
7. dft - Single-point DFT energy refinement
```

The `all` command runs steps 1-7 automatically. You can also run each step individually for debugging or custom workflows.

---

## Command line basics

The main entry point is the `mlmm` command, installed via `pip`. Internally it uses the **Click** library, and the default subcommand is `all`.

This is equivalent to:

```bash
mlmm [OPTIONS]...
# is equivalent to
mlmm all [OPTIONS]...
```

The `all` command runs the full pipeline -- ML region extraction, MM parameterization, layer definition, MEP search, TS optimization, vibrational analysis, and optional DFT -- in a single invocation.

All high-level workflows share two important options when you use ML region extraction:

- `-i/--input`: one or more **full structures** (reactant, intermediate(s), product).
- `-c/--center`: how to define the **substrate / extraction center** (e.g., residue names or residue IDs).

If you omit `--center/-c`, ML region extraction is skipped and the **full input structure** is used directly.

---

## Main workflow modes

### Multi-structure MEP workflow (reactant -> product)

Use this when you already have several full PDB structures along a putative reaction coordinate (e.g., R -> I1 -> I2 -> P).

**Minimal example**

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

**Richer example**

```bash
mlmm -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --out-dir ./result_all --tsopt --thermo --dft
```

Behavior:

- takes two or more **full systems** in reaction order,
- defines the ML region for each structure,
- generates Amber parm7/rst7 topology and assigns 3-layer ML/MM partitioning,
- performs a **recursive MEP search** via `path-search` by default (outputs under `path_search/`),
- optionally switches to a **single-pass** `path-opt` run with `--no-refine-path`,
- when PDB templates are available, merges the ML-region MEP back into the **full system**,
- optionally runs TS optimization, vibrational analysis, and single-point DFT calculations for each segment.

This is the recommended mode when you can generate reasonably spaced intermediates (e.g., from docking, MD, or manual modeling).

```{important}
`mlmm-toolkit` assumes that multiple input PDBs contain **exactly the same atoms in the same order** (only coordinates may differ). If any non-coordinate fields differ across inputs, an error is raised. Input PDB files must also contain **hydrogen atoms**.
```

---

### Single-structure + staged scan (feeds MEP refinement)

Use this when you prefer to define reaction coordinates yourself, rather than providing multiple endpoint structures.

Provide a single `-i` together with `--scan-lists`:

**Minimal example**

```bash
mlmm -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]'
```

**Richer example**

```bash
mlmm -i SINGLE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' '[("TYR 285 CB","MMT 309 C11",1.20)]' --multiplicity 1 --out-dir ./result_scan_all --tsopt --thermo --dft
```

Key points:

- `--scan-lists` describes **staged distance scans** on the extracted ML region.
- Each tuple `(i, j, target_A)` is:
 - a PDB atom selector string like `'TYR,285,CA'` (**delimiters can be: space/comma/slash/backtick/backslash ` ` `,` `/` `` ` `` `\`**) **or** a 1-based atom index,
 - automatically remapped to the ML region indices.
- Supplying one `--scan-lists` literal runs a single scan stage; multiple literals run sequential stages. Pass multiple literals after a single flag (repeated flags are not accepted).
- Each stage writes a `stage_XX/result.pdb`, which is treated as a candidate intermediate or product.
- The default `all` workflow refines the concatenated stages with recursive `path-search`.
- With `--no-refine-path`, it instead performs a single-pass `path-opt` chain and skips the recursive refiner.

This mode is useful for building reaction paths starting from a single structure.

---

### Single-structure TSOPT-only mode

Use this when you already have a **transition state candidate** and only want to optimize it and proceed to IRC calculations.

Provide exactly one PDB and enable `--tsopt`:

**Minimal example**

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt
```

**Richer example**

```bash
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft --out-dir ./result_tsopt_only
```

Behavior:

- skips the MEP/path search entirely,
- optimizes the **ML/MM TS** with TS optimization,
- runs an **IRC** in both directions and optimizes both ends to relax down to R and P minima,
- can then perform `freq` and `dft` on the R/TS/P,
- produces MLIP, Gibbs, and DFT//MLIP energy diagrams.

```{important}
Single-input runs require **either** `--scan-lists` (staged scan -> GSM) **or** `--tsopt` (TSOPT-only). Supplying only a single `-i` without one of these will not trigger a full workflow.
```

---

## Important CLI options and behaviors

Below are the most commonly used options across workflows.

| Option | Description |
|--------|-------------|
| `-i, --input PATH...` | Input structures. **>= 2 PDBs** -> MEP search; **1 PDB + `--scan-lists`** -> staged scan -> GSM; **1 PDB + `--tsopt`** -> TSOPT-only mode. |
| `-c, --center TEXT` | Defines the substrate / extraction center. Supports residue names (`'SAM,GPP'`), residue IDs (`A:123,B:456`), or PDB paths. |
| `-l, --ligand-charge TEXT` | Charge info: mapping (`'SAM:1,GPP:-3'`) or single integer. |
| `-q, --charge INT` | Hard override of total system charge. |
| `-m, --multiplicity INT` | Spin multiplicity (e.g., `1` for singlet). |
| `-s, --scan-lists TEXT...` | Staged distance scans for single-input runs (YAML/JSON file or inline literals). |
| `-o, --out-dir PATH` | Top-level output directory. |
| `--tsopt/--no-tsopt` | Enable TS optimization and IRC. |
| `--thermo/--no-thermo` | Run vibrational analysis and thermochemistry. |
| `--dft/--no-dft` | Perform single-point DFT calculations. |
| `--refine-path/--no-refine-path` | Recursive MEP refinement (default) vs single-pass. |
| `--mep-mode gsm\|dmf` | MEP method: Growing String Method or Direct Max Flux. | `gsm` |
| `--opt-mode grad\|hess` | Workflow preset in `all`: `grad` (LBFGS/Dimer, default) or `hess` (RFO/RS-I-RFO). |
| `-b, --backend uma\|orb\|mace\|aimnet2` | MLIP backend for the ML region (default: `uma`). |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction (default: off). |
| `--hessian-calc-mode Analytical\|FiniteDifference` | ML Hessian calculation mode. `Analytical` is available for the UMA backend only (recommended when VRAM is available); other backends use `FiniteDifference`. |

For a full matrix of options and YAML schemas, see [YAML Reference](yaml-reference.md).

---

## Run summaries

Every `mlmm all` run writes:

- `summary.log` -- formatted summary for quick inspection, and
- `summary.yaml` -- YAML version summary.

They typically contain:

- the exact CLI command invoked,
- global MEP statistics (e.g. maximum barrier, path length),
- per-segment barrier heights and key bond changes,
- energies from the MLIP backend, thermochemistry, and DFT post-processing (where enabled).

Each segment directory under `path_search/` also gets its own `summary.log` and `summary.yaml`, so you can inspect local refinements independently.

---

## CLI commands

Most users will primarily call `mlmm all`. The CLI also exposes individual subcommands -- each supports `-h/--help`.
`mlmm all --help` shows core options and `mlmm all --help-advanced` shows the complete list.
`scan`, `scan2d`, `scan3d`, the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`), and selected utility commands (`mm-parm`, `define-layer`, `add-elem-info`, `trj2fig`, `energy-diagram`, `oniom-export`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).

| Subcommand | Role | Documentation |
|------------|------|---------------|
| `all` | End-to-end workflow | [all](all.md) |
| `extract` | Define ML region (QM region) | [extract](extract.md) |
| `mm-parm` | Generate Amber parm7/rst7 topology | [mm_parm](mm-parm.md) |
| `define-layer` | Assign 3-layer ML/MM partitioning | [define_layer](define-layer.md) |
| `opt` | Geometry optimization | [opt](opt.md) |
| `tsopt` | Transition state optimization | [tsopt](tsopt.md) |
| `path-opt` | MEP optimization (GSM/DMF) | [path_opt](path-opt.md) |
| `path-search` | Recursive MEP search | [path_search](path-search.md) |
| `scan` | 1D bond-length scan | [scan](scan.md) |
| `scan2d` | 2D distance scan | [scan2d](scan2d.md) |
| `scan3d` | 3D distance scan | [scan3d](scan3d.md) |
| `irc` | IRC calculation | [irc](irc.md) |
| `freq` | Vibrational analysis | [freq](freq.md) |
| `dft` | Single-point DFT | [dft](dft.md) |
| `oniom-export` | Export to Gaussian ONIOM / ORCA QM/MM (`--mode g16|orca`) | [oniom_export](oniom-export.md) |
| `oniom-import` | Import Gaussian/ORCA ONIOM input and reconstruct XYZ + layered PDB | [oniom_import](oniom-import.md) |
| `trj2fig` | Plot energy profiles | [trj2fig](trj2fig.md) |
| `energy-diagram` | Draw state energy diagram from numeric values | [energy-diagram](energy-diagram.md) |
| `add-elem-info` | Repair PDB element columns | [add_elem_info](add-elem-info.md) |
| `fix-altloc` | Remove alternate location (altLoc) indicators from PDB files | [fix_altloc](fix-altloc.md) |
| `pysis` | Run pysisyphus YAML workflow (v0.1.x compat) | [pysis](pysis.md) |

```{tip}
In `all`, `tsopt`, `freq` and `irc`, setting **`--hessian-calc-mode Analytical`** (for the ML region) is strongly recommended when you have enough VRAM. Note: `Analytical` mode is only available with the UMA backend; other backends automatically use `FiniteDifference`.
```

---

## Quick reference

**Common command patterns:**

```bash
# Basic MEP search (2+ structures)
mlmm -i R.pdb P.pdb -c 'SUBSTRATE' -l 'SUB:-1'

# Full workflow with post-processing
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft

# Single structure with staged scan
mlmm -i SINGLE.pdb -c 'LIG' -l 'LIG:-1' --scan-lists '[("RES1,100,CA","LIG,200,C1",2.0)]'

# TS-only optimization
mlmm -i TS.pdb -c 'LIG' -l 'LIG:-1' --tsopt --thermo

# Individual subcommands (after running extract + mm-parm + define-layer)
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
```

**Essential options:**

| Option | Purpose |
|--------|---------|
| `-i` | Input structure(s) |
| `-c` | Substrate definition for ML region extraction |
| `-l, --ligand-charge` | Substrate charges (e.g., `'SAM:1,GPP:-3'`) |
| `--parm` | Amber parm7 topology file (required for subcommands) |
| `--model-pdb` | ML region PDB file (required for subcommands) |
| `--tsopt` | Enable TS optimization + IRC |
| `--thermo` | Run vibrational analysis |
| `--dft` | Run single-point DFT |
| `-b, --backend` | MLIP backend (`uma`, `orb`, `mace`, `aimnet2`) |
| `--embedcharge` | Enable xTB point-charge embedding correction |
| `-o, --out-dir` | Output directory |

---

## Getting help

For any subcommand:

```bash
mlmm <subcommand> --help
mlmm <subcommand> --help-advanced
mlmm all --help-advanced
```

For `all`, `--help` is intentionally short. Use `--help-advanced` to see every option.
