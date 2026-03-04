# **ML/MM toolkit** — Toward Accelerated Mechanistic Investigation of Enzymatic Reactions

## Overview

<img src="./docs/mlmm_toolkit_overview.png" alt="Overview of ML/MM toolkit" width="100%">

Quantum mechanics/molecular mechanics (QM/MM) methods have long been used to analyze enzymatic reaction mechanisms *in silico*. While treating the active site with QM and the remainder with MM reduces the computational cost, the inherently high cost of the QM calculation remains a major limitation. Replacing QM with Machine Learning (ML) interatomic potentials yields ML/MM approaches that further reduce computational cost while retaining accuracy.

**ML/MM toolkit** is an open-source command-line toolkit for ML/MM calculations. It streamlines the workflow necessary for enzymatic reaction mechanism analyses — energy minimization, transition-state (TS) search, and vibrational analysis — to calculate the reaction free energy (&Delta;G) and activation free energy (&Delta;G<sup>&ddagger;</sup>) for a given enzyme–substrate complex structure. A link atom boundary treatment is implemented to include amino-acid residues in the ML region, and full-system Hessians are available for TS searches and vibrational analyses. To accelerate TS searches in systems comprising around 10,000 atoms, we developed a *Partial Hessian Guided Dimer (PHG-Dimer)* method that uses the active-site Hessian to determine initial dimer orientation. We also integrated a mass-scaled flattening loop to suppress spurious imaginary modes.

A **single command** can generate a first-pass enzymatic reaction path with ML/MM accuracy:

```bash
mlmm all -i R.pdb P.pdb -c PRE --ligand-charge "PRE:-2"
```

The full workflow — **MEP search → TS optimization → IRC → thermochemistry → single-point DFT** — can be run in one command:

```bash
mlmm all -i R.pdb P.pdb -c PRE --ligand-charge "PRE:-2" \
    --tsopt --thermo --dft
```

Given **(i) two or more PDB files** (R → ... → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `mlmm_toolkit` automatically:

- extracts an **active-site pocket** around user-defined substrates,
- assigns **3-layer ONIOM regions** (ML / Movable MM / Frozen) via B-factor encoding,
- generates **MM parameters** (parm7/rst7) using AmberTools,
- explores **minimum-energy paths (MEPs)** with GSM or DMF,
- *optionally* optimizes **transition states**, runs **vibrational analysis**, **IRC**, and **single-point DFT**.

> **Expectation setting for TS search**
> - Treat single-command outputs as a strong initial guess, not guaranteed final TS validation.
> - Always validate TS candidates with frequency analysis and IRC before mechanistic interpretation.

> **Important (prerequisites):**
> - Input PDB files must already contain **hydrogen atoms**.
> - When providing multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ).
> - A **parm7 topology** file (AmberTools) is required for MM calculations; use `mlmm mm-parm` to generate one.

### Supported ML potentials

| Potential | Repository |
|-----------|------------|
| **UMA** | <https://github.com/facebookresearch/fairchem> |
| **AIMNet2** | <https://github.com/isayevlab/aimnetcentral> |

The MM layer uses **OpenMM**. Any force field capable of generating parameters (e.g., AMBER ff14SB + TIP3P + GAFF2) can be used.

> **UMA-only workflow**
> If you wish to perform chemical-reaction-mechanism analysis using **UMA alone** (without the ML/MM hybrid layer), a dedicated **UMA – Pysisyphus Interface** is available at **<https://github.com/t-0hmura/uma_pysis>**.

### Key difference from pdb2reaction

`mlmm_toolkit` is based on [`pdb2reaction`](https://github.com/t-0hmura/pdb2reaction), which uses only a machine-learning potential (UMA). `mlmm_toolkit` extends this to the **ONIOM (ML/MM) scheme**:

| | pdb2reaction | mlmm_toolkit |
|---|---|---|
| Calculator | UMA only | ONIOM: UMA + hessian_ff |
| Energy | E_ML(full) | E_MM(real) + E_ML(model) - E_MM(model) |
| Boundary | None | Link atom (H cap) + Jacobian chain rule |
| Atom layers | freeze/active (2 layers) | ML / Movable MM / Frozen (3 layers) |
| Layer encoding | None | B-factor (0/10/20) |
| Additional CLI | None | `--parm`, `--model-pdb`, `--detect-layer` |

---

## Installation

### Quick install

For CUDA 12.6:
```bash
pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu126
pip install git+https://github.com/t-0hmura/mlmm_toolkit.git
huggingface-cli login
```

For CUDA 12.8 (recommended for RTX 50 series):
```bash
pip install torch==2.7.0 --index-url https://download.pytorch.org/whl/cu128
pip install git+https://github.com/t-0hmura/mlmm_toolkit.git
huggingface-cli login
```

### Prerequisites

| Requirement | Notes |
|-------------|-------|
| **Python 3.11** | — |
| **CUDA runtime >= 12.6** | CUDA 12.8 recommended for RTX 50 series |
| Linux / WSL 2 | — |

### Full setup (conda)

```bash
conda create -n mlmm python=3.11 -y
conda activate mlmm
conda install -c conda-forge ambertools pdbfixer -y
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/mlmm_toolkit.git
plotly_get_chrome -y
huggingface-cli login
```

If you are on an HPC cluster that uses *environment modules*, load CUDA **before** installing PyTorch:

```bash
module load cuda/12.6
```

> `mlmm_toolkit` internally installs `fairchem-core` and `aimnet` from forked repositories.

### Avoid AmberTools conflicts (optional)

If `AmberTools` is loaded, unload it before installing to prevent conflicts for **ParmEd**:

```bash
module unload amber
```

### For DMF method

```bash
conda install -c conda-forge cyipopt -y
```

### Authenticate for UMA (Hugging Face)

UMA model weights are on Hugging Face Hub. You need to log in once (see <https://github.com/facebookresearch/fairchem>):

```bash
huggingface-cli login
```

---

## Preparing an Enzyme–Substrate System

1. **Build a structural model of the complex.**
   Download coordinates from Protein Data Bank. If an experimental structure is not available, use structure prediction programs such as **AlphaFold3**, docking simulation programs, or GUI software such as **PyMOL**.

2. **Generate parameter/topology and coordinate files.**
   Create `.pdb`, `.parm7`, and `.rst7` files of the complex (see the [OpenMM tutorial](https://openmm.github.io/openmm-cookbook/latest/tutorials)).
   To mimic aqueous conditions, solvate the complex, then remove water molecules beyond ~6 A.
   > Note: elemental information (columns 77–78) is omitted in PDB files generated by tleap. Use `mlmm add-elem-info` to fix this.

3. **Define the ML region.**
   Use `mlmm extract` or any molecular viewer to select residues around the substrate:

   ```bash
   mlmm extract -i complex.pdb -c PRE -r 6.0 --ligand-charge "PRE:-2" -o ml_region.pdb
   ```

   **Important:** atom order, residue names, and residue numbers must match between the *full* PDB and the *ML-region* PDB. (In PyMOL, tick **"Original atom order"** when exporting.)

---

## Quick Examples

### Full workflow (multi-structure)
```bash
mlmm all -i R.pdb P.pdb -c PRE --ligand-charge "PRE:-2" \
    --tsopt --thermo --dft
```

### Scan mode (single structure)
```bash
mlmm all -i R.pdb -c PRE --ligand-charge "PRE:-2" \
    --scan-lists "[('PRE 353 O1\'','PRE 353 C3',1.2)]"
```

### TS optimization only
```bash
mlmm tsopt -i TS_candidate_layered.pdb --parm complex.parm7 \
    -q 1 --opt-mode light
```

### Generate a starter config template
```bash
mlmm init --out mlmm_all.config.yaml
mlmm all --config mlmm_all.config.yaml --dry-run
```

### Step-by-step workflow
```bash
# 1. Generate MM parameters
mlmm mm-parm -i complex.pdb --ligand-charge "PRE:-2"

# 2. Extract active-site pocket
mlmm extract -i complex.pdb -c '353' --ligand-charge "PRE:-2" -r 6.0

# 3. Assign 3-layer ONIOM regions
mlmm define-layer -i complex.pdb --model-pdb pocket.pdb --radius-freeze 10.0

# 4. Optimize geometry
mlmm opt -i complex_layered.pdb --parm complex.parm7 -q 1 --opt-mode heavy

# 5. MEP search
mlmm path-search -i R_layered.pdb P_layered.pdb --parm complex.parm7 -q 1

# 6. TS optimization
mlmm tsopt -i hei.pdb --parm complex.parm7 -q 1 --opt-mode light

# 7. Frequency analysis
mlmm freq -i ts_optimized.pdb --parm complex.parm7 -q 1

# 8. IRC
mlmm irc -i ts_optimized.pdb --parm complex.parm7 -q 1

# 9. DFT single-point
mlmm dft -i optimized.pdb --parm complex.parm7 -q 1
```

Fully working example scripts are provided in the [`examples/`](examples/) directory.
Start with the minimal [`examples/toy_system/`](examples/toy_system/) example, then explore realistic enzyme cases in [`examples/chorismate_mutase/`](examples/chorismate_mutase/) and [`examples/methyltransferase/`](examples/methyltransferase/).
For a step-by-step walkthrough, see [`examples/tutorial.md`](examples/tutorial.md).

---

## CLI Subcommands

### Workflow

| Subcommand | Role | Documentation |
|---|---|---|
| `all` | End-to-end: extraction → MEP → TS → IRC → freq → DFT | [docs/all.md](docs/all.md) |
| `init` | Generate a starter YAML template for `mlmm all` | [docs/init.md](docs/init.md) |

### Structure Preparation

| Subcommand | Role | Documentation |
|---|---|---|
| `mm-parm` | Generate parm7/rst7 topology via AmberTools | [docs/mm_parm.md](docs/mm_parm.md) |
| `extract` | Extract active-site pocket (cluster model) | [docs/extract.md](docs/extract.md) |
| `define-layer` | Assign 3-layer B-factor encoding | [docs/define_layer.md](docs/define_layer.md) |
| `add-elem-info` | Add/repair PDB element columns (77-78) | [docs/add_elem_info.md](docs/add_elem_info.md) |
| `fix-altloc` | Resolve alternate conformations | [docs/fix_altloc.md](docs/fix_altloc.md) |

### Optimization & Path Search

| Subcommand | Role | Documentation |
|---|---|---|
| `opt` | Geometry optimization (L-BFGS or RFO) | [docs/opt.md](docs/opt.md) |
| `tsopt` | TS optimization (Dimer or RS-I-RFO) | [docs/tsopt.md](docs/tsopt.md) |
| `path-opt` | MEP optimization via GSM or DMF | [docs/path_opt.md](docs/path_opt.md) |
| `path-search` | Recursive MEP search with refinement | [docs/path_search.md](docs/path_search.md) |
| `scan` | 1D bond-length driven scan | [docs/scan.md](docs/scan.md) |
| `scan2d` | 2D distance grid scan | [docs/scan2d.md](docs/scan2d.md) |
| `scan3d` | 3D distance grid scan | [docs/scan3d.md](docs/scan3d.md) |

### Analysis

| Subcommand | Role | Documentation |
|---|---|---|
| `freq` | Vibrational frequency analysis + thermochemistry | [docs/freq.md](docs/freq.md) |
| `irc` | IRC calculation (EulerPC) | [docs/irc.md](docs/irc.md) |
| `dft` | Single-point DFT (GPU4PySCF / PySCF) | [docs/dft.md](docs/dft.md) |

### Visualization & Export

| Subcommand | Role | Documentation |
|---|---|---|
| `trj2fig` | Energy plot from XYZ trajectory | [docs/trj2fig.md](docs/trj2fig.md) |
| `energy-diagram` | Energy diagram from numeric values | [docs/energy_diagram.md](docs/energy_diagram.md) |
| `oniom-export` | Generate Gaussian ONIOM / ORCA QM/MM input | [docs/oniom_export.md](docs/oniom_export.md) |
| `oniom-import` | Import ONIOM input and reconstruct XYZ/layered PDB | [docs/oniom_import.md](docs/oniom_import.md) |

---

## 3-Layer ONIOM System

`mlmm_toolkit` uses a 3-layer system encoded in PDB B-factor columns:

| Layer | B-factor | Description | Hessian |
|---|---|---|---|
| ML | 0.0 | Active site (UMA energy/force/Hessian) | Yes |
| Movable MM | 10.0 | MM atoms allowed to move | No (by B-factor alone) |
| Frozen | 20.0 | Distant protein (coordinates fixed) | No |

Covalent bonds cut at the ML/MM boundary are capped with hydrogen *link atoms*.
The ML region can therefore include entire amino-acid side chains when necessary.

Use `mlmm define-layer` to assign layers automatically based on a model PDB (the ML region).

---

## Units

| Interface | Energy | Distance | Force | Hessian |
|---|---|---|---|---|
| **Core & ASE** | eV | A | eV A<sup>-1</sup> | eV A<sup>-2</sup> |
| **PySisyphus (CLI)** | Hartree | Bohr | Ha Bohr<sup>-1</sup> | Ha Bohr<sup>-2</sup> |

Unit conversions are handled automatically inside the ML/MM calculator.

---

## Python API

Interfaces are available for **ASE** and **PySisyphus**.

### MLMMCore (Base-Level)

```python
from mlmm_toolkit.mlmm_calc import MLMMCore

core = MLMMCore(
    input_pdb    = "complex.pdb",    # Full system PDB (protein + substrate + solvent)
    real_parm7   = "complex.parm7",  # Amber topology for the full system
    model_pdb    = "ml_region.pdb",  # ML region only (trimmed PDB)
    model_charge = -1,               # Formal charge of the ML region including link H atoms
    model_mult   = 1,                # Spin multiplicity of the ML region
)

from ase.io import read
atoms = read("structure.pdb")
results = core.compute(atoms.get_positions(), return_forces=True, return_hessian=True)

energy  = results["energy"]    # float, eV
forces  = results["forces"]    # ndarray (N, 3), eV/A
hessian = results["hessian"]   # torch.Tensor (3N, 3N), eV/A^2
```

### ASE Interface

```python
from mlmm_toolkit.mlmm_calc import MLMMASECalculator

atoms.calc = MLMMASECalculator(
    input_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1, model_mult=1,
)
energy = atoms.get_potential_energy()  # eV
```

### PySisyphus Interface

```python
from pysisyphus.io.pdb import geom_from_pdb
from pysisyphus.optimizers.LBFGS import LBFGS
from mlmm_toolkit.mlmm_calc import mlmm

geom = geom_from_pdb("structure.pdb")
geom.set_calculator(mlmm(
    input_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1, model_mult=1,
))

opt = LBFGS(geom, max_cycles=10000, thresh="gau")
opt.run()
```

> **Notes**
> - `complex.pdb` is the reference PDB used when the Amber parameters were generated, whereas `structure.pdb` can contain any starting geometry.
> - If `model_charge` is omitted, the charge is estimated with **RDKit**; set it explicitly for safety.

---

## Getting Help

```bash
mlmm --help
mlmm <subcommand> --help
mlmm <subcommand> --help-advanced
```

`--help` shows core options. `--help-advanced` shows the full option list including advanced parameters.

---

## Documentation

- [**Getting Started**](docs/getting_started.md) — Installation and first steps
- [**Concepts**](docs/concepts.md) — 3-layer system, ONIOM, link atoms
- [**CLI Command Reference (generated)**](docs/reference/commands/index.md)
- [**YAML Schema (generated)**](docs/reference/yaml.md)
- [**Troubleshooting**](docs/troubleshooting.md) — Common errors and fixes
- **Full command index**: [docs/index.md](docs/index.md)

---

## References

[1] Ohmura, T., Inoue, S., Terada, T. (2025). ML/MM toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions. ChemRxiv. https://doi.org/10.26434/chemrxiv-2025-jft1k
[2] Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. http://arxiv.org/abs/2506.23971
[3] Steinmetzer, J. et al. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). https://doi.org/10.1002/qua.26390

---

## License

**ML/MM toolkit** is distributed under the **GNU General Public License version 3 (GPL-3.0)**.

***This software is still under development. Please use it at your own risk.***
