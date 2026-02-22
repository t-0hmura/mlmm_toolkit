# mlmm_toolkit: ML/MM (ONIOM) enzyme reaction mechanism analysis

## Overview

<img src="./docs/mlmm_toolkit_overview.png" alt="Overview of ML/MM toolkit" width="100%"> 

`mlmm_toolkit` is a Python CLI toolkit for **enzymatic reaction mechanism analysis** using the **ML/MM (ONIOM) hybrid method**. It combines a machine-learning interatomic potential (UMA) for the active site with a classical force field (hessian_ff) for the surrounding protein environment.

A **single command** can generate a first-pass enzymatic reaction path with ML/MM accuracy:

```bash
mlmm all -i R.pdb P.pdb -c PRE --ligand-charge "PRE:-2"
```

---

The full workflow — **MEP search → TS optimization → IRC → thermochemistry → single-point DFT** — can be run in one command:

```bash
mlmm all -i R.pdb P.pdb -c PRE --ligand-charge "PRE:-2" \
    --tsopt --thermo --dft
```

---

Given **(i) two or more PDB files** (R → ... → P), **or (ii) one PDB with `--scan-lists`**, **or (iii) one TS candidate with `--tsopt`**, `mlmm_toolkit` automatically:

- extracts an **active-site pocket** around user-defined substrates,
- assigns **3-layer ONIOM regions** (ML / Movable MM / Frozen) via B-factor encoding,
- generates **MM parameters** (parm7/rst7) using AmberTools,
- explores **minimum-energy paths (MEPs)** with GSM or DMF,
- *optionally* optimizes **transition states**, runs **vibrational analysis**, **IRC**, and **single-point DFT**,

using Meta's **UMA** machine-learning interatomic potential for the ML region and **hessian_ff** for the MM region.

> **Expectation setting for TS search**
> - Treat single-command outputs as a strong initial guess, not guaranteed final TS validation.
> - Always validate TS candidates with frequency analysis and IRC before mechanistic interpretation.

> **Important (prerequisites):**
> - Input PDB files must already contain **hydrogen atoms**.
> - When providing multiple PDBs, they must contain **the same atoms in the same order** (only coordinates may differ).
> - Boolean CLI options are toggle-style: use `--flag` / `--no-flag` (e.g., `--tsopt`, `--no-dft`).
> - A **parm7 topology** file (AmberTools) is required for MM calculations; use `mlmm mm-parm` to generate one.

## Key difference from pdb2reaction

`mlmm_toolkit` is based on [`pdb2reaction`](https://github.com/t-0hmura/pdb2reaction), which uses only a machine-learning potential (UMA). `mlmm_toolkit` extends this to the **ONIOM (ML/MM) scheme**:

| | pdb2reaction | mlmm_toolkit |
|---|---|---|
| Calculator | UMA only | ONIOM: UMA + hessian_ff |
| Energy | E_ML(full) | E_MM(real) + E_ML(model) - E_MM(model) |
| Boundary | None | Link atom (H cap) + Jacobian chain rule |
| Atom layers | freeze/active (2 layers) | ML / Movable MM / Frozen (3 layers) |
| Layer encoding | None | B-factor (0/10/20) |
| Additional CLI | None | `--real-parm7`, `--model-pdb`, `--detect-layer` |

## Documentation

- [**Getting Started**](docs/getting_started.md) — Installation and first steps
- [**Concepts**](docs/concepts.md) — 3-layer system, ONIOM, link atoms
- [**CLI Command Reference (generated)**](docs/reference/commands/index.md)
- [**YAML Schema (generated)**](docs/reference/yaml.md)
- [**YAML Reference**](docs/yaml_reference.md) — Advanced configuration (`--args-yaml`, or layered `--config`/`--override-yaml` in `all`)
- [**Troubleshooting**](docs/troubleshooting.md) — Common errors and fixes
- **Full command index**: [docs/index.md](docs/index.md)

---

## Installation

`mlmm_toolkit` requires Linux with a CUDA-capable GPU.

### Prerequisites

- Python >= 3.11
- CUDA 12.x
- AmberTools (for MM parameter generation)

### Minimal setup (CUDA 12.9, torch 2.8.0)

```bash
conda create -n mlmm python=3.11 -y
conda activate mlmm
conda install -c conda-forge ambertools pdbfixer -y
pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cu129
pip install git+https://github.com/t-0hmura/mlmm_toolkit.git
plotly_get_chrome -y
huggingface-cli login
```

### For DMF method

```bash
conda install -c conda-forge cyipopt -y
```

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
mlmm tsopt -i TS_candidate_layered.pdb --real-parm7 complex.parm7 \
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
mlmm define-layer -i complex.pdb --model-pdb pocket.pdb

# 4. Optimize geometry
mlmm opt -i complex_layered.pdb --real-parm7 complex.parm7 -q 1 --opt-mode heavy

# 5. MEP search
mlmm path-search -i R_layered.pdb P_layered.pdb --real-parm7 complex.parm7 -q 1

# 6. TS optimization
mlmm tsopt -i hei.pdb --real-parm7 complex.parm7 -q 1 --opt-mode light

# 7. Frequency analysis
mlmm freq -i ts_optimized.pdb --real-parm7 complex.parm7 -q 1

# 8. IRC
mlmm irc -i ts_optimized.pdb --real-parm7 complex.parm7 -q 1

# 9. DFT single-point
mlmm dft -i optimized.pdb --real-parm7 complex.parm7 -q 1
```

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
| `oniom-gaussian` | Generate Gaussian ONIOM input | [docs/oniom_export.md](docs/oniom_export.md) |
| `oniom-orca` | Generate ORCA ONIOM input | [docs/oniom_export.md](docs/oniom_export.md) |

---

## 3-Layer ONIOM System

`mlmm_toolkit` uses a 3-layer system encoded in PDB B-factor columns:

| Layer | B-factor | Description | Hessian |
|---|---|---|---|
| ML | 0.0 | Active site (UMA energy/force/Hessian) | Yes |
| Movable MM | 10.0 | MM atoms allowed to move | No (by B-factor alone) |
| Frozen | 20.0 | Distant protein (coordinates fixed) | No |

Hessian-target MM atoms are selected by calculator settings (`hess_cutoff`, explicit `hess_mm_atoms`, or equivalent YAML), not by a dedicated B-factor layer.

Use `mlmm define-layer` to assign layers automatically based on a model PDB (the ML region).

---

## Getting Help

```bash
mlmm --help
mlmm <subcommand> --help
mlmm <subcommand> --help-advanced
mlmm all --help-advanced
# Equivalent module invocation
python -m mlmm_toolkit --help
```

`mlmm all --help` shows core options. Use `mlmm all --help-advanced` for the full option list.
`scan`, `scan2d`, `scan3d`, the calculation commands (`opt`, `path-opt`, `path-search`, `tsopt`, `freq`, `irc`, `dft`), and selected utility commands (`mm-parm`, `define-layer`, `add-elem-info`, `trj2fig`, `energy-diagram`, `oniom-gaussian`, `oniom-orca`) now follow the same progressive-help pattern (`--help` core, `--help-advanced` full). `extract` and `fix-altloc` also support progressive help (`--help` core, `--help-advanced` full parser options).
For advanced configuration, use `--args-yaml` on individual subcommands, or `--config`/`--override-yaml` on `mlmm all`. See [YAML Reference](docs/yaml_reference.md).

---

## Docs / Smoke Checks

To keep docs and CLI behavior in sync:

```bash
python scripts/check_intro_template.py
python scripts/check_markdown_links.py
python scripts/smoke_docs_commands.py
python scripts/check_all_scan_contract.py
```

To verify trajectory dump behavior across `opt` / `tsopt` routes:

```bash
# Optional: override per-case timeout in seconds (default: 300)
export MLMM_DUMP_CASE_TIMEOUT_SEC=300
python scripts/smoke_dump_trajectories.py
```

`smoke_dump_trajectories.py` auto-detects a toy fixture (or uses `MLMM_DUMP_FIXTURE_DIR`) and checks both dump-on and dump-off paths, including legacy bool compatibility (`--dump True`).

---

## Units

| Interface | Energy | Distance | Force | Hessian |
|---|---|---|---|---|
| **Core & ASE** | eV | A | eV/A | eV/A^2 |
| **PySisyphus (CLI)** | Hartree | Bohr | Ha/Bohr | Ha/Bohr^2 |

Unit conversions are handled automatically inside the ML/MM calculator.

---

## Python API

### MLMMCore (Base-Level)

```python
from mlmm import MLMMCore

core = MLMMCore(
    real_pdb     = "complex.pdb",
    real_parm7   = "complex.parm7",
    real_rst7    = "complex.rst7",
    model_pdb    = "ml_region.pdb",
    model_charge = -1,
    model_mult   = 1,
    ml_device    = "auto",
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
from mlmm import mlmm_ase

calc = mlmm_ase(
    real_pdb="complex.pdb", real_parm7="complex.parm7",
    model_pdb="ml_region.pdb", model_charge=-1,
)
atoms.calc = calc
energy = atoms.get_potential_energy()  # eV
```

---

## References

[1] Ohmura, T., Inoue, S., Terada, T. (2025). ML/MM toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions. ChemRxiv. https://doi.org/10.26434/chemrxiv-2025-jft1k
[2] Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. http://arxiv.org/abs/2506.23971
[3] Steinmetzer, J. et al. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). https://doi.org/10.1002/qua.26390

---

## License

`mlmm_toolkit` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.
