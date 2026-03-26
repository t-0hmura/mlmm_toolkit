# mlmm-toolkit Documentation

*Version: v0.2.5*

**mlmm-toolkit** is a Python CLI toolkit for automated enzymatic reaction-path modeling using ML/MM (machine learning / molecular mechanics) methods.

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting_started
quickstart_all
quickstart_scan_spec
quickstart_tsopt_freq
concepts
recipes_common_errors
troubleshooting
cli_conventions
ja/getting_started
ja/quickstart_all
ja/quickstart_scan_spec
ja/quickstart_tsopt_freq
ja/concepts
ja/recipes_common_errors
ja/troubleshooting
ja/cli_conventions
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
extract
add_elem_info
mm_parm
define_layer
opt
tsopt
path_opt
path_search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
oniom_export
oniom_import
fix_altloc
energy_diagram
device_hpc
oniom_gaussian
oniom_orca
ja/all
ja/extract
ja/add_elem_info
ja/mm_parm
ja/define_layer
ja/opt
ja/tsopt
ja/path_opt
ja/path_search
ja/scan
ja/scan2d
ja/scan3d
ja/freq
ja/irc
ja/dft
ja/trj2fig
ja/oniom_export
ja/oniom_import
ja/fix_altloc
ja/energy_diagram
ja/device_hpc
ja/oniom_gaussian
ja/oniom_orca
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

reference/commands/index
reference/yaml
yaml_reference
mlmm_calc
python_api
pysis
glossary
ja/yaml_reference
ja/mlmm_calc
ja/python_api
ja/pysis
ja/glossary
```

```{toctree}
:maxdepth: 1
:caption: Language
:hidden:

日本語 <ja/index>
```


---

## Documentation Guide

| Topic | Page |
|-------|------|
| **Installation & first run** | [Getting Started](getting_started.md) |
| **Key terms & workflow overview** | [Concepts & Workflow](concepts.md) |
| **Symptom-first failure routing** | [Common Error Recipes](recipes_common_errors.md) |
| **Common errors & fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI conventions & input requirements** | [CLI Conventions](cli_conventions.md) |

---

## Quick Start by Goal

| Objectives | Command | Guide |
|-------------------------|---------|-------|
| First run (end-to-end) | `mlmm all` | [Quickstart: all](quickstart_all.md) |
| Single-structure staged scan (`-s`) | `mlmm scan` | [Quickstart: scan](quickstart_scan_spec.md) |
| TS validation (`tsopt` + vibrational analysis) | `mlmm tsopt` | [Quickstart: tsopt](quickstart_tsopt_freq.md) |
| Run complete reaction path search from PDB | `mlmm all` | [all.md](all.md) |
| View current configuration | `mlmm opt --show-config` | [YAML Reference](yaml_reference.md) |
| Extract QM region from protein-ligand complex | `mlmm extract` | [extract.md](extract.md) |
| Build MM topology (parm7/rst7) | `mlmm mm-parm` | [mm_parm.md](mm_parm.md) |
| Define ML/MM layers | `mlmm define-layer` | [define_layer.md](define_layer.md) |
| Optimize a single structure | `mlmm opt` | [opt.md](opt.md) |
| Find transition state candidate via MEP search | `mlmm path-search` | [path_search.md](path_search.md) |
| Optimize a transition state candidate | `mlmm tsopt` | [tsopt.md](tsopt.md) |
| Run IRC from a transition state | `mlmm irc` | [irc.md](irc.md) |
| Visualize energy profile | `mlmm trj2fig` | [trj2fig.md](trj2fig.md) |
| Export to Gaussian ONIOM / ORCA QM/MM | `mlmm oniom-export` | [oniom_export.md](oniom_export.md) |
| Rebuild XYZ/layered PDB from ONIOM input | `mlmm oniom-import` | [oniom_import.md](oniom_import.md) |
| Draw state energy diagram from numeric values | `mlmm energy-diagram` | [energy_diagram.md](energy_diagram.md) |
| Follow worked tutorials | -- | [Tutorial](getting_started.md) |
| Diagnose failures by symptom | -- | [Common Error Recipes](recipes_common_errors.md) |
| Understand the big picture (concepts & terms) | -- | [Concepts & Workflow](concepts.md) |
| Resolve common errors | -- | [Troubleshooting](troubleshooting.md) |
| Look up abbreviations and terms | -- | [Glossary](glossary.md) |

---

## CLI Subcommands

### Main Workflow
| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: ML/MM model setup -> MEP search -> TS optimization -> IRC -> freq -> DFT |

### Structure Preparation
| Subcommand | Description |
|------------|-------------|
| [`extract`](extract.md) | Define ML region (QM region) from protein-ligand complex |
| [`add-elem-info`](add_elem_info.md) | Repair PDB element columns (77-78) |
| [`mm-parm`](mm_parm.md) | Build AMBER topology (parm7/rst7) with tleap + GAFF2 |
| [`define-layer`](define_layer.md) | Define 3-layer ML/MM regions via B-factor annotation |

### Geometry Optimization
| Subcommand | Description |
|------------|-------------|
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS / RFO) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer / RS-I-RFO) |

### Path Search & Optimization
| Subcommand | Description |
|------------|-------------|
| [`path-opt`](path_opt.md) | MEP optimization via GSM or DMF (two structures) |
| [`path-search`](path_search.md) | Recursive MEP search with automatic refinement (2+ structures) |

### Scans
| Subcommand | Description |
|------------|-------------|
| [`scan`](scan.md) | 1D bond-length driven scan with restraints |
| [`scan2d`](scan2d.md) | 2D distance grid scan |
| [`scan3d`](scan3d.md) | 3D distance grid scan |

### Analysis & Post-processing
| Subcommand | Description |
|------------|-------------|
| [`irc`](irc.md) | Intrinsic Reaction Coordinate calculation |
| [`freq`](freq.md) | Vibrational frequency analysis & thermochemistry |
| [`dft`](dft.md) | Single-point DFT calculations (GPU4PySCF / PySCF) |
| [`trj2fig`](trj2fig.md) | Plot energy profiles from XYZ trajectories |
| [`energy-diagram`](energy_diagram.md) | Build an energy diagram from numeric input values |

### Export
| Subcommand | Description |
|------------|-------------|
| [`oniom-export`](oniom_export.md) | Export to Gaussian ONIOM / ORCA QM/MM (`--mode g16` or `orca`) |
| [`oniom-import`](oniom_import.md) | Import Gaussian/ORCA ONIOM input and reconstruct XYZ + layered PDB |

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML schema** | [YAML Schema](reference/yaml.md) |
| **YAML configuration options** | [YAML Reference](yaml_reference.md) |
| **ML/MM calculator architecture** | [ML/MM Calculator](mlmm_calc.md) |
| **Terminology** | [Glossary](glossary.md) |

---

## System Requirements

### Hardware
- **OS**: Linux (Ubuntu 20.04+ or CentOS 8+ tested)
- **GPU**: CUDA 12.x compatible
- **VRAM**: Minimum 8 GB (16 GB+ recommended for 1000+ atoms)
- **RAM**: 16 GB+ recommended

### Software
- Python >= 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit
- AmberTools (for `mm-parm`)

---

## Quick Examples

### Basic ML/MM MEP search
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

### Full workflow with TS optimization
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft
```

### Single-structure scan mode
```bash
mlmm scan -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -s scan.yaml --print-parsed
```

### TS-only optimization
```bash
mlmm -i TS_candidate.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
 --tsopt
```

---

## Key Concepts

### ML/MM 3-Layer System
mlmm uses a 3-layer partitioning scheme encoded via PDB B-factors:
- **ML region** (B=0.0): Treated with the selected MLIP backend (default: UMA)
- **Movable-MM** (B=10.0): MM atoms that move during optimization
- **Frozen** (B=20.0): Fixed MM atoms that provide a static potential field. Their coordinates do not change during optimization, but their non-bonded interactions (electrostatics and van der Waals) with the Movable-MM and ML regions are still included in the MM energy evaluation.

Hessian-target MM atoms are selected by calculator options (`hess_cutoff` / explicit lists), not by a dedicated B-factor layer.

### Charge and spin
- Use `--ligand-charge` to specify unknown residue charges: `'SAM:1,GPP:-3'`
- Use `-q/--charge` to set the ML-region net charge
- Spin multiplicity is set with `-m/--multiplicity` (default `1`)

### Boolean options
Boolean CLI options use toggle form (`--flag` / `--no-flag`):
```bash
--tsopt --thermo --no-dft
```

### YAML configuration
See the [YAML Reference](yaml_reference.md) for detailed options.

---

## Output Structure

Typical `mlmm all` output:
```
result_all/
├── ml_region.pdb # ML-region definition
├── summary.log # Human-readable summary
├── summary.yaml # Machine-readable summary
├── pockets/ # ML region structures determined by extract
├── mm_parm/ # AMBER topology files
├── scan/ # (Optional) scan results
├── path_search/ # MEP trajectories and diagrams
│ ├── mep_trj.xyz # MEP trajectory
│ ├── mep.pdb # MEP in PDB format
│ └── seg_*/ # Per-segment details
└── path_search/post_seg_*/ # Post-processing outputs
 ├── tsopt/ # TS optimization results
 ├── irc/ # IRC trajectories
 ├── freq/ # Vibrational modes
 └── dft/ # DFT results
```

---

## Citation

If you use this software in your research, please cite:

[1] Ohmura, T., Inoue, S., Terada, T. (2025). ML/MM toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions. ChemRxiv. https://doi.org/10.26434/chemrxiv-2025-jft1k

## License

`mlmm-toolkit` is distributed under the **GNU General Public License version 3 (GPL-3.0)**.

---

## Getting Help

```bash
# General help
mlmm --help

# Command help
mlmm <subcommand> --help
```

---

*Note: This documentation is under active development. Some sections may be incomplete or subject to change.*
