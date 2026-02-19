# mlmm_toolkit Documentation

**Version: {{ version }}**

**mlmm_toolkit** is a Python CLI toolkit for automated enzymatic reaction-path modeling using ML/MM (machine-learning / molecular mechanics) methods.

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
concepts
troubleshooting
cli-conventions
ja/getting-started
ja/concepts
ja/troubleshooting
ja/cli-conventions
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
fix_altloc
energy_diagram
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
ja/fix_altloc
ja/energy_diagram
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

yaml-reference
mlmm_calc
glossary
ja/yaml-reference
ja/mlmm_calc
ja/glossary
```

```{toctree}
:maxdepth: 1
:caption: Language
:hidden:

ja/index
```


---

## Quick Start by Goal

| Use Case | Command | Guide |
|-------------------------|---------|-------|
| Run complete reaction path search from PDB | `mlmm all` | [all.md](all.md) |
| Extract QM region from protein-ligand complex | `mlmm extract` | [extract.md](extract.md) |
| Build MM topology (parm7/rst7) | `mlmm mm-parm` | [mm_parm.md](mm_parm.md) |
| Define ML/MM layers | `mlmm define-layer` | [define_layer.md](define_layer.md) |
| Optimize a single structure | `mlmm opt` | [opt.md](opt.md) |
| Find and optimize a transition state | `mlmm tsopt` | [tsopt.md](tsopt.md) |
| Search for minimum energy path | `mlmm path-search` | [path_search.md](path_search.md) |
| Run IRC from a transition state | `mlmm irc` | [irc.md](irc.md) |
| Visualize energy profile | `mlmm trj2fig` | [trj2fig.md](trj2fig.md) |
| Export to Gaussian/ORCA ONIOM | `mlmm oniom-gaussian` / `oniom-orca` | [oniom_export.md](oniom_export.md) |
| Follow worked tutorials | -- | [Tutorial](getting-started.md) |
| Understand the big picture (concepts & terms) | -- | [Concepts & Workflow](concepts.md) |
| Resolve common errors | -- | [Troubleshooting](troubleshooting.md) |
| Look up abbreviations and terms | -- | [Glossary](glossary.md) |

---

## Documentation Guide

| Topic | Page |
|-------|------|
| **Installation & first run** | [Getting Started](getting-started.md) |
| **Key terms & workflow overview** | [Concepts & Workflow](concepts.md) |
| **Common errors & fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI conventions & input requirements** | [CLI Conventions](cli-conventions.md) |

---

## CLI Subcommands

### Main Workflow
| Subcommand | Description |
|------------|-------------|
| [`all`](all.md) | End-to-end workflow: extraction -> MM parm -> MEP -> TS optimization -> IRC -> freq -> DFT |

### Structure Preparation
| Subcommand | Description |
|------------|-------------|
| [`extract`](extract.md) | Extract active-site pocket (cluster model) from protein-ligand complex |
| [`add-elem-info`](add_elem_info.md) | Repair PDB element columns (77-78) |
| [`mm-parm`](mm_parm.md) | Build AMBER topology (parm7/rst7) with tleap + GAFF2 |
| [`define-layer`](define_layer.md) | Define 4-layer ML/MM regions via B-factor annotation |

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

### Export
| Subcommand | Description |
|------------|-------------|
| [`oniom-gaussian`](oniom_export.md) | Export to Gaussian ONIOM format |
| [`oniom-orca`](oniom_export.md) | Export to ORCA ONIOM format |

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
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
- Python 3.11
- PyTorch with CUDA support
- CUDA 12.x toolkit
- AmberTools (for `mm-parm`)

---

## Quick Examples

### Basic ML/MM MEP search
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3'
```

### Full workflow with TS optimization
```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True
```

### Single-structure scan mode
```bash
mlmm -i R.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### TS-only optimization
```bash
mlmm -i TS_candidate.pdb -c 'SAM,GPP' --ligand-charge 'SAM:1,GPP:-3' \
    --tsopt True
```

---

## Key Concepts

### ML/MM 4-Layer System
mlmm_toolkit uses a 4-layer partitioning scheme encoded via PDB B-factors:
- **ML region** (B=10.0): Treated with UMA machine-learning potential
- **Hessian-MM** (B=20.0): MM atoms included in Hessian calculation
- **Movable-MM** (B=30.0): MM atoms that move during optimization
- **Frozen** (B=40.0): Fixed MM atoms

### Charge and spin
- Use `--ligand-charge` to specify unknown residue charges: `'SAM:1,GPP:-3'`
- Use `-q/--charge` to set the ML-region total charge
- Spin multiplicity is set with `-m/--multiplicity` (default `1`)

### Boolean options
All boolean CLI options must be explicitly set to `True` or `False`:
```bash
--tsopt True --thermo True --dft False
```

### YAML configuration
Advanced settings can be provided with `--args-yaml`.
```bash
mlmm all -i R.pdb P.pdb -c 'LIG' --args-yaml config.yaml
```
See the [YAML Reference](yaml-reference.md) for all options.

---

## Output Structure

Typical `mlmm all` output:
```
result_all/
├── ml_region.pdb              # ML-region definition
├── summary.log                # Human-readable summary
├── summary.yaml               # Machine-readable summary
├── pockets/                   # Extracted cluster models
├── mm_parm/                   # AMBER topology files
├── scan/                      # (Optional) scan results
├── path_search/               # MEP trajectories and diagrams
│   ├── mep.trj               # MEP trajectory
│   ├── mep.pdb               # MEP in PDB format
│   ├── mep_w_ref.pdb         # MEP merged with full system
│   └── seg_*/                 # Per-segment details
└── path_search/post_seg_*/    # Post-processing outputs
    ├── tsopt/                 # TS optimization results
    ├── irc/                   # IRC trajectories
    ├── freq/                  # Vibrational modes
    └── dft/                   # DFT results
```

---

## License

`mlmm_toolkit` is distributed under the **GNU General Public License version 3 (GPL-3.0)** and is derived from Pysisyphus.

---

## References

1. Wood, B. M. et al. (2025). UMA: A Family of Universal Models for Atoms. [arXiv:2506.23971](http://arxiv.org/abs/2506.23971)
2. Steinmetzer, J., Kupfer, S., & Grafe, S. (2021). pysisyphus: Exploring potential energy surfaces in ground and excited states. *Int. J. Quantum Chem.*, 121(3). [DOI:10.1002/qua.26390](https://doi.org/10.1002/qua.26390)

---

## Getting Help

```bash
# General help
mlmm --help

# Command help
mlmm <subcommand> --help
```
