# mlmm-toolkit Documentation

*Version: v0.2.8*

**mlmm-toolkit** is a Python CLI toolkit for automated enzymatic reaction-path modeling using **ML/MM** (machine learning / molecular mechanics) methods.

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
concepts
quickstart-all
quickstart-scan-spec
quickstart-tsopt-freq
recipes-common-errors
troubleshooting
cli-conventions
ja/getting-started
ja/concepts
ja/quickstart-all
ja/quickstart-scan-spec
ja/quickstart-tsopt-freq
ja/recipes-common-errors
ja/troubleshooting
ja/cli-conventions
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
extract
add-elem-info
mm-parm
define-layer
opt
tsopt
path-opt
path-search
scan
scan2d
scan3d
freq
irc
dft
trj2fig
oniom-export
oniom-import
fix-altloc
energy-diagram
bond-summary
device-hpc
oniom-gaussian
oniom-orca
ja/all
ja/extract
ja/add-elem-info
ja/mm-parm
ja/define-layer
ja/opt
ja/tsopt
ja/path-opt
ja/path-search
ja/scan
ja/scan2d
ja/scan3d
ja/freq
ja/irc
ja/dft
ja/trj2fig
ja/oniom-export
ja/oniom-import
ja/fix-altloc
ja/energy-diagram
ja/bond-summary
ja/device-hpc
ja/oniom-gaussian
ja/oniom-orca
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

reference/commands/index
reference/yaml
yaml-reference
json-output
mlmm-calc
python-api
pysis
glossary
ja/yaml-reference
ja/json-output
ja/mlmm-calc
ja/python-api
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

## Where to start

| Goal | Page |
|---|---|
| Install and run a first end-to-end pipeline | [Getting Started](getting-started.md) |
| Concepts (3-layer ONIOM, microiteration, link atoms) | [Concepts & Workflow](concepts.md) |
| End-to-end pipeline from a PDB (`mlmm all`) | [Quickstart: all](quickstart-all.md) |
| Single-structure staged scan (`mlmm scan`) | [Quickstart: scan](quickstart-scan-spec.md) |
| TS validation (`tsopt` + frequencies) | [Quickstart: tsopt](quickstart-tsopt-freq.md) |
| CLI conventions and input requirements | [CLI Conventions](cli-conventions.md) |
| Symptom-first failure routing | [Common Error Recipes](recipes-common-errors.md) |
| Common errors and fixes | [Troubleshooting](troubleshooting.md) |
| Abbreviations and terminology | [Glossary](glossary.md) |

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
| [`add-elem-info`](add-elem-info.md) | Repair PDB element columns (77-78) |
| [`mm-parm`](mm-parm.md) | Build AMBER topology (parm7/rst7) with tleap + GAFF2 |
| [`define-layer`](define-layer.md) | Define 3-layer ML/MM regions via B-factor annotation |

### Geometry Optimization
| Subcommand | Description |
|------------|-------------|
| [`opt`](opt.md) | Single-structure geometry optimization (L-BFGS / RFO) |
| [`tsopt`](tsopt.md) | Transition state optimization (Dimer / RS-I-RFO) |

### Path Search & Optimization
| Subcommand | Description |
|------------|-------------|
| [`path-opt`](path-opt.md) | MEP optimization via GSM or DMF (two structures) |
| [`path-search`](path-search.md) | Recursive MEP search with automatic refinement (2+ structures) |

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
| [`energy-diagram`](energy-diagram.md) | Build an energy diagram from numeric input values |
| [`bond-summary`](bond-summary.md) | Compare structures and report bond changes |

### Utilities
| Subcommand | Description |
|------------|-------------|
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate conformations (altloc) |
| [`device-hpc`](device-hpc.md) | Check GPU device information on HPC |

### Export
| Subcommand | Description |
|------------|-------------|
| [`oniom-export`](oniom-export.md) | Export to Gaussian ONIOM / ORCA QM/MM (`--mode g16` or `orca`) |
| [`oniom-import`](oniom-import.md) | Import Gaussian/ORCA ONIOM input and reconstruct XYZ + layered PDB |

---

## Configuration & Reference

| Topic | Page |
|-------|------|
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML schema** | [YAML Schema](reference/yaml.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) |
| **ML/MM calculator architecture** | [ML/MM Calculator](mlmm-calc.md) |
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
See the [YAML Reference](yaml-reference.md) for detailed options.

---

## Output Structure

Typical `mlmm all` output:
```
result_all/
├── ml_region.pdb # ML-region definition
├── summary.log # Human-readable summary
├── summary.json # Machine-readable summary
├── pockets/ # ML region structures determined by extract
├── mm_parm/ # AMBER topology files
├── scan/ # (Optional) scan results
├── path_search/ # MEP trajectories and diagrams (--no-refine-path uses path_opt/ instead)
│ ├── mep_trj.xyz # MEP trajectory
│ ├── mep.pdb # MEP in PDB format
│ └── seg_*/ # Per-segment details
└── path_search/post_seg_*/ # Post-processing outputs (--no-refine-path uses path_opt/post_seg_*/)
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

## Agent Skills

`mlmm-toolkit` ships AI-agent instructions under [`.claude/skills/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/.claude/skills) so your agent can drive ML/MM ONIOM mechanism studies via Claude Code, Cursor, etc. The bundle covers design overview, the 22 CLI subcommands with canonical recipes, structure I/O (PDB B-factor encoding, XYZ / GJF / Amber parm7+rst7), backend installation (UMA / Orb / MACE / AIMNet2 / AmberTools / DFT / xTB), workflows and `summary.json` parsing, and HPC operation (PBS / SLURM, dynamic dispatch). Copy `.claude/skills/` into your project repository or home directory.

