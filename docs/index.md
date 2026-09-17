# mlmm-toolkit Documentation

*Version: v{{ release }}*

**mlmm-toolkit** is a Python CLI for modeling enzymatic reaction paths from PDB structures using ML/MM, combining machine-learning interatomic potentials and molecular mechanics through ONIOM.

<img src="./mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
cif
concepts
quickstart-all
quickstart-scan-spec
quickstart-tsopt-freq
recipes-common-errors
troubleshooting
cli-conventions
reproducibility
ja/getting-started
ja/cif
ja/concepts
ja/quickstart-all
ja/quickstart-scan-spec
ja/quickstart-tsopt-freq
ja/recipes-common-errors
ja/troubleshooting
ja/cli-conventions
ja/reproducibility
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
sp
trj2fig
oniom-export
oniom-import
fix-altloc
energy-diagram
bond-summary
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
ja/sp
ja/trj2fig
ja/oniom-export
ja/oniom-import
ja/fix-altloc
ja/energy-diagram
ja/bond-summary
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
backends
device-hpc
architecture
output-layout
mcp_server
glossary
ja/yaml-reference
ja/json-output
ja/mlmm-calc
ja/python-api
ja/device-hpc
ja/glossary
```

```{toctree}
:maxdepth: 1
:caption: Language
:hidden:

日本語 <ja/index>
```

## Quick start

See [Getting Started](getting-started.md) for installation and input preparation.

| Goal | Guide |
|---|---|
| First end-to-end run | [Quickstart: all](quickstart-all.md) |
| Start from a single structure and bond scans | [Quickstart: scan](quickstart-scan-spec.md) |
| Validate an existing TS candidate | [Quickstart: tsopt](quickstart-tsopt-freq.md) |
| Use the interactive GPU GUI | [Open Colab](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb) |
| Diagnose a failed run | [Common error recipes](recipes-common-errors.md) |

## CLI subcommands

### Main workflow

| Subcommand | Description |
|---|---|
| [`all`](all.md) | ML/MM model setup and MEP search; optional TS, IRC, thermochemistry, and DFT |

### Structure preparation

| Subcommand | Description |
|---|---|
| [`extract`](extract.md) | Define the ML region from a protein–ligand complex |
| [`add-elem-info`](add-elem-info.md) | Fill PDB element columns 77–78 |
| [`mm-parm`](mm-parm.md) | Build Amber parm7/rst7 topology and coordinates |
| [`define-layer`](define-layer.md) | Assign ML / movable-MM / frozen-MM B-factor layers |

### Geometry optimization

| Subcommand | Description |
|---|---|
| [`opt`](opt.md) | Optimize a geometry with L-BFGS or RFO |
| [`tsopt`](tsopt.md) | Optimize a TS candidate with RS-P-RFO, Dimer, or another supported TS optimizer |

### Path search and optimization

| Subcommand | Description |
|---|---|
| [`path-opt`](path-opt.md) | Optimize a two-endpoint MEP with GSM or DMF |
| [`path-search`](path-search.md) | Search and recursively refine an MEP |

### Scans

| Subcommand | Description |
|---|---|
| [`scan`](scan.md) | Restrained distance scans; concerted coordinates and sequential stages |
| [`scan2d`](scan2d.md) | Two-dimensional energy landscapes |
| [`scan3d`](scan3d.md) | Three-dimensional energy landscapes |

### Analysis and post-processing

| Subcommand | Description |
|---|---|
| [`irc`](irc.md) | Trace the intrinsic reaction coordinate |
| [`freq`](freq.md) | Vibrational analysis and thermochemistry |
| [`dft`](dft.md) | Single-point DFT with GPU4PySCF or PySCF |
| [`sp`](sp.md) | ML/MM ONIOM energy and forces; optional Hessian |
| [`trj2fig`](trj2fig.md) | Plot an XYZ trajectory's energy profile |
| [`energy-diagram`](energy-diagram.md) | Draw a state-energy diagram from numeric values |
| [`bond-summary`](bond-summary.md) | Report covalent bond changes between structures |

### Utilities

| Subcommand | Description |
|---|---|
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate conformations |

### Export and import

| Subcommand | Description |
|---|---|
| [`oniom-export`](oniom-export.md) | Generate Gaussian ONIOM or ORCA QM/MM input |
| [`oniom-import`](oniom-import.md) | Read an ONIOM input into XYZ / layered PDB |

## Configuration and reference

| Topic | Page |
|---|---|
| CLI conventions and input formats | [CLI conventions](cli-conventions.md) · [mmCIF](cif.md) |
| Concepts and terminology | [Concepts](concepts.md) · [Glossary](glossary.md) |
| YAML options | [YAML reference](yaml-reference.md) |
| Output files and JSON | [Output layout](output-layout.md) · [JSON schema](json-output.md) |
| Backends and reproducibility | [Backends](backends.md) · [Reproducibility](reproducibility.md) |
| Devices and HPC | [Device and HPC setup](device-hpc.md) |
| Python API and architecture | [Python API](python-api.md) · [ML/MM calculator](mlmm-calc.md) · [Architecture](architecture.md) |
| MCP server | [MCP server](mcp_server.md) |
| Troubleshooting | [Troubleshooting](troubleshooting.md) |
| Generated CLI reference (English) | [Command reference](reference/commands/index.md) |
| Starter configuration (English) | [Curated YAML subset](reference/yaml.md) |

## System requirements

See [Getting Started](getting-started.md#installation) for installation and backend compatibility.
The GPU and driver must meet the selected backend's requirements.
Estimate VRAM, RAM, and runtime with a representative calculation.
`mm-parm` requires AmberTools.

## Key concepts

- **Layers:** B=0 marks ML, B=10 movable MM, and B=20 frozen MM. Frozen coordinates still contribute MM nonbonded interactions. `hess_cutoff` / `hess_mm_atoms` separately select MM atoms for the Hessian.
- **Charge and spin:** `--ligand-charge` assigns residue charges (for example `'SAM:1,GPP:-3'`); `-q/--charge` overrides the ML-region net charge; `-m/--multiplicity` sets multiplicity (default 1).
- **Boolean options:** use `--flag` / `--no-flag`, for example `--tsopt --thermo --no-dft`.
- **Configuration:** see [YAML Reference](yaml-reference.md). Preview the resolved settings without optimization:

```bash
mlmm opt -i layered.pdb --parm system.parm7 -q 0 --show-config --dry-run
```

## Output layout

In MEP mode, `all` writes `summary.log` and `summary.json`, the MEP (`mep.pdb` / `mep_trj.xyz`, plus `mep.cif` for bridged input), and `energy_diagram_MEP.png`.
Reusable preparation outputs are `ml_region.pdb`, `mm_parm/`, and `layered/`.
`segments/seg_NN/` holds R/TS/P structures and requested TS/IRC/freq/DFT results; `_work/` holds intermediate preparation, scan, and path outputs.
TS-only runs use E1/TS/E2 labels and produce no MEP.
See [all](all.md#outputs) for the full tree and [Output Layout](output-layout.md) for file conventions.

## Citation

Ohmura, T., Inoue, S., Terada, T. (2025). *ML/MM toolkit — Toward Accelerated Mechanistic Investigation of Enzymatic Reactions.* [ChemRxiv](https://doi.org/10.26434/chemrxiv-2025-jft1k).

## License

GNU General Public License version 3 or later (GPL-3.0-or-later).

## Help

```bash
mlmm --help
mlmm <command> --help
```
