# mlmm-toolkit Documentation

[GitHub](https://github.com/t-0hmura/mlmm_toolkit) · [ChemRxiv preprint](https://doi.org/10.26434/chemrxiv-2025-jft1k) · [Open in Google Colab](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb)

*Version: v{{ release }}*

---

<img src="./mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

**mlmm-toolkit** is a Python CLI for modeling enzymatic reaction paths from PDB structures using ML/MM, combining machine-learning interatomic potentials and molecular mechanics through ONIOM.

New to mlmm-toolkit? Start with [Getting Started](getting-started.md).

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

getting-started
installation
quickstart-all
quickstart-scan
quickstart-tsopt
model-setup
mechanism-tips
dft-backend
troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: Commands
:hidden:

all
fix-altloc
add-elem-info
mm-parm
extract
define-layer
opt
scan
scan2d
scan3d
path-opt
path-search
tsopt
irc
freq
dft
sp
bond-summary
trj2fig
energy-diagram
oniom-export
oniom-import
```

```{toctree}
:maxdepth: 2
:caption: Reference
:hidden:

cli-conventions
reference/commands/index
reference/yaml
yaml-reference
json-output
output-layout
backends
mlmm-calc
device-hpc
mcp_server
glossary
architecture
```

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

ja/index
ja/getting-started
ja/installation
ja/quickstart-all
ja/quickstart-scan
ja/quickstart-tsopt
ja/model-setup
ja/mechanism-tips
ja/dft-backend
ja/troubleshooting
```

```{toctree}
:maxdepth: 2
:caption: コマンド
:hidden:

ja/all
ja/fix-altloc
ja/add-elem-info
ja/mm-parm
ja/extract
ja/define-layer
ja/opt
ja/scan
ja/scan2d
ja/scan3d
ja/path-opt
ja/path-search
ja/tsopt
ja/irc
ja/freq
ja/dft
ja/sp
ja/bond-summary
ja/trj2fig
ja/energy-diagram
ja/oniom-export
ja/oniom-import
```

```{toctree}
:maxdepth: 2
:caption: リファレンス
:hidden:

ja/cli-conventions
ja/yaml-reference
ja/json-output
ja/output-layout
ja/backends
ja/mlmm-calc
ja/device-hpc
ja/mcp_server
ja/glossary
ja/architecture
```

## Quick start

| Goal | Guide |
|---|---|
| Run the whole pathway from R and P | [Quickstart: all](quickstart-all.md) |
| Start from one structure (no product structure) | [Quickstart: scan](quickstart-scan.md) |
| Optimize and check a TS candidate | [Quickstart: TS-only mode](quickstart-tsopt.md) |
| Choose the ML region and layers, or make a run lighter | [Building the ML region and layers](model-setup.md) |
| Study a mechanism, or the TS search fails | [Tips for studying reaction mechanisms](mechanism-tips.md) |
| Check the TS with DFT | [Refine an MLIP TS with DFT](dft-backend.md) |
| A run failed | [Troubleshooting](troubleshooting.md) |

See [Installation](installation.md) for prerequisites.

## CLI subcommands

| Subcommand | Description |
|---|---|
| [`all`](all.md) | ML/MM model setup and MEP search; optional TS, IRC, thermochemistry, and DFT |
| [`fix-altloc`](fix-altloc.md) | Resolve PDB alternate conformations |
| [`add-elem-info`](add-elem-info.md) | Fill PDB element columns 77–78 |
| [`mm-parm`](mm-parm.md) | Build Amber parm7/rst7 topology and coordinates |
| [`extract`](extract.md) | Define the ML region from a protein–ligand complex |
| [`define-layer`](define-layer.md) | Assign ML / movable-MM / frozen-MM B-factor layers |
| [`opt`](opt.md) | Optimize a geometry with L-BFGS or RFO |
| [`scan`](scan.md) | Restrained distance scans; concerted coordinates and sequential stages |
| [`scan2d`](scan2d.md) | Two-dimensional energy landscapes |
| [`scan3d`](scan3d.md) | Three-dimensional energy landscapes |
| [`path-opt`](path-opt.md) | Optimize a two-endpoint MEP with GSM or DMF |
| [`path-search`](path-search.md) | Search and recursively refine an MEP |
| [`tsopt`](tsopt.md) | Optimize a TS candidate with RS-P-RFO, Dimer, or another supported TS optimizer |
| [`irc`](irc.md) | Trace the intrinsic reaction coordinate |
| [`freq`](freq.md) | Vibrational analysis and thermochemistry |
| [`dft`](dft.md) | Single-point DFT with GPU4PySCF or PySCF |
| [`sp`](sp.md) | ML/MM ONIOM energy and forces; optional Hessian |
| [`bond-summary`](bond-summary.md) | Report covalent bond changes between structures |
| [`trj2fig`](trj2fig.md) | Plot an XYZ trajectory's energy profile |
| [`energy-diagram`](energy-diagram.md) | Draw a state-energy diagram from numeric values |
| [`oniom-export`](oniom-export.md) | Generate Gaussian ONIOM or ORCA QM/MM input |
| [`oniom-import`](oniom-import.md) | Read an ONIOM input into XYZ / layered PDB |

## Configuration and reference

| Topic | Page |
|---|---|
| CLI conventions and input formats | [Common options and selectors](cli-conventions.md) |
| Frozen atoms and distance restraints (`--freeze-atoms`, `--distance-restraint`) | {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>` |
| Terminology | [Glossary](glossary.md) |
| YAML options | [YAML reference](yaml-reference.md) |
| Output files and JSON | [Output layout](output-layout.md) · [JSON schema](json-output.md) |
| Backends and reproducibility | [Backends](backends.md) |
| Devices and HPC | [Device and HPC setup](device-hpc.md) |
| Python API and architecture | [ML/MM calculator](mlmm-calc.md) · [Architecture](architecture.md) |
| MCP server | [MCP server](mcp_server.md) |
| Troubleshooting | [Troubleshooting](troubleshooting.md) |
| Generated CLI reference (English) | [Command reference](reference/commands/index.md) |
| Starter configuration (English) | [Curated YAML subset](reference/yaml.md) |

## System requirements

See [Installation](installation.md) for installation and backend compatibility.
The GPU and driver must meet the selected backend's requirements.
Estimate VRAM, RAM, and runtime with a representative calculation.
`mm-parm` requires AmberTools.

## Key concepts

The three layers (ML, movable MM, frozen MM) and how ONIOM combines them are explained in [Getting Started](getting-started.md#overview); how to choose the ML region and layers is in [Building the ML region and layers](model-setup.md).

## Agent skills

`skills/` contains guides for AI agents on the CLI commands, structure I/O, backends, workflows and outputs, and HPC use.
See the [Skills index](https://github.com/t-0hmura/mlmm_toolkit/blob/main/skills/README.md) for installation and the full list.

## Citation

Ohmura, T., Inoue, S., Terada, T. (2025). *ML/MM toolkit — Toward Accelerated Mechanistic Investigation of Enzymatic Reactions.* [ChemRxiv](https://doi.org/10.26434/chemrxiv-2025-jft1k).

## License

GNU General Public License version 3 or later (GPL-3.0-or-later).

## Help

```bash
mlmm --help
mlmm <command> --help
mlmm <command> --help-advanced
```

Report problems and feature requests on [GitHub Issues](https://github.com/t-0hmura/mlmm_toolkit/issues).
