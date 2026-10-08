# [mlmm-toolkit]{.p2r-wordmark} Documentation

:::{container} p2r-hero-meta
[Version: v{{ release }}]{.p2r-pill} [GitHub](https://github.com/t-0hmura/mlmm_toolkit){.p2r-meta-gh} [ChemRxiv preprint](https://doi.org/10.26434/chemrxiv-2025-jft1k){.p2r-meta-paper}
:::

:::{container} p2r-hero
<img src="./mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" class="p2r-hero-figure">

{.p2r-tagline}
**mlmm-toolkit** is a Python CLI toolkit for reaction-mechanism analysis from structures such as PDB files of enzyme complexes, using ML/MM, which combines machine-learning interatomic potentials and molecular mechanics through ONIOM.

{.p2r-lead}
New to mlmm-toolkit? Start with [Getting Started](getting-started.md).

{.p2r-cta}
[Getting Started](getting-started.md){.p2r-btn .p2r-btn-primary} [Installation](installation.md){.p2r-btn .p2r-btn-install} [Open in Google Colab](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb){.p2r-btn .p2r-btn-colab}
:::

```{toctree}
:maxdepth: 2
:caption: Introduction
:hidden:

getting-started
installation
```

```{toctree}
:maxdepth: 2
:caption: Quickstart
:hidden:

Endpoint mode <quickstart-all>
Scan-list mode <quickstart-scan>
TS-only mode <quickstart-tsopt>
```

```{toctree}
:maxdepth: 2
:caption: Guides
:hidden:

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
:caption: 導入
:hidden:

ja/getting-started
ja/installation
```

```{toctree}
:maxdepth: 2
:caption: クイックスタート
:hidden:

Endpoint モード <ja/quickstart-all>
Scan-list モード <ja/quickstart-scan>
TS-only モード <ja/quickstart-tsopt>
```

```{toctree}
:maxdepth: 2
:caption: ガイド
:hidden:

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
:caption: 参照資料
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

::::{container} p2r-cards
:::{container} p2r-card p2r-card-endpoint
**Analyze the mechanism end to end from the structures before and after the reaction**

<!-- p2r-mode-stages endpoint -->

[Quickstart: all in Endpoint mode](quickstart-all.md)
:::

:::{container} p2r-card p2r-card-scan
**Analyze the mechanism end to end from one structure**

<!-- p2r-mode-stages scan -->

[Quickstart: all in Scan-list mode](quickstart-scan.md)
:::

:::{container} p2r-card p2r-card-tsonly
**Analyze the mechanism end to end from a TS structure**

<!-- p2r-mode-stages tsonly -->

[Quickstart: TS-only mode](quickstart-tsopt.md)
:::
::::

| Goal | Page |
|------|------|
| **Choose the ML region and layers, or make a run lighter** | [Building the ML region and layers](model-setup.md) |
| **Study a mechanism, or the TS search fails** | [Tips for studying reaction mechanisms](mechanism-tips.md) |
| **Optimize the TS structure with DFT** | [Optimize the TS structure with DFT](dft-backend.md) |
| **A run failed** | [Troubleshooting](troubleshooting.md) |

## Subcommands

<!-- p2r-stage-strip -->

| Subcommand | Description |
|---------|------|
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
|-------|------|
| **Common options and input requirements** | [Common options and selectors](cli-conventions.md) |
| **Frozen atoms and distance restraints (`--freeze-atoms`, `--distance-restraint`)** | {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>` |
| **Common errors and fixes** | [Troubleshooting](troubleshooting.md) |
| **CLI command reference** | [Command Reference](reference/commands/index.md) |
| **YAML configuration options** | [YAML Reference](yaml-reference.md) · [Curated YAML subset](reference/yaml.md) |
| **MLIP backend settings** | [MLIP Backends](backends.md) |
| **Files each command writes** | [Output Directory Layout](output-layout.md) |
| **Keys of `result.json` and `summary.json`** | [JSON Output Reference](json-output.md) |
| **GPU and CPU assignment, HPC** | [Device Configuration & HPC Setup](device-hpc.md) |
| **Use the ML/MM calculator from Python** | [ML/MM calculator](mlmm-calc.md) |
| **Calling mlmm-toolkit from an AI agent (MCP)** | [MCP server](mcp_server.md) |
| **Code structure (for developers)** | [Architecture](architecture.md) |
| **Terminology** | [Glossary](glossary.md) |

## System requirements

### Hardware

- **OS:** Linux (on Windows, install it in Linux under WSL2).
- **GPU:** an NVIDIA driver compatible with the backend and PyTorch wheel. CPU execution is also supported but slower.
- **VRAM / RAM:** depends on the model, the system size, and the Hessian mode; measure the peak on a representative run.

### Software

- Python >= 3.11.
- CPU or CUDA-enabled PyTorch. Prebuilt wheels include their CUDA runtime; a local toolkit is normally needed only for source builds.
- AmberTools, which `mm-parm` uses to build the Amber topology.

See [Installation](installation.md) for setup.

## Key concepts

The three layers (ML, movable MM, frozen MM) and how ONIOM combines them are explained in [Getting Started](getting-started.md); how to choose the ML region and layers is in [Building the ML region and layers](model-setup.md).

## Agent skills

`skills/` contains guides for AI agents on the CLI commands, structure I/O, backends, workflows and outputs, and HPC use.
To install them, tell your AI agent:

> Import `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` as skills, and install mlmm-toolkit by following `mlmm-install-backends`.

If you cloned the GitHub repository, you can give the local `skills/` path instead. Then you can ask, for example:

> Read *the paper*, build a model from the PDB structure *PDB ID*, and study the mechanism of *the reaction step* with the mlmm-toolkit skills.

## Citation

```bibtex
@article{ohmura2025mlmm,
  author  = {Ohmura, Takuto and Inoue, Sei and Terada, Tohru},
  title   = {ML/MM Toolkit -- Toward Accelerated Mechanistic Investigation of Enzymatic Reactions},
  journal = {ChemRxiv}, year = {2025}, doi = {10.26434/chemrxiv-2025-jft1k}
}
```

To cite the software or a specific release, use the Zenodo record:

```bibtex
@software{ohmura2026mlmm_software,
  author       = {Ohmura, Takuto},
  title        = {mlmm-toolkit},
  year         = {2026},
  version      = {0.4.0},
  url          = {https://github.com/t-0hmura/mlmm_toolkit},
  license      = {GPL-3.0-or-later},
  doi          = {10.5281/zenodo.19197863}
}
```

## License

GNU General Public License v3 or later (GPL-3.0-or-later).

## Getting Help

```bash
# General help
mlmm --help

# Command help
mlmm <subcommand> --help

# Advanced options (internal tuning)
mlmm <subcommand> --help-advanced
```

Report problems and feature requests on [GitHub Issues](https://github.com/t-0hmura/mlmm_toolkit/issues).
