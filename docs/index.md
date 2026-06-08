# mlmm-toolkit Documentation

*Version: v0.3.0* — Python CLI for ML/MM ONIOM analyses of enzymatic reactions.

<img src="./mlmm_toolkit_overview.png" alt="mlmm-toolkit workflow overview" width="90%">

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
reproducibility
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

| Goal | Page |
|---|---|
| Install + run a first end-to-end pipeline | [Getting Started](getting-started.md) |
| 3-layer ONIOM, microiteration, link atoms | [Concepts & Workflow](concepts.md) |
| End-to-end pipeline from a PDB | [Quickstart: all](quickstart-all.md) |
| Single-structure staged scan | [Quickstart: scan](quickstart-scan-spec.md) |
| TS validation (`tsopt` + freq) | [Quickstart: tsopt](quickstart-tsopt-freq.md) |
| Symptom-first failure routing | [Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) |
| CLI conventions, YAML schema | [CLI Conventions](cli-conventions.md) · [YAML Reference](yaml-reference.md) |
| Python API / ML/MM calculator architecture | [Python API](python-api.md) · [ML/MM Calculator](mlmm-calc.md) |
| Configure GPU/CPU devices and submit on HPC | [Device & HPC Setup](device-hpc.md) |
| Terminology | [Glossary](glossary.md) |

## Subcommands

| Subcommand | Role |
|---|---|
| [`all`](all.md) | End-to-end: ML/MM model setup → MEP → TS → IRC → freq → DFT |
| [`extract`](extract.md) · [`mm-parm`](mm-parm.md) · [`define-layer`](define-layer.md) · [`add-elem-info`](add-elem-info.md) · [`fix-altloc`](fix-altloc.md) | Structure preparation |
| [`opt`](opt.md) · [`tsopt`](tsopt.md) | Geometry / TS optimisation |
| [`path-opt`](path-opt.md) · [`path-search`](path-search.md) | MEP optimisation / recursive refinement |
| [`scan`](scan.md) · [`scan2d`](scan2d.md) · [`scan3d`](scan3d.md) | 1D / 2D / 3D bond-distance scans |
| [`freq`](freq.md) · [`irc`](irc.md) | Vibrational analysis + thermochemistry / IRC (EulerPC) |
| [`dft`](dft.md) · [`sp`](sp.md) | Single-point DFT / single-point ML/MM ONIOM |
| [`bond-summary`](bond-summary.md) | Bond-change report between consecutive structures |
| [`trj2fig`](trj2fig.md) · [`energy-diagram`](energy-diagram.md) | Energy plot / R→TS→P diagram |
| [`oniom-export`](oniom-export.md) · [`oniom-import`](oniom-import.md) | Gaussian / ORCA QM/MM round-trip |

## Citation

Ohmura, T., Inoue, S., Terada, T. (2025). *ML/MM toolkit — Towards Accelerated Mechanistic Investigation of Enzymatic Reactions.* ChemRxiv. <https://doi.org/10.26434/chemrxiv-2025-jft1k>

## License

GNU General Public License v3 (GPL-3.0).
