# **mlmm-toolkit**: An End-to-End ML/MM ONIOM Platform for Automated Enzymatic Reaction Mechanism Analysis

[![PyPI](https://img.shields.io/pypi/v/mlmm-toolkit.svg)](https://pypi.org/project/mlmm-toolkit/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb)

`mlmm-toolkit` is an open-source CLI for **ML/MM ONIOM** analyses of enzymatic reactions.

## Overview

<img src="https://raw.githubusercontent.com/t-0hmura/mlmm_toolkit/main/docs/mlmm_toolkit_overview.png" alt="Overview of ML/MM toolkit" width="90%">

`mlmm-toolkit` replaces the QM region of conventional QM/MM with a machine-learning interatomic potential (MLIP, default: UMA) while keeping the surrounding protein under an analytical Amber force field (`hessian_ff`), and chains **ML-region selection → MM topology/layer preparation → MEP search → TS optimization → IRC → thermochemical correction → DFT single-point** in one command.

Test a reaction mechanism in a single command:

```bash
# Multi-structure MEP (R + P endpoints → MEP, with TS optimization + thermo)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

The run writes candidate R / TS / P structures, an energy diagram, and `summary.log` / `summary.json` to `result_all/`, plus a reusable ML region (`ml_region.pdb`) and Amber topology (parm7).

### What it is for

- **Trial and error on reaction mechanisms**: screen mechanisms in the full enzyme, where QM/MM with DFT alone would take too long
- **Starting structures for QM/MM**: build R, TS, and P of the full system and write Gaussian ONIOM or ORCA QM/MM input with `oniom-export`
- **High-throughput calculations over many systems**: explore reaction pathways across substrate variants and enzyme mutants

Once MLIP/MM finds a plausible path, mlmm-toolkit can take that TS straight into DFT/MM TS optimization: the TS optimization → IRC → endpoint optimization → frequency workflow runs as GPU-accelerated DFT/MM with GPU4PySCF. See [DFT backend](docs/dft-backend.md).

## Installation

- **OS:** Linux. Native Windows is not supported because AmberTools (`tleap`) is not available there.
- **Python:** 3.11 or later, 3.12 recommended. ORB requires 3.11 or 3.12.
- **GPU:** an NVIDIA GPU.
- **AmberTools (`tleap`):** builds the topology in `mm-parm` and `all`; not needed when you reuse a matching `--parm7`.
- **C++20 compiler:** for the default `hessian_ff` MM backend; `cxx-compiler` in the conda command below installs it.
- **pdbfixer:** needed only for `mm-parm --add-h`.

Full requirements and step-by-step setup: [docs/installation.md](docs/installation.md).

```bash
# 1. New env + AmberTools + CUDA-enabled PyTorch
conda create -n mlmm-toolkit python=3.12 -y && conda activate mlmm-toolkit
conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" pdbfixer cxx-compiler -y
#    (choose the official 2.13 wheel for your driver/GPU)
pip install torch==2.13.0 --index-url https://download.pytorch.org/whl/cu130

# 2. Install mlmm-toolkit
pip install mlmm-toolkit

# 3. Authenticate Hugging Face once (only required for the default UMA backend)
#    Accept the FAIR Chemistry License v1 at https://huggingface.co/facebook/UMA, then:
hf auth login                               # interactive
# OR: export HF_TOKEN=hf_xxx && hf auth login --token "$HF_TOKEN"   # CI / HPC
```

> **Avoid AmberTools conflicts:** on clusters with a system AmberTools module loaded, run `module unload amber` before installing to prevent a ParmEd conflict with the conda-installed AmberTools.

**Optional extras** (install only what you need):

| Extra | Adds |
|---|---|
| `[orb]` / `[aimnet]` | Orb / AIMNet2 MLIP backend — *not* HF-gated |
| `[dft]` / `[dft-cuda12]` | DFT calculator and standalone command with native CUDA 13 / CUDA 12 GPU4PySCF |
| `[mcp]` | Model Context Protocol server (`mlmm-mcp`) for agent clients |
| `[pdbfixer]` | PDBFixer extra (alternative to the conda install above) |
| `[openmm]` | OpenMM low-level backend, including virtual-site water models |

The MACE backend (`-b mace`) does not install into the same environment as UMA; create a dedicated environment as described in [docs/installation.md](docs/installation.md).

CUDA module-load recipes, alternative-backend installs, DMF / `cyipopt`, Plotly Chromium, and HPC job-script templates: [docs/installation.md](docs/installation.md) and [docs/device-hpc.md](docs/device-hpc.md).

## Preparing an Enzyme-Substrate System

For most systems the only hard requirement is a **PDB with explicit hydrogens** (at the intended protonation state). `mlmm all` then builds the MM topology, selects the ML region, and runs the whole pipeline in one command — see [Quick Examples](#quick-examples) for the three input modes (multi-structure R → P, single-structure scan, TS-only). The preparation steps below are **optional**.

1. **Build a structural model of the complex.**
   Download coordinates from the Protein Data Bank. If an experimental structure is not available, use structure-prediction programs such as **AlphaFold3**, **Boltz2**, or **Chai**; docking programs; or GUI software such as **PyMOL**. Add hydrogens at the intended protonation state (or let `mm-parm --add-h --ph 7` add them), and match `-l RES:CHARGE` to the H count actually present (e.g. SAM with 23 H = `SAM:1`, 22 H = `SAM:0`). For multi-structure (R → P) runs, every PDB must share the same atoms in the same order.

2. **(Optional) Build the MM topology yourself — it is automatic by default.**
   `mlmm all` (via [`mlmm mm-parm`](docs/mm-parm.md)) generates the Amber `.parm7` / `.rst7` from the PDB automatically; unknown residues (ligands, cofactors) are parameterized with GAFF2 / AM1-BCC — pass formal charges with `-l 'RES:CHARGE'`. Build the topology by hand when it helps — a custom force field, special solvation, or a system the automatic route cannot handle — then pass it with `--parm7`. To mimic aqueous conditions, solvate the complex and remove water molecules beyond ~6 Å (see the [OpenMM cookbook](https://openmm.github.io/openmm-cookbook/latest/tutorials) / tleap).

3. **(Optional) Define the ML region yourself.**
   `mlmm all` extracts the ML region from `-c/--center` and `-r/--radius` automatically. To define it yourself instead, build an ML-region PDB — with [`mlmm extract`](docs/extract.md) or any molecular viewer — and feed it to `mlmm all` (or the per-stage subcommands) with `--model-pdb`; this skips the automatic extraction:

   ```bash
   mlmm extract -i complex.pdb -c 'SAM,GPP' -r 6.0 -l 'SAM:1,GPP:-3' -o ml_region.pdb
   ```

   **Important:** the ML-region PDB is a subset of the full-system atoms: keep their order, names, residue IDs, and chain IDs unchanged, do not add link H, and use the same selection for every state (in PyMOL, tick **"Original atom order"** when exporting). See [How to construct a reliable `model.pdb`](docs/model-setup.md#how-to-construct-a-reliable-modelpdb).

## Quick Examples

The source repository includes full-system [COMT](examples/comt/README.md) and [BezA](examples/beza/README.md) endpoint mechanisms, a methyltransferase scan, and a 122-atom ML/MM fixture; see [`examples/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples). The commands below use the BezA structures and run from the repository root.

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermochemistry)
mlmm all -i examples/beza/1.R.pdb examples/beza/3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --out-dir result_mep

# Scan mode (single structure → staged bond scan → MEP)
mlmm all -i examples/beza/1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists "[('SAM 320 CS1','GPP 321 C7',1.60)]" --tsopt --thermo --out-dir result_scan

# TS-only validation (existing TS candidate)
mlmm all -i TS_candidate_layered.pdb --parm7 complex.parm7 -q 1 --tsopt --thermo --out-dir result_tsonly
```

For Gaussian ONIOM / ORCA QM/MM input-deck export and import, use [`oniom-export`](docs/oniom-export.md) / [`oniom-import`](docs/oniom-import.md). For a walkthrough, see [docs/getting-started.md](docs/getting-started.md) and [docs/quickstart-all.md](docs/quickstart-all.md).

Each stage (`mm-parm` → `extract` → `define-layer` → `opt` → `path-opt` → `tsopt` → `irc` → `freq` → `dft`) also runs as its own subcommand; see [CLI Subcommands](#cli-subcommands) for the per-stage pages.

## Output

A run writes its deliverables to `--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.pdb` for MEP-oriented and TS-only segments
- `mep_trj.pdb` / `mep_trj.xyz` — the merged reaction path; `energy_diagram_MEP.png` — barrier diagram
- `summary.log` / `summary.json`
- Reusable inputs for follow-up runs: `ml_region.pdb` (`--model-pdb`), `mm_parm/*.parm7` (`--parm7`), `layered/` (B-factor-annotated full-system PDBs)
- Directly inspectable model systems before/after link-H insertion:
  `ml_region_without_linkH.xyz` and `ml_region_with_linkH.xyz`,
  plus matching .pdb files for PDB input

Pipeline scratch lives under `_work/` (safe to delete). Full layout and filename conventions: [docs/output-layout.md](docs/output-layout.md).

## Colab GUI workspace

**An interactive GUI workspace is available in Google Colab.** It brings full-system coordinates and topology input, ML-region setup, Mol* visualization and atom picking, controls generated from the live CLI, execution, and linked MEP/IRC/result inspection into one notebook. Choose a GPU runtime and [open the Colab GUI workspace](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb).

<img src="https://raw.githubusercontent.com/t-0hmura/mlmm_toolkit/main/docs/colab_workspace.png" alt="mlmm-toolkit Colab GUI workspace showing Mol* structure setup and ML/MM controls" width="90%">

## CLI Subcommands

| Subcommand | Role | Doc |
|---|---|---|
| `all` (default) | End-to-end: extract → MM topology/layers → MEP → TS → IRC → freq → DFT | [all](docs/all.md) |
| `mm-parm` | Generate parm7/rst7 via AmberTools | [mm-parm](docs/mm-parm.md) |
| `extract` | Extract active-site pocket | [extract](docs/extract.md) |
| `define-layer` | Assign 3-layer ML/MM B-factor encoding | [define-layer](docs/define-layer.md) |
| `fix-altloc` | Resolve PDB altlocs | [fix-altloc](docs/fix-altloc.md) |
| `add-elem-info` | Repair PDB element columns | [add-elem-info](docs/add-elem-info.md) |
| `opt` | Geometry optimization | [opt](docs/opt.md) |
| `tsopt` | TS optimization | [tsopt](docs/tsopt.md) |
| `path-opt` | MEP via GSM/DMF | [path-opt](docs/path-opt.md) |
| `path-search` | Recursive MEP refinement | [path-search](docs/path-search.md) |
| `scan` / `scan2d` / `scan3d` | 1D / 2D / 3D bond-distance scans | [scan](docs/scan.md) · [scan2d](docs/scan2d.md) · [scan3d](docs/scan3d.md) |
| `freq` | Vibrational analysis + thermo | [freq](docs/freq.md) |
| `irc` | IRC (EulerPC) | [irc](docs/irc.md) |
| `dft` | Single-point DFT | [dft](docs/dft.md) |
| `sp` | Single-point ML/MM ONIOM | [sp](docs/sp.md) |
| `bond-summary` | Compare structures, report bond changes | [bond-summary](docs/bond-summary.md) |
| `trj2fig` / `energy-diagram` | Energy plot / R→TS→P diagram | [trj2fig](docs/trj2fig.md) · [energy-diagram](docs/energy-diagram.md) |
| `oniom-export` | Gaussian ONIOM / ORCA QM/MM input-deck export | [oniom-export](docs/oniom-export.md) |
| `oniom-import` | Gaussian ONIOM / ORCA QM/MM input-deck import | [oniom-import](docs/oniom-import.md) |

How the ML/MM calculator works (ONIOM energy, link atoms, units, and the Python API): [docs/mlmm-calc.md](docs/mlmm-calc.md).

## Documentation

- [Getting Started](docs/getting-started.md) · [Installation](docs/installation.md) · [Quickstart: all](docs/quickstart-all.md) · [Building the ML region and layers](docs/model-setup.md) · [DFT backend](docs/dft-backend.md) · [Troubleshooting](docs/troubleshooting.md)
- Full site: <https://t-0hmura.github.io/mlmm_toolkit/>

## Agent Skills

`skills/` holds Agent Skills that let an AI coding agent run `mlmm-toolkit` workflows and subcommands. Copy the skill folders into `.claude/skills/` or `~/.claude/skills/` for Claude Code, or into `.agents/skills/` or `~/.agents/skills/` for Codex. The list and an example copy command are in [`skills/README.md`](skills/README.md).

## Getting Help

```bash
mlmm --help                       # top-level
mlmm <subcmd> --help              # core options
mlmm <subcmd> --help-advanced     # full option set
```

Issues: <https://github.com/t-0hmura/mlmm_toolkit/issues>.

## Related tools

| Tool | Use case |
|---|---|
| [**pdb2reaction**](https://github.com/t-0hmura/pdb2reaction) | Pure-MLIP reaction paths for **cluster models and small molecules** from PDB / XYZ / GJF. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | Lightweight **YAML-driven UMA–pysisyphus interface** for quick/exploratory reaction-mechanism studies (GS / TS / IRC / ΔG). |

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point** cost and practical region size depend on method, basis,
  hardware, memory, and system; benchmark the intended setup before production.
- **Check every TS.** When `tsopt` converges or stops on an energy plateau, it computes the Hessian of the final geometry and reports n_imag; a successful TS optimization shows one imaginary mode along the reaction. IRC and endpoint optimization then confirm that it connects R and P (reaching the right R and P counts even if IRC stops early). ORB runs in fp64 by default; validate frequencies and IRC yourself if you switch to fp32.
- **CPU-only execution** may be substantially slower than GPU depending on the
  backend and system.
- `mlmm-toolkit` bundles a GPU-optimized pysisyphus fork that is **not** compatible with upstream pysisyphus — do not install it into an environment that already has upstream pysisyphus.

## Citation

```bibtex
@article{ohmura2025mlmm,
  author = {Ohmura, Takuto and Inoue, Sei and Terada, Tohru},
  title  = {ML/MM Toolkit -- Toward Accelerated Mechanistic Investigation of Enzymatic Reactions},
  year   = {2025}, journal = {ChemRxiv}, doi = {10.26434/chemrxiv-2025-jft1k}
}
```

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License version 3 or later (GPL-3.0-or-later).
