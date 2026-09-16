# **mlmm-toolkit**: An End-to-End ML/MM ONIOM Platform for Automated Enzymatic Reaction Mechanism Analysis

[![PyPI](https://img.shields.io/pypi/v/mlmm-toolkit.svg)](https://pypi.org/project/mlmm-toolkit/) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb)

## Overview

<img src="https://raw.githubusercontent.com/t-0hmura/mlmm_toolkit/main/docs/mlmm_toolkit_overview.png" alt="Overview of ML/MM toolkit" width="90%">

`mlmm-toolkit` is an open-source CLI for **ML/MM ONIOM** analyses of enzymatic reactions. It replaces the QM region of conventional QM/MM with a machine-learning interatomic potential (MLIP, default: UMA) while keeping the surrounding protein under an analytical Amber force field (`hessian_ff`), and chains **ML-region selection → MM topology/layer preparation → MEP search → TS optimization → IRC → thermochemical correction → DFT single-point** in one command. A link-atom boundary handles amino-acid residues straddling the ML/MM cut, and a microiteration scheme separates ML and MM relaxation in large systems.

Test a reaction mechanism in a single command:

```bash
# Multi-structure MEP (R + P endpoints → MEP, with TS optimization + thermo)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

The source repository includes full-system COMT and [BezA](examples/beza/README.md) endpoint mechanisms, a methyltransferase scan, and a 122-atom ML/MM fixture; see [`examples/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples). Each stage is also exposed as an [individual subcommand](#cli-subcommands).

> **Prerequisites:** input PDB/mmCIF structures must already contain hydrogens; multiple reaction states must share the same atoms in the same order (only coordinates differ). `mlmm all` runs `mm-parm` automatically. Match `-l RES:CHARGE` to the H count actually present (e.g. SAM with 23 H = `SAM:1` cation, 22 H = `SAM:0` neutral) — full input-prep checklist in [docs/getting-started.md](docs/getting-started.md).

## Colab GUI workspace

**An interactive GUI workspace is available in Google Colab.** It brings full-system coordinates and topology input, ML-region setup, Mol* visualization and atom picking, controls generated from the live CLI, execution, and linked MEP/IRC/result inspection into one notebook. Choose a GPU runtime and [open the Colab GUI workspace](https://colab.research.google.com/github/t-0hmura/mlmm_toolkit/blob/main/examples/mlmm_colab.ipynb).

<img src="https://raw.githubusercontent.com/t-0hmura/mlmm_toolkit/main/docs/colab_workspace.png" alt="mlmm-toolkit Colab GUI workspace showing Mol* structure setup and ML/MM controls" width="90%">

## Related tools

| Tool | Use case |
|---|---|
| [**pdb2reaction**](https://github.com/t-0hmura/pdb2reaction) | Pure-MLIP reaction paths for **cluster models and small molecules** from PDB / XYZ / GJF. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | Lightweight **YAML-driven UMA–pysisyphus interface** for quick/exploratory reaction-mechanism studies (GS / TS / IRC / ΔG). |

> `mlmm-toolkit` bundles a GPU-optimized pysisyphus fork that is **not** compatible with upstream pysisyphus — do not install it into an environment that already has upstream pysisyphus.

## Documentation

- [Getting Started](docs/getting-started.md) · [mmCIF and Large Structures](docs/cif.md) · [Concepts](docs/concepts.md) · [Installation](docs/getting-started.md#installation) · [Troubleshooting](docs/troubleshooting.md)
- [Python API](docs/python-api.md) · [CLI Conventions](docs/cli-conventions.md) · [YAML Reference](docs/yaml-reference.md) · [JSON Output Schema](docs/json-output.md)
- Full command index: [docs/index.md](docs/index.md)

## System requirements

| Component | Requirement |
|---|---|
| OS / Python | Linux recommended; native Windows unsupported (AmberTools/`tleap` unavailable). Python 3.11–3.12. |
| GPU / CUDA / VRAM | A backend-compatible NVIDIA GPU/driver for GPU execution; size VRAM from a representative target-system pilot. |
| RAM / Disk | Size RAM and disk for the selected backend, model cache, topology tools, and expected artifacts. |

**AmberTools** (`tleap`) is required to generate a topology with `mm-parm` or `all`; reuse an existing matching `--parm` to skip that preparation. The default `hessian_ff` backend needs **a C++20-capable compiler** (validated with GCC 13.3) to JIT-compile native kernels on first use. **pdbfixer** is optional — only `mm-parm --add-h` needs it — `conda install -c conda-forge ambertools pdbfixer -y` installs both. CPU-only ML/MM execution is supported but can be substantially slower than GPU execution; benchmark the selected backend and system. Full requirement and tuning details: [docs/getting-started.md#installation](docs/getting-started.md#installation).

## Installation

```bash
# 1. New env + AmberTools + CUDA-enabled PyTorch
conda create -n mlmm-toolkit python=3.12 -y && conda activate mlmm-toolkit
conda install -c conda-forge ambertools pdbfixer -y
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
| `[dft]` | PySCF + GPU4PySCF single-point DFT (`--dft` / `mlmm dft`); cost and memory depend on the system and method |
| `[mcp]` | Model Context Protocol server (`mlmm-mcp`) for agent clients |
| `[pdbfixer]` | PDBFixer extra (alternative to the conda install above) |
| `[openmm]` | OpenMM low-level backend, including virtual-site water models |

The MACE backend (`-b mace`) is **not** a pip extra: `mace-torch` pins `e3nn==0.4.4`, which conflicts with `fairchem-core`'s `e3nn>=0.5` (UMA), so it needs a dedicated environment — `pip uninstall -y fairchem-core && pip install mace-torch` (see [docs/getting-started.md#installation](docs/getting-started.md#installation)).

CUDA module-load recipes, alternative-backend installs, DMF / `cyipopt`, Plotly Chromium, and HPC job-script templates: [docs/getting-started.md](docs/getting-started.md#installation) and [docs/device-hpc.md](docs/device-hpc.md).

## Preparing an Enzyme-Substrate System

For most systems the only hard requirement is a **PDB with explicit hydrogens** (at the intended protonation state). `mlmm all` then builds the MM topology, selects the ML region, and runs the whole pipeline in one command — see [Quick Examples](#quick-examples) for the three input modes (multi-structure R → P, single-structure scan, TS-only). The preparation steps below are **optional**.

1. **Build a structural model of the complex.**
   Download coordinates from the Protein Data Bank. If an experimental structure is not available, use structure-prediction programs such as **AlphaFold3**, **Boltz2**, or **Chai**; docking programs; or GUI software such as **PyMOL**. Add hydrogens at the intended protonation state (or let `mm-parm --add-h --ph 7` add them). For multi-structure (R → P) runs, every PDB must share the same atoms in the same order.

2. **(Optional) Build the MM topology yourself — it is automatic by default.**
   `mlmm all` (via [`mlmm mm-parm`](docs/mm-parm.md)) generates the Amber `.parm7` / `.rst7` from the PDB automatically; unknown residues (ligands, cofactors) are parameterized with GAFF2 / AM1-BCC — pass formal charges with `-l 'RES:CHARGE'`. Build the topology by hand when it helps — a custom force field, special solvation, or a system the automatic route cannot handle — then pass it with `--parm`. To mimic aqueous conditions, solvate the complex and remove water molecules beyond ~6 Å (see the [OpenMM cookbook](https://openmm.github.io/openmm-cookbook/latest/tutorials) / tleap).
   With an explicit `--out-prefix` (or `--add-h`), `mm-parm` also exports LEaP's topology-matched PDB and fills missing element columns while preserving its atom records and order.

3. **(Optional) Define the ML region yourself.**
   `mlmm all` extracts the ML region from `-c/--center` and `-r/--radius` automatically. To define it yourself instead, build an ML-region PDB — with [`mlmm extract`](docs/extract.md) or any molecular viewer — and feed it to `mlmm all` (or the per-stage subcommands) with `--model-pdb`; this skips the automatic extraction:

   ```bash
   mlmm extract -i complex.pdb -c 'SAM,GPP' -r 6.0 -l 'SAM:1,GPP:-3' -o ml_region.pdb
   ```

   **Important:** `model.pdb` must be an unchanged subset of the full PDB/`parm7` atom topology. Preserve original atom order, names, residue IDs, and chain IDs; do not add link H manually. Terminate retained backbone fragments consistently at Cα (`CA`), put other ML/MM boundaries on aliphatic C–C single bonds whenever possible, and avoid peptide/polar/conjugated/metal bonds. Use the same selection for every R/IM/P state. (In PyMOL, tick **"Original atom order"** when exporting.) See [How to construct `model.pdb`](docs/concepts.md#how-to-construct-a-reliable-modelpdb).

   ML-region precedence is explicit: `--model-pdb` wins; otherwise `all` uses
   `-c/--center` extraction, while per-stage commands may use
   `--model-indices` or `--detect-layer` (B factors 0/10/20). You never specify
   link H for the normal case. The calculator finds every `parm7` bond crossing
   the ML selection and inserts one link H there; Cartesian distance is used
   only to place that H along the known bond, not to decide whether a bond
   exists. `--link-atom-method` selects scaled or fixed placement.

## Quick Examples

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermo + DFT)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft

# Scan mode (single structure → staged bond scans → MEP)
mlmm all -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --scan-lists "[('SAM 359 CS1','GPP 360 C8',1.3)]"

# TS-only validation (existing TS candidate)
mlmm all -i TS_candidate_layered.pdb --parm complex.parm7 -q 1 --tsopt --opt-mode grad
```

For Gaussian-ONIOM / ORCA-QM/MM input-deck export and import use [`oniom-export`](docs/oniom-export.md) / [`oniom-import`](docs/oniom-import.md). Per-stage walkthrough (`mm-parm` → `extract` → `define-layer` → `opt` → `path-search` → `tsopt` → `freq` → `irc` → `dft`): [docs/getting-started.md](docs/getting-started.md) and [docs/quickstart-all.md](docs/quickstart-all.md). Working examples (COMT, BezA, methyltransferase, and toy system): [examples/](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples).

## Output

A run writes its deliverables to `--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.pdb` for MEP-oriented segments; TS-only mode writes chemically unassigned `{e1,ts,e2}.pdb`
- `mep.pdb` / `mep_trj.xyz` — the merged reaction path; `energy_diagram_MEP.png` — barrier diagram
- `summary.log` / `summary.json`
- Reusable inputs for follow-up runs: `ml_region.pdb` (`--model-pdb`), `mm_parm/*.parm7` (`--parm`), `layered/` (B-factor-annotated full-system PDBs)
- Directly inspectable model systems before/after link-H insertion:
  `ml_region_without_linkH.xyz` and `ml_region_with_linkH.xyz`, plus matching
  PDB companions for PDB input

Pipeline scratch lives under `_work/` (safe to delete). Full layout and filename conventions: [docs/output-layout.md](docs/output-layout.md).

## CLI Subcommands

| Subcommand | Role | Doc |
|---|---|---|
| `all` (default) | End-to-end: extract → MM topology/layers → MEP → TS → IRC → freq → DFT | [all](docs/all.md) |
| `mm-parm` | Generate parm7/rst7 via AmberTools | [mm-parm](docs/mm-parm.md) |
| `extract` | Extract active-site pocket | [extract](docs/extract.md) |
| `define-layer` | Assign 3-layer ML/MM B-factor encoding | [define-layer](docs/define-layer.md) |
| `add-elem-info` / `fix-altloc` | Repair PDB element columns / resolve altlocs | [add-elem-info](docs/add-elem-info.md) · [fix-altloc](docs/fix-altloc.md) |
| `opt` / `tsopt` | Geometry / TS optimization | [opt](docs/opt.md) · [tsopt](docs/tsopt.md) |
| `path-opt` / `path-search` | MEP via GSM/DMF; recursive refinement | [path-opt](docs/path-opt.md) · [path-search](docs/path-search.md) |
| `scan` / `scan2d` / `scan3d` | 1D / 2D / 3D bond-distance scans | [scan](docs/scan.md) · [scan2d](docs/scan2d.md) · [scan3d](docs/scan3d.md) |
| `freq` / `irc` | Vibrational analysis + thermo / IRC (EulerPC) | [freq](docs/freq.md) · [irc](docs/irc.md) |
| `dft` / `sp` | Single-point DFT / single-point ML/MM ONIOM | [dft](docs/dft.md) · [sp](docs/sp.md) |
| `bond-summary` | Compare structures, report bond changes | [bond-summary](docs/bond-summary.md) |
| `trj2fig` / `energy-diagram` | Energy plot / R→TS→P diagram | [trj2fig](docs/trj2fig.md) · [energy-diagram](docs/energy-diagram.md) |
| `oniom-export` / `oniom-import` | Gaussian ONIOM / ORCA QM/MM input-deck exchange | [oniom-export](docs/oniom-export.md) · [oniom-import](docs/oniom-import.md) |

3-layer system (ML / Movable-MM / Frozen-MM, B-factor encoded), link-atom treatment, and units (energy: eV or Hartree; coordinates: Å or Bohr; forces: eV/Å or Hartree/Bohr): [docs/concepts.md](docs/concepts.md). Python API (`MLMMCore`, `MLMMASECalculator`, pysisyphus `mlmm` calculator): [docs/python-api.md](docs/python-api.md).

## Getting Help

```bash
mlmm --help                       # top-level
mlmm <subcmd> --help              # core options
mlmm <subcmd> --help-advanced     # full option set
```

Issues: <https://github.com/t-0hmura/mlmm_toolkit/issues>.

## Citation

```bibtex
@article{ohmura2025mlmm,
  author = {Ohmura, Takuto and Inoue, Sei and Terada, Tohru},
  title  = {ML/MM Toolkit -- Toward Accelerated Mechanistic Investigation of Enzymatic Reactions},
  year   = {2025}, journal = {ChemRxiv}, doi = {10.26434/chemrxiv-2025-jft1k}
}
```

## Agent Skills

Agent Skills for Claude Code / Codex / Cursor etc. in [`skills/`](skills/) — copy into your project's skill location (e.g. `.claude/skills/`) to let an agent drive `mlmm-toolkit` workflows and subcommands.

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point** cost and practical region size depend on method, basis,
  hardware, memory, and system; benchmark the intended setup before production.
- **MLIP backends** can differ in stationary-point curvature; validate the
  selected backend on representative structures and inspect the modes.
- **CPU-only execution** may be substantially slower than GPU depending on the
  backend and system; AmberTools (`tleap`) is required for `mm-parm`.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License version 3 or later (GPL-3.0-or-later).
