# **ML/MM toolkit** — Towards Accelerated Mechanistic Investigation of Enzymatic Reactions

## Overview

<img src="./docs/mlmm_toolkit_overview.png" alt="Overview of ML/MM toolkit" width="90%">

`mlmm-toolkit` is an open-source CLI for **ML/MM ONIOM** analyses of enzymatic reactions. It replaces the QM region of conventional QM/MM with a machine-learning interatomic potential (MLIP, default: UMA) while keeping the surrounding protein under an analytical Amber force field (`hessian_ff`), and chains **MM parametrization → ML-region selection → MEP search → TS optimization → IRC → frequencies → DFT single-point** in one command. A link-atom boundary handles amino-acid residues straddling the ML/MM cut, and a microiteration scheme makes TS optimization and Hessian-based methods tractable on ~10 000-atom systems.

A useful initial reaction path is one command:

```bash
# Multi-structure MEP (R + P endpoints → MEP, with TS optimisation + thermo)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

For scan-mode on a single structure and the bundled methyltransferase walk-through, see [`examples/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples). Each stage is also exposed as an [individual subcommand](#cli-subcommands).

> **Prerequisites:** input PDBs must already contain hydrogens; multiple PDBs must share the same atoms in the same order (only coordinates differ). `mlmm all` runs `mm-parm` automatically. Match `-l RES:CHARGE` to the H count actually present (e.g. SAM with 23 H = `SAM:1` cation, 22 H = `SAM:0` neutral) — full input-prep checklist in [docs/getting-started.md](docs/getting-started.md).

## Related tools

| Tool | Use case |
|---|---|
| **`mlmm-toolkit`** (this repo) | **ML/MM ONIOM** with the full protein environment; automates MM parameterization and ML-region assignment from a single PDB. |
| [**pdb2reaction**](https://github.com/t-0hmura/pdb2reaction) | Pure-MLIP reaction paths for **cluster models and small molecules** from PDB / XYZ / GJF — no MM force field required. |
| [**uma_pysis**](https://github.com/t-0hmura/uma_pysis) | Lightweight **YAML-driven UMA–pysisyphus interface** for single PES jobs (GS / TS / IRC / ΔG). |

`mlmm-toolkit` and `pdb2reaction` bundle the same GPU-optimized pysisyphus fork; it is **not** compatible with upstream pysisyphus — do not install them side by side.

## Documentation

- [Getting Started](docs/getting-started.md) · [Concepts](docs/concepts.md) · [Installation](docs/getting-started.md#installation) · [Troubleshooting](docs/troubleshooting.md)
- [Python API](docs/python-api.md) · [CLI Conventions](docs/cli-conventions.md) · [YAML Reference](docs/yaml-reference.md) · [JSON Output Schema](docs/json-output.md)
- Full command index: [docs/index.md](docs/index.md)

## System requirements

| Component | Requirement |
|---|---|
| OS / Python | Linux x86_64 (validated); WSL 2 also works. macOS untested; native Windows unsupported (AmberTools/`tleap` unavailable). Python >= 3.11 (3.12 tested). |
| GPU / CUDA / VRAM | NVIDIA GPU, CUDA >= 12.6 (12.9 recommended; required for RTX 50-series, matched to the PyTorch wheel). 8 GB VRAM minimum, 16 GB recommended (24 GB for analytical Hessian / `mm_backend: openmm`). |
| RAM / Disk | 32 GB RAM minimum (120 GB recommended for enzyme active-site models); 20 GB free disk for the conda env, AmberTools, UMA cache, and artifacts. |

Required external tools: **AmberTools** (`tleap`) and **pdbfixer** — `conda install -c conda-forge ambertools pdbfixer -y`. CPU-only execution works for setup commands but is 10–100× slower for any ML/MM dynamics or Hessian step. Full requirement and tuning details: [docs/getting-started.md#installation](docs/getting-started.md#installation).

## Installation

```bash
# 1. New env + AmberTools + CUDA-enabled PyTorch
conda create -n mlmm-toolkit python=3.11 -y && conda activate mlmm-toolkit
conda install -c conda-forge ambertools pdbfixer -y
pip install torch --index-url https://download.pytorch.org/whl/cu129

# 2. mlmm-toolkit (editable from a local clone, or `pip install mlmm-toolkit`)
pip install -e .

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
| `[dft]` | PySCF + GPU4PySCF single-point DFT (`--dft` / `mlmm dft`) — practical to ~500 ML-region atoms |
| `[mcp]` | Model Context Protocol server (`mlmm-mcp`) for agent clients |
| `[pdbfixer]` | PDBFixer extra (alternative to the conda install above) |

The MACE backend (`-b mace`) is **not** a pip extra: `mace-torch` pins `e3nn==0.4.4`, which conflicts with `fairchem-core`'s `e3nn>=0.5` (UMA), so it needs a dedicated environment — `pip uninstall -y fairchem-core && pip install mace-torch` (see [docs/getting-started.md#installation](docs/getting-started.md#installation)).

CUDA module-load recipes, alternative-backend installs, DMF / `cyipopt`, Plotly Chromium, and HPC job-script templates: [docs/getting-started.md](docs/getting-started.md#installation) and [docs/device-hpc.md](docs/device-hpc.md).

## Preparing an Enzyme-Substrate System

1. **Build a structural model of the complex.**
   Download coordinates from Protein Data Bank. If an experimental structure is not available, use structure prediction programs such as **AlphaFold3**, docking simulation programs, or GUI software such as **PyMOL**.

2. **Generate parameter/topology and coordinate files.**
   Create `.pdb`, `.parm7`, and `.rst7` files of the complex (see the [OpenMM tutorial](https://openmm.github.io/openmm-cookbook/latest/tutorials)).
   To mimic aqueous conditions, solvate the complex, then remove water molecules beyond ~6 Å.
   > Note: elemental information (columns 77–78) is omitted in PDB files generated by tleap. Use [`mlmm add-elem-info`](docs/add-elem-info.md) to fix this.

3. **Define the ML region.**
   Use [`mlmm extract`](docs/extract.md) or any molecular viewer to select residues around the substrate:

   ```bash
   mlmm extract -i complex.pdb -c 'SAM,GPP' -r 6.0 -l 'SAM:1,GPP:-3' -o ml_region.pdb
   ```

   **Important:** atom order, residue names, and residue numbers must match between the *full* PDB and the *ML-region* PDB. (In PyMOL, tick **"Original atom order"** when exporting.)

## Quick Examples

```bash
# Multi-structure MEP (R + P → MEP, with TS + thermo + DFT)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft

# Scan mode (single structure → staged bond scans → MEP)
mlmm all -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --scan-lists "[('SAM 359 CS1','GPP 360 C8',1.3)]"

# TS-only validation (existing TS candidate)
mlmm tsopt -i TS_candidate_layered.pdb --parm complex.parm7 -q 1 --opt-mode grad
```

For Gaussian-ONIOM / ORCA-QM/MM round-trips use [`oniom-export`](docs/oniom-export.md) / [`oniom-import`](docs/oniom-import.md). Per-stage walkthrough (`mm-parm` → `extract` → `define-layer` → `opt` → `path-search` → `tsopt` → `freq` → `irc` → `dft`): [docs/getting-started.md](docs/getting-started.md) and [docs/quickstart-all.md](docs/quickstart-all.md). Working scripts (methyltransferase + toy_system): [examples/](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples).

## Output

A run writes its deliverables to `--out-dir` (default `./result_all/`):

- `segments/seg_NN/{reactant,ts,product}.pdb` — the canonical R / TS / P structures to cite
- `mep.pdb` / `mep_trj.xyz` — the merged reaction path; `energy_diagram_MEP.png` — barrier diagram
- `summary.log` (human-readable) / `summary.json` (machine-readable)
- Reusable inputs for follow-up runs: `ml_region.pdb` (`--model-pdb`), `mm_parm/*.parm7` (`--parm`), `layered/` (B-factor-annotated full-system PDBs)

Pipeline scratch lives under `_work/` (safe to delete). Full layout and filename conventions: [docs/output-layout.md](docs/output-layout.md).

## CLI Subcommands

| Subcommand | Role | Doc |
|---|---|---|
| `all` (default) | End-to-end: `mm-parm` → extract → MEP → TS → IRC → freq → DFT | [all](docs/all.md) |
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
| `oniom-export` / `oniom-import` | Gaussian ONIOM / ORCA QM/MM round-trip | [oniom-export](docs/oniom-export.md) · [oniom-import](docs/oniom-import.md) |

3-layer system (ML / Movable-MM / Frozen-MM, B-factor encoded), link-atom treatment, units (eV·Å in core / Ha·Bohr in pysisyphus CLI): [docs/concepts.md](docs/concepts.md). Python API (`MLMMCore`, `MLMMASECalculator`, pysisyphus `mlmm` calculator): [docs/python-api.md](docs/python-api.md).

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
  title  = {ML/MM Toolkit -- Towards Accelerated Mechanistic Investigation of Enzymatic Reactions},
  year   = {2025}, journal = {ChemRxiv}, doi = {10.26434/chemrxiv-2025-jft1k}
}
```

## Agent Skills

Agent instructions for Claude Code / Cursor live in [`skills/`](skills/) — copy into your project's skill location (e.g. `.claude/skills/`) to let an agent drive `mm-parm` / `extract` / `tsopt` / `irc` / `dft` end-to-end.

## Known limitations

- **MACE + UMA cannot coexist** (`e3nn` version conflict). Use separate conda envs.
- **DFT single-point** is practical to ~500 ML-region atoms; larger regions need fragmentation or an external QM program.
- **ORB backend** sometimes converges TS with extra soft imaginary modes — prefer UMA / MACE or re-score with DFT for clean single-saddle spectra.
- **CPU-only execution** is 10–100× slower than GPU; AmberTools (`tleap`) is required for `mm-parm`.

## Contributing

Issues and pull requests are welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

GNU General Public License v3 (GPL-3.0).
