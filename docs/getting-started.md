# Getting Started

## Overview

`mlmm-toolkit` is a Python CLI for **ML/MM ONIOM** analyses of enzymatic reactions. It replaces the QM region of conventional QM/MM with a machine-learning interatomic potential (MLIP, default UMA; `orb` / `mace` / `aimnet2` via `-b`), with the surrounding protein under the bundled Amber force field `hessian_ff`. The ONIOM decomposition is:

```
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

A single command generates a useful initial reaction path:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'                  # MEP only
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo --dft   # full
```

Given (i) ≥ 2 PDBs (R → ... → P), (ii) one PDB with `--scan-lists`, or (iii) one TS candidate with `--tsopt`, `mlmm all` defines the ML region, runs `mm-parm` + `define-layer`, performs MEP search (GSM), and optionally chains TS optimisation, IRC, frequencies, and single-point DFT.

```{important}
- Input PDBs must already contain **hydrogen atoms**. The "Input prep checklist" below covers the common pitfalls.
- Multiple PDBs must share the same atoms in the same order (only coordinates differ).
- Most per-stage subcommands require `--parm` (from `mm-parm`) and `--model-pdb` (from `extract` / `define-layer`); `mlmm all` generates both automatically.
```

For background concepts (3-layer system, link atoms, microiteration, units), read [Concepts & Workflow](concepts.md). For symptom-first diagnosis, jump to [Troubleshooting](troubleshooting.md) or [Common Error Recipes](recipes-common-errors.md).

### CLI conventions

| Convention | Example |
|---|---|
| Residue selector | `'SAM,GPP'` or `'A:123,B:456'` |
| Charge mapping | `-l 'SAM:1,GPP:-3'` |
| Atom selector | `'TYR,285,CA'` or `'TYR 285 CA'` |

Full table: [CLI Conventions](cli-conventions.md).

### Input prep checklist

- **Hydrogens present.** `mlmm` does not auto-protonate. Add with AmberTools `reduce`, OpenMM `Modeller.addHydrogens`, `pdb2pqr --ff=AMBER`, Open Babel `obabel -h`, or `mlmm mm-parm --add-h` (PDBFixer wrapper). Apply the same tool to every input to keep atom order consistent.
- **Match `-l RES:CHARGE` to the H count actually in the file** (e.g. SAM with 23 H = `SAM:1` cation, 22 H = `SAM:0` neutral). Mismatch breaks `antechamber` with an odd-electron sqm failure — do not re-protonate "to look canonical".
- **R/P atom order must match.** In PyMOL, tick *Original atom order* on export.
- **Multi-chain enzymes need `TER` records** between chains so `tleap` segments them correctly.
- **Per-stage subcmds**: `-q/--charge` is the **ML-region (ONIOM model-system) net charge**, not the full-system charge. Passing the whole-enzyme charge silently builds a wrong ML region.

---

## Installation

```bash
# 0. Clone the repo (skip if you only want `pip install mlmm-toolkit` once published)
git clone https://github.com/t-0hmura/mlmm_toolkit.git && cd mlmm_toolkit

# 1. New env + AmberTools + CUDA-enabled PyTorch (match your CUDA runtime)
conda create -n mlmm-toolkit python=3.11 -y && conda activate mlmm-toolkit
conda install -c conda-forge ambertools pdbfixer -y
pip install torch --index-url https://download.pytorch.org/whl/cu129

# 2. mlmm-toolkit (editable from a local clone, or `pip install mlmm-toolkit` once published)
pip install -e .
# Optional MLIP extras: pip install -e ".[orb]"  /  ".[aimnet]"  /  ".[dft]"  /  ".[mcp]"
# MACE: install in a dedicated env (incompatible with UMA via e3nn==0.4.4 vs >=0.5)

# 3. (UMA backend only) Authenticate Hugging Face once
#    Accept the FAIR Chemistry License v1 at https://huggingface.co/facebook/UMA, then:
hf auth login                                                # interactive
# OR: export HF_TOKEN=hf_xxx && hf auth login --token "$HF_TOKEN" --add-to-git-credential   # CI / HPC

# 4. Verify
mlmm --version
```

### Optional components

| Component | When to add | Install |
|---|---|---|
| `hessian_ff` native build | If you see a "native extension not available" warning (JIT compilation usually handles it) | `cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make` (install `ninja` first on most clusters: `conda install -c conda-forge ninja -y`) |
| `cyipopt` | DMF MEP backend for the standalone `path-search` / `path-opt` subcommands (`--mep-mode dmf`); `mlmm all` is GSM-only | `conda install -c conda-forge cyipopt -y` |
| xTB | `--embedcharge` (xTB point-charge embedding) | `conda install -c conda-forge xtb -y` (custom binary: set `xtb_cmd` in YAML) |
| Plotly Chrome | Static PNG export beyond default `kaleido` | `plotly_get_chrome -y` (~150 MB) |
| HPC `cuda/<X.Y>` module | HPC clusters using environment modules | Load **before** `pip install torch`; match X.Y to your wheel (`cu126` ↔ 12.6, `cu129` ↔ 12.9) |

If you switch runtime environments (node / container / Python / PyTorch), rebuild `hessian_ff` in the new env. Detailed HPC job-script templates: [docs/device-hpc.md](device-hpc.md).

If no usable conda `xtb` package is available, build xTB from source (requires GCC >= 10):

```bash
git clone --depth 1 https://github.com/grimme-lab/xtb.git
cd xtb
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
make -C build -j8
```

---

## Quickstart routes

- [Quickstart: `mlmm all`](quickstart-all.md) — multi-structure MEP
- [Quickstart: `mlmm` scan-spec route](quickstart-scan-spec.md) — single structure with staged bond scans
- [Quickstart: validate TS with `mlmm tsopt`](quickstart-tsopt-freq.md) — TS-only mode

## Typical workflow

```text
1. extract       — Define ML region from full protein-ligand PDB
2. mm-parm       — Generate Amber parm7/rst7 topology (requires AmberTools)
3. define-layer  — Assign 3-layer ML/MM partitioning (B-factor encoding)
4. path-search   — MEP search (recursive; `--no-refine-path` for single-pass `path-opt`)
5. tsopt         — Transition state optimisation
6. freq          — Vibrational analysis + thermochemistry
7. dft           — Single-point DFT energy refinement
```

`mlmm all` orchestrates all seven; each step is also a standalone subcommand for debugging or custom flows.

## Main workflow modes

| Mode | Trigger | Use when |
|---|---|---|
| Multi-structure MEP | `-i R.pdb P.pdb [I1.pdb ...]` | You have ≥ 2 endpoints / intermediates (docking, MD, manual modelling). |
| Staged scan | `-i ONE.pdb --scan-lists '[...]' [ '[...]' ...]` | You'd rather define reaction coordinates than provide multiple endpoints. |
| TS-only | `-i TS_CANDIDATE.pdb --tsopt` | You already have a TS guess and want `tsopt → IRC → freq`. |

`mlmm [OPTIONS]` is equivalent to `mlmm all [OPTIONS]` — `all` is the default subcommand, so the bare `mlmm -i ...` examples below run the full `all` workflow.

```bash
# Multi-structure MEP (richer)
mlmm -i R.pdb I1.pdb I2.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
     --out-dir ./result_all --tsopt --thermo --dft

# Staged scan
mlmm -i R.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
     --scan-lists '[("TYR 285 CA","MMT 309 C10",2.20),("TYR 285 CB","MMT 309 C11",1.80)]' \
                  '[("TYR 285 CB","MMT 309 C11",1.20)]'

# TS-only
mlmm -i TS_CANDIDATE.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

Each tuple `(i, j, target_Å)` accepts a PDB atom selector (`'TYR,285,CA'` — space/comma/slash/backtick/backslash) or a 1-based atom index. Pass multiple stages as multiple literals after a single `--scan-lists` flag.

```{important}
Single-input runs require **either** `--scan-lists` (staged scan → GSM) **or** `--tsopt` (TS-only). A bare `-i ONE.pdb` will not trigger a full workflow.
```

## Multi-backend examples

```bash
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b orb         # ORB
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 -b mace        # MACE
mlmm opt -i ml_region.pdb --parm real.parm7 --model-pdb ml.pdb -q 0 --embedcharge  # xTB embedding
```

## DFT refinement via Gaussian / ORCA (hand-off)

`mlmm-toolkit` provides a round-trip hand-off — Gaussian or ORCA must be licensed and on `PATH` separately:

```bash
# 1. ML/MM TS refinement
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1

# 2. Export to Gaussian ONIOM (.com)
mlmm oniom-export --mode g16 --parm real.parm7 -i result_tsopt/final_geometry.pdb \
     --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1 --method "wB97XD/def2-TZVPD"

# 3. Run externally (ORCA via --mode orca also supported)
g16 < ts_refine.com > ts_refine.log

# 4. Pull DFT-refined geometry back in
mlmm oniom-import -i ts_refine.com --ref-pdb ml_region.pdb -o ts_dft

# 5. Continue freq / IRC on the refined geometry
mlmm freq -i ts_dft_layered.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1
mlmm irc  -i ts_dft_layered.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1
```

Full flag references: [oniom-export](oniom-export.md), [oniom-import](oniom-import.md), [oniom-gaussian](oniom-gaussian.md), [oniom-orca](oniom-orca.md).

## Common options

| Option | Description |
|---|---|
| `-i, --input PATH...` | Input structures. Dispatch trichotomy is the same as the "Main workflow modes" table above. |
| `-c, --center TEXT` | Substrate / extraction center (residue names `'SAM,GPP'`, residue IDs `A:123,B:456`, or PDB paths). |
| `-l, --ligand-charge TEXT` | Charge mapping (`'SAM:1,GPP:-3'`) or single integer. |
| `-q, --charge INT` / `-m, --multiplicity INT` | ML-region net charge / spin multiplicity. |
| `-s, --scan-lists TEXT...` | Staged distance scans for single-input runs (literals or YAML/JSON file). |
| `-o, --out-dir PATH` | Top-level output directory. |
| `--tsopt` / `--thermo` / `--dft` | TS optimisation + IRC / vibrational analysis / single-point DFT. |
| `--refine-path` / `--no-refine-path` | Recursive `path-search` (default) vs single-pass `path-opt`. |
| `-b, --backend uma\|orb\|mace\|aimnet2` | MLIP backend (default `uma`). |
| `--embedcharge` | xTB point-charge embedding correction (default off). |
| `--hessian-calc-mode Analytical\|FiniteDifference` | ML Hessian mode. `Analytical` is UMA-only; recommended when VRAM allows. |

DMF MEP is selectable only via the standalone `path-search` / `path-opt` subcommands (`--mep-mode dmf`); `mlmm all` always uses GSM (passing `--mep-mode` to `mlmm all` is silently ignored).

Full option matrix and YAML schema: [YAML Reference](yaml-reference.md). Subcommand-by-subcommand table: [README "CLI Subcommands"](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md#cli-subcommands).

## Run summaries

Every `mlmm all` run writes `summary.log` (human) + `summary.json` (machine) with the CLI command, global MEP statistics, per-segment barriers / bond changes, and MLIP / thermo / DFT energies (when enabled). Per-segment `segments/seg_NN/` subdirectories carry their own summaries.

## Getting help

```bash
mlmm --help                            # top-level
mlmm <subcommand> --help               # core options
mlmm <subcommand> --help-advanced      # full option set
```

## Driving from an AI coding agent

`mlmm-toolkit` ships `skills/` with agent-readable instructions. Copy `skills/` into your project as `.claude/skills/` (or merge into `~/.claude/skills/`) for Claude Code / Cursor / OpenCode pickup.

```{warning}
This software is still under development. Please use it at your own risk.
```
