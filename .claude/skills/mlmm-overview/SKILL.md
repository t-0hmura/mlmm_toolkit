---
name: mlmm-overview
description: What mlmm-toolkit is, when to use it, and the design choices that distinguish it from cluster-only or generic QM/MM workflows (3-layer ML/movable-MM/frozen ONIOM, analytical hessian_ff, microiteration, AmberTools-driven MM parameterization).
---

# mlmm-toolkit Overview

## What it is

`mlmm-toolkit` is a command-line toolkit for ML/MM ONIOM workflows on
solvated enzyme systems. It chains active-site definition, ML region
optimization with MM-environment relaxation, MEP search, TS
optimization, IRC validation, vibrational analysis, and an optional
DFT single-point — all driven by the same MLIP backends as
`pdb2reaction` plus an MM force-field layer.

The design choices that make it distinct:

1. **3-layer ONIOM via PDB B-factor encoding.** The B-factor field
   classifies every atom into ML (0.0), movable-MM (10.0), or frozen
   (20.0). One PDB → one defined system, no separate topology files
   for the partitioning.
2. **`hessian_ff` analytical Hessian for MM.** Custom analytical
   Hessian for the MM force field (instead of finite-difference),
   used both in microiteration and in the ML/MM-coupled freq step.
3. **Microiteration outer/inner loop.** ML region geometry update
   alternates with MM relaxation; outer ML steps see a relaxed MM
   environment.
4. **Bundled GPU pysisyphus fork (same as `pdb2reaction`).**
   Geometry / TS / IRC stay on the same device as the MLIP.
5. **AmberTools-driven MM parameterization.** `mlmm mm-parm` builds
   `parm7`/`rst7` from a PDB; `define-layer` assigns the ML / movable
   / frozen labels.

## When to use it

| Goal | Fit |
|---|---|
| Solvated enzyme reaction with explicit MM environment | Primary use case |
| Cluster model only (no MM region) | Use `pdb2reaction` instead |
| Need link-atom + microiteration coupling | This toolkit |
| Need recursive multistep path search | Both share the same `path-search` engine |
| Reuse a Gaussian g16 ONIOM input | `mlmm oniom-import` (then run downstream stages) |

## When *not* to use it

- Pure QM (DFT-only) cluster: `pdb2reaction dft` or a direct ORCA /
  Gaussian / Q-Chem workflow.
- Free-energy simulations (umbrella sampling, metadynamics): out of
  scope.
- ML-only enzyme reaction, no MM: `pdb2reaction` is leaner.

## Quick check

```bash
mlmm --version
mlmm --help              # 22 subcommands listed
mlmm all --help          # end-to-end pipeline
```

If `mlmm` is not on PATH or imports fail, see
`mlmm-install-backends/SKILL.md`.

## Pipeline at a glance

```
PDB(s)          (B-factor: 0.0=ML, 10.0=movable-MM, 20.0=frozen)
  │
  ▼
[mm-parm]       AmberTools tleap → parm7 / rst7
  │
  ▼
[define-layer]  expand / refine / verify the ML/MM/Frozen labels
  │
  ▼
[path-search]   recursive MEP with ONIOM gradients (ML + MM coupling)
  │
  ▼
[tsopt]         TS refinement per segment
  │
  ▼
[irc]           forward/backward IRC + endpoint LBFGS
  │
  ▼
[freq]          analytical-Hessian ONIOM frequencies + QRRHO thermo
  │
  ▼
[dft]           (optional) single-point DFT on ML region only
```

Each step is also available as its own subcommand. `mlmm all` chains
the whole pipeline.

## Backend choices

Same MLIP family as `pdb2reaction`:

| `-b` | Model | Notes |
|---|---|---|
| `uma` (default) | UMA-s-1.1 / UMA-m-1.1 (config strings: `uma-s-1p1` / `uma-m-1p1`) | Default for ML region |
| `mace` | MACE-OMOL-0 | Separate env (e3nn conflict) |
| `orb` | `orb_v3_conservative_omol` (Orb-v3-omol in papers) | Fast screening |
| `aimnet2` | AIMNet2 | Limited element coverage |

MM backend is fixed: `hessian_ff` (CPU only, with analytical Hessian).
DFT (optional) uses PySCF / GPU4PySCF.

## ML/MM-aware CLI conventions

Every ML/MM-evaluating subcommand (`opt`, `tsopt`, `path-search`,
`scan`, `freq`, `irc`, `dft`, `all`, …) takes:

| flag | purpose |
|---|---|
| `-i, --input` | Full-enzyme PDB (or XYZ + `--ref-pdb`) |
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (default on) |
| `--model-indices` | Alternative to `--model-pdb`: comma-separated atom indices (e.g. `1-50,75,100-110`) |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` / `-l, --ligand-charge` / `-m, --multiplicity` | ML region charge / spin |
| `-b, --backend` | ML backend (uma / orb / mace / aimnet2) |

These supersede the cluster-only conventions used by `pdb2reaction`.
See `mlmm-cli/SKILL.md` for per-subcommand specifics.

## Sibling project: pdb2reaction

`pdb2reaction` shares the same MLIP backends, `pysisyphus` fork, and
`thermoanalysis` core but targets **cluster models** (no MM region).

| Use case | Toolkit |
|---|---|
| Cluster, no MM | `pdb2reaction` |
| Solvated enzyme with MM around ML | `mlmm-toolkit` |
| Need automatic Amber parameterization | `mlmm mm-parm` |
| Recursive multistep path search | Both |

The two toolkits **cannot share a Python environment** (incompatible
`pysisyphus` forks and `e3nn` versions). Keep separate `conda env`s.

## Where the code lives

| File | What's there |
|---|---|
| `mlmm/cli.py` | Click entry point, subcommand registry |
| `mlmm/defaults.py` | All default kwarg dicts (MLMM_CALC_KW, MICROITER_KW, BFACTOR_*, IRC_KW, …) |
| `mlmm/mlmm_calc.py` | The ONIOM ASE calculator (ML + MM gradient assembly, link-atom math) |
| `mlmm/extract.py` | Active-site extraction with layer assignment |
| `mlmm/define_layer.py` | B-factor → layer mapping helpers |
| `mlmm/mm_parm.py` | AmberTools tleap driver (parm7 / rst7) |
| `mlmm/oniom_export.py` / `oniom_import.py` | Gaussian g16 / ORCA ONIOM round-trip |
| `mlmm/all.py` | End-to-end pipeline |
| bundled `hessian_ff/` | Analytical-Hessian MM force field |
| bundled `pysisyphus/` | GPU-tensor pysisyphus fork |
| bundled `thermoanalysis/` | QRRHO thermochemistry |

## Navigation map of the skill set

| You want to … | Read |
|---|---|
| Pick a subcommand and run it | `mlmm-cli/SKILL.md` then the per-subcommand md |
| Read or edit a `.pdb` / `.xyz` / `.gjf` / `.parm7` | `mlmm-structure-io/{SKILL,pdb,xyz,gjf,parm7}.md` |
| Decide charge / multiplicity for a substrate | `mlmm-structure-io/charge-multiplicity.md` |
| Install the toolkit, an MLIP backend, AmberTools, or DFT | `mlmm-install-backends/` |
| Build an analytical recipe (full ONIOM / scan-list / ts-only) | `mlmm-workflows-output/SKILL.md` |
| Submit on PBS / SLURM | `mlmm-hpc/SKILL.md` |
| Detect the cluster / GPU / scheduler you're on | `mlmm-env-detect/SKILL.md` |