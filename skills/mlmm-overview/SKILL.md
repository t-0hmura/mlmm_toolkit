---
name: mlmm-overview
description: Orientation for mlmm-toolkit — what it is, when to use it, and how it differs from generic ML/MM MD packages (three-layer ML/movable-MM/frozen ONIOM via PDB B-factor encoding, analytical hessian_ff full-system Hessian, microiteration with link-atom Jacobian coupling, AmberTools-driven MM parameterization). TRIGGER on first-touch / "what is mlmm-toolkit" / "should I use it" / "how does it compare to OpenMM / GROMACS / Sire" questions. SKIP when the user has already named a subcommand, an install issue, an output file, a structure format, or a cluster — sibling skills cover those.
---

# mlmm-toolkit Overview

## Purpose

`mlmm-toolkit` is a command-line toolkit for ML/MM ONIOM workflows on
solvated enzyme systems. It chains active-site definition, ML region
optimization with MM-environment relaxation, MEP search, TS
optimization, IRC validation, vibrational analysis, and an optional
DFT single-point — all driven by GPU-resident MLIP backends together
with an MM force-field layer.

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
4. **Bundled GPU pysisyphus fork.** Geometry / TS / IRC stay on the
   same device as the MLIP.
5. **AmberTools-driven MM parameterization.** `mlmm mm-parm` builds
   `parm7`/`rst7` from a PDB; `define-layer` assigns the ML / movable
   / frozen labels.

## When to use it

| Goal | Fit |
|---|---|
| Solvated enzyme reaction with explicit MM environment | Primary use case |
| Need link-atom + microiteration coupling | This toolkit |
| Need recursive multistep path search | `path-search` engine |
| Reuse a Gaussian g16 ONIOM input | `mlmm oniom-import` (then run downstream stages) |

## When *not* to use it

- Pure QM (DFT-only) cluster: a direct ORCA / Gaussian / Q-Chem
  workflow is leaner.
- Free-energy simulations (umbrella sampling, metadynamics): out of
  scope.

## Quick check

```bash
mlmm --version
mlmm --help              # lists the available subcommands
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
[path-opt]      single-pass MEP with ONIOM gradients (ML + MM coupling);
                recursive [path-search] with --refine-path
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

Supported MLIP backends:

| `-b` | Model | Notes |
|---|---|---|
| `uma` (default) | UMA-s-1.1 / UMA-m-1.1 (config strings: `uma-s-1p1` / `uma-m-1p1`) | Default for ML region |
| `mace` | MACE-OMOL-0 | Separate env (e3nn conflict) |
| `orb` | `orb_v3_conservative_omol` (Orb-v3-omol in papers) | Fast screening |
| `aimnet2` | AIMNet2 | Limited element coverage |

MM backend defaults to `hessian_ff` (CPU, analytical Hessian); the
finite-difference `openmm` backend is selectable via `--mm-backend openmm`.
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

See `mlmm-cli/SKILL.md` for per-subcommand specifics.

## Where the code lives

| File | What's there |
|---|---|
| `mlmm/cli/app.py` | Click entry point, subcommand registry |
| `mlmm/core/defaults.py` | All default kwarg dicts (MLMM_CALC_KW, MICROITER_KW, BFACTOR_*, IRC_KW, …) |
| `mlmm/backends/mlmm_calc.py` | The ONIOM ASE calculator (ML + MM gradient assembly, link-atom math) |
| `mlmm/workflows/extract.py` | Active-site extraction with layer assignment |
| `mlmm/workflows/define_layer.py` | B-factor → layer mapping helpers |
| `mlmm/workflows/mm_parm.py` | AmberTools tleap driver (parm7 / rst7) |
| `mlmm/workflows/oniom_export.py` / `oniom_import.py` | Gaussian g16 / ORCA ONIOM round-trip |
| `mlmm/workflows/all.py` | End-to-end pipeline |
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
