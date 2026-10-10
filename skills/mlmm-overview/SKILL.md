---
name: mlmm-overview
description: "Orientation for mlmm-toolkit: which `all` mode fits the available structures (Endpoint mode for two or more full-system structures, Scan-list mode for one structure with -s, TS-only mode for one TS candidate), when to run stage by stage instead, and how to judge each stage. Covers the three-layer ML/MM model (ML, Movable-MM, and Frozen-MM, set by PDB B-factors 0/10/20), TS-candidate strategy and retries (`ts-strategy.md`), reading `summary.json`, `result.json`, and the output tree (`outputs.md`), and where the source code lives. TRIGGER on first-touch questions, choosing an all mode or workflow, barrier or imaginary-frequency questions, mutant comparisons, extracting numbers for a paper, or locating code. SKIP when the user already named a subcommand, an install issue, a structure format, ML-region or layer design (mlmm-model-setup), or a cluster job; sibling skills cover those."
---

# mlmm-toolkit overview

`mlmm all` picks its mode from the inputs: two or more full-system structures in reaction order → Endpoint mode (`all-endpoint-mep.md`); one structure with `-s` → Scan-list mode (`all-scan-list.md`); one TS candidate with `--tsopt` → TS-only mode (`all-ts-only.md`). `-c` sets the ML region, and `all` builds the parm7 and layers unless you pass them; run the stages one by one when you want to check each result first.

## Pick an all mode

| Structures you have | Mode | Read |
|---|---|---|
| R and P, with any intermediates between them, in reaction order | Endpoint mode | [all-endpoint-mep.md](../mlmm-cli/all-endpoint-mep.md) |
| R only, plus the bonds to drive, given with `-s` | Scan-list mode | [all-scan-list.md](../mlmm-cli/all-scan-list.md) |
| One TS candidate (`-i`), with `--tsopt` and no `-s` | TS-only mode | [all-ts-only.md](../mlmm-cli/all-ts-only.md) |

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -o result_mep
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.50),("CS1 SAM 320","SD SAM 320",3.30)]' \
       '[("C7 GPP 321","H11 GPP 321",2.90),("OE2 GLU 186","H11 GPP 321",1.00)]' \
    --tsopt --thermo -o result_scan
mlmm all -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo -o result_ts
```

A run finished every requested stage when each TS prints `[Imaginary modes] n=1 (...)` and the console ends with `Scientific status: success` under `====== Pipeline summary ======`. Then check that `segments/seg_NN/reactant.*` and `product.*` are the R and P you intended ([outputs.md](outputs.md)).

- **Endpoint mode**: `path-opt` finds one MEP per neighbouring pair of inputs. `--refine-path` runs the recursive `path-search` instead, which splits a multistep reaction where bonds change, so `n_segments` can exceed the number of pairs. Either way, a segment is a candidate step until its TS and IRC are checked.
- **Scan-list mode**: each literal after `-s` is one stage; the stage ends become the inputs of the MEP search.
- **TS-only mode**: for a candidate from another code or an earlier run, `all` runs `tsopt`, the IRC, and both endpoint optimizations, without an MEP search. The same chain can be run by hand ([Run stage by stage](#run-stage-by-stage)).
- **DFT//MLIP/MM**: `--dft` (with `--tsopt`) adds DFT single points on the ML region of R, TS, and P; for standalone runs see [DFT//MLIP/MM on the TS candidate](../mlmm-cli/dft.md#dftmlipmm-on-the-ts-candidate).

Inputs are full-system `.pdb`, `.cif`, `.mmcif`, or `.xyz` files; an XYZ needs `--ref-pdb`. With `-c`, `all` cuts the ML region around the given residues; without `-c`, it uses the layers in the input B-factors. `--parm7` skips `mm-parm`, and `--model-pdb` takes priority over `-c` and the B-factors. The bundled `1.R.pdb` has an empty chain column, so atoms are written as residue name, number, and atom name; with chains, write `A:SAM:320:CS1`. What to put in the ML region and the layers is in [mlmm-model-setup](../mlmm-model-setup/SKILL.md).

Pitfalls: two or more structures with `-s` is an error, and so is one structure with neither `-s` nor `--tsopt`. One structure with both `-s` and `--tsopt` runs Scan-list mode. `--thermo` and `--dft` need `--tsopt`. Without `-c`, `--no-detect-layer` with no `--model-pdb` stops with an error.

## Pipeline at a glance

```text
full-system structure(s)   B-factor layers optional: 0 = ML, 10 = Movable-MM, 20 = Frozen-MM
  │
[extract]        ML region around -c (skipped without -c: B-factor layers or --model-pdb)
[mm-parm]        AmberTools tleap → parm7 / rst7 (skipped with --parm7)
[define-layer]   ML / Movable-MM / Frozen-MM written into the B-factors
[scan]           one structure with -s: staged restrained scan
[path-opt]       MEP with ONIOM gradients; recursive [path-search] with --refine-path
[tsopt]          TS optimization of each segment's HEI (--tsopt)
[irc]            IRC in both directions, then optimization of both ends (--tsopt)
[freq]           partial-Hessian (PHVA) frequencies and QRRHO thermochemistry (--thermo)
[dft]            DFT single points on the ML region only (--dft)
```

TS-only mode skips the scan and the MEP. Each step is also its own subcommand. Run alone, the setup order is `mm-parm → extract → define-layer`, so the ML region is cut from the PDB that matches the parm7 ([cli/extract.md](../mlmm-cli/extract.md)).

## Run stage by stage

Run the subcommands one by one instead of a single `mlmm all` when you want to judge each stage before spending GPU time on the next: confirm the MEP found the right bond changes before optimizing a TS, and confirm the TS before thermochemistry or DFT. Every ML/MM stage needs the same `--parm7`, the same ML region (`--model-pdb`, or the B-factor layers), and the same `-l` / `-q` / `-m`; pass them on every command. After each stage, read `execution_status` and `scientific_status` in its `result.json` or `summary.json`.

**Stage 0, setup** (only when starting from a raw full-system PDB; most campaigns start from prepared, layered R and P PDBs and a parm7):

```bash
mlmm mm-parm -i input.pdb -l 'SAM:1,GPP:-3' --out-prefix real
mlmm extract -i real.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o ml_region.pdb
mlmm define-layer -i real.pdb --model-pdb ml_region.pdb -o R_layered.pdb
```

GATE: `[mm-parm] Wrote:` lists `real.pdb` and `real.parm7`; `real.pdb` has filled element columns and the same atoms in the same order as `real.parm7`; the layered PDB carries 0/10/20 on the intended atoms.

**Pre-optimization.** `all`, `path-opt`, and `path-search` optimize each endpoint without restraints before the MEP (`--preopt`, on by default; in Scan-list mode `--scan-preopt` follows it), and this can already move a proton, break a weak bond, or complete part of the reaction. Before reading the MEP, run `bond-summary` between each input and its pre-optimized structure (`_work/scan/preopt/result.*` in Scan-list mode; otherwise `final_geometry.*` in the `initNN_*_opt/` directories that `path-search` writes, under `_work/path_search/` in `all --refine-path`, or `preopt/endNN/final_geometry.xyz` of each `path-opt` run, under `_work/path_opt/seg_NN_mep/` in `all`) and compare the moving H atoms; judge later bond changes against this optimized structure, not the raw MD frame, whose short contacts can count as bonds. A non-reacting bond that breaks points to the ML region or the backend: enlarge the ML region or change the backend.

If the chemistry you meant to start from changed, hold the bonds to keep with `opt --distance-restraint '[(i, j)]'` (a pair without a target keeps its current distance); start a scan from that result with `all --no-scan-preopt`, and for an MEP endpoint run an unrestrained `opt` from it and check the bonds again before passing it to `all`. If you continue from the changed R instead, treat it as another chemical state: check the ML/MM boundary, the layers, and protonation, and do not rank its barriers with candidates that kept the original state.

**Stage 1, MEP**:

```bash
mlmm path-search -i R_layered.pdb P_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -l 'SAM:1,GPP:-3' -o ps
```

GATE: `ps/summary.json` has `"scientific_status": "success"`, and every required `stage_outcomes[]` entry is usable. Read `n_segments` and each segment's `bond_changes`: the intended bonds must form and break on the right atoms. Fix the chemistry or the inputs before any TS work if the segmentation or the bond changes are wrong.

**Stage 2, TS and IRC** for each reactive segment, starting from `ps/hei_seg_NN.xyz`:

```bash
mlmm tsopt -i ps/hei_seg_NN.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -l 'SAM:1,GPP:-3' --out-json -o seg_NN/tsopt
mlmm irc -i seg_NN/tsopt/final_geometry.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -l 'SAM:1,GPP:-3' --out-json -o seg_NN/irc
```

GATE for the TS: in `seg_NN/tsopt/result.json`, `optimization_status` is `converged`, `hessian_status` is `completed`, `saddle_validation` is `first_order` (`n_imaginary_modes` 1), and the imaginary mode in `vib/imag_*_trj.xyz` moves the reacting atoms. Standalone `tsopt` sets `scientific_status` from convergence alone, so read n_imag yourself. A run that stops at max cycles without converging computes no Hessian and reports no n_imag; a run stopped on an energy plateau (`--stop-plateau`, `stalled`) computes the Hessian and reports n_imag. If n_imag is not 1, see [ts-strategy.md](ts-strategy.md).

GATE for the IRC: the console line `Transition vector is mode 0 with wavenumber … cm⁻¹.` shows a negative wavenumber. The IRC has no success verdict of its own: read the frame counts and the stop reason of each branch in its `result.json`, then optimize both ends:

```bash
mlmm opt -i seg_NN/irc/forward_first.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -l 'SAM:1,GPP:-3' --out-json -o seg_NN/end_forward
mlmm opt -i seg_NN/irc/backward_last.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -l 'SAM:1,GPP:-3' --out-json -o seg_NN/end_backward
```

Both optimizations must converge. Even if the IRC did not converge, the result is usable when the optimized ends are the intended R and P; forward and backward do not tell which is which, so compare the bonds and the moving H atoms. Use each `final_geometry.xyz` downstream; standalone commands do not write the `segments/seg_NN/reactant.*` and `product.*` that `all` writes. A TS with n_imag = 1 does not by itself establish the elementary step.

**Stage 3, thermochemistry** (optional, as `all --thermo`): run `freq` on R, TS, and P for the Gibbs profile.

```bash
mlmm freq -i seg_NN/tsopt/final_geometry.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -l 'SAM:1,GPP:-3' --out-json -o seg_NN/freq_TS
```

**Stage 4, DFT//MLIP/MM** (optional, as `all --dft`): repeat for both optimized ends with the same settings.

```bash
mlmm dft -i seg_NN/tsopt/final_geometry.xyz --ref-pdb R_layered.pdb --parm7 real.parm7 \
    --model-pdb ml_region.pdb -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' --out-json -o seg_NN/dft_TS
```

GATE: each DFT `result.json` has `"converged": true`; an SCF that does not converge exits with code 1.

**Stage 5, energy diagram**:

```bash
mlmm energy-diagram -i 0.0 -i 21.5 -i -0.7 --label-x R --label-x TS --label-x P -o diagram.png
```

Pitfalls and recovery:

- After a walltime stop, rerun the same commands; for `all`, repeat the original command with `--resume-segment N` ([all.md](../mlmm-cli/all.md)).
- On any status other than `success`, read `summary.log`, then the `result.json` of the failed stage. Inside an `all` run these exist under `ts/`, `irc/`, and `endpoint_opt/`, not under `freq/` or `dft/`.
- If the Bofill Hessian update of an IRC runs out of GPU memory, rerun with `PYSIS_BOFILL_CPU_OFFLOAD=1`. It does the update on the CPU at the cost of host memory and two full-matrix transfers, and does not help a frequency Hessian that runs out of memory.
- `all` reuses the final tsopt Hessian in IRC and in the TS `freq` through an in-process cache, so separate `tsopt` → `irc` → `freq` commands build the same dense Hessian up to three times. Add `--dump-hess ts_hess.npy` to `tsopt` (not with `--skip-final-freq`) and pass `--read-hess ts_hess.npy` to `irc` (which needs `irc.hessian_init: calc`, the default) and to `freq`, both run on the tsopt `final_geometry.*`.

## What it does

`mlmm-toolkit` runs ML/MM ONIOM reaction studies on solvated enzymes: ML-region selection, MEP search, TS optimization, IRC, vibrational analysis, and optional DFT single points, with an MLIP for the ML region and an Amber force field for the rest.

1. **Three layers in the PDB B-factors**: every atom is ML (0), Movable-MM (10), or Frozen-MM (20), so one PDB and one parm7 define the system. The energy is `E_MM(real) + E_ML(model) − E_MM(model)`, and link hydrogens cap each parm7 bond that crosses the ML/MM boundary.
2. **Microiteration**: an optimizer step on the ML atoms (with the MM atoms bonded to them) alternates with an L-BFGS relaxation of the other movable MM atoms under the force field alone. It is on by default in `opt --opt-mode hess` and in the Hessian TS optimizers of `tsopt`; `--no-microiter` turns it off.
3. **MM Hessian**: the MM engine is `hessian_ff` on the CPU by default (`--mm-backend openmm` for OpenMM). The MM Hessian uses finite differences by default; YAML `calc.mm_fd: false` switches to the analytical `hessian_ff` Hessian.
4. **Bundled pysisyphus**: a GPU-capable copy of pysisyphus runs the geometry optimizations, TS searches, and IRC integrations.
5. **AmberTools topology**: `mm-parm` builds the parm7 and rst7 from a PDB with tleap, and `define-layer` writes the three layers.

## When to use it

- A reaction in a solvated enzyme with an explicit MM environment: this is the main use.
- A study that needs link atoms and microiteration between the ML and MM regions.
- A multistep reaction whose steps must be found by a recursive path search (`--refine-path`).
- An existing Gaussian or ORCA ONIOM input to continue from: `mlmm oniom-import`, then the later stages ([cli/oniom.md](../mlmm-cli/oniom.md)).

## When not to use it

- A pure QM cluster model with DFT only: an ORCA or Gaussian workflow on its own is leaner.
- Free-energy sampling (umbrella sampling, metadynamics): out of scope. The Gibbs barrier from `--thermo` is the static TS of one ML/MM structure with QRRHO corrections, not a potential of mean force: keep its method label when you quote it, and for a condensed-phase free-energy barrier, pass the validated R, TS, and P to a free-energy sampling tool.

## Quick check

```bash
mlmm --version
mlmm --help              # lists the subcommands
mlmm all --help          # the end-to-end pipeline
```

If `mlmm` is not on PATH or an import fails, see [Verify the install](../mlmm-install/SKILL.md#verify-the-install).

## ML/MM layers in one paragraph

Every ML/MM command takes the full system with `-i`, the topology with `--parm7`, and the ML region from `--model-pdb` or the B-factors; the shared options are in [Shared ML/MM conventions](../mlmm-cli/SKILL.md#shared-mlmm-conventions), the B-factor encoding in [ML region and layers](../mlmm-model-setup/SKILL.md#ml-region-and-layers), and how to choose, trim, and enlarge the region in [mlmm-model-setup](../mlmm-model-setup/SKILL.md).

## Backends

`-b` selects the ML backend (`uma` by default, `orb`, `mace`, `aimnet2`, or `dft`); the table and the install steps are in [Choose a backend](../mlmm-install/SKILL.md#choose-a-backend), and DFT single points in [cli/dft.md](../mlmm-cli/dft.md).

## Where the code lives

The package body `mlmm/` has one directory per layer; `pysisyphus/`, `thermoanalysis/`, and `hessian_ff/` install as separate top-level packages next to it.

| Concern | Open |
|---|---|
| Subcommand list and entry point | `mlmm/cli/app.py` |
| Shared option decorators | `mlmm/cli/common_options.py`, or the subcommand file itself |
| Default of a flag | its Click definition, and shared values in `mlmm/core/defaults.py` |
| Body of a subcommand | `mlmm/workflows/<subcommand>.py` for the calculation commands (`all.py`, `extract.py`, `mm_parm.py`, `define_layer.py`, `oniom_export.py`, `oniom_import.py`, …); the utilities are in `mlmm/io/` and `mlmm/domain/` (`_LAZY_SUBCOMMANDS` in `mlmm/cli/app.py` maps each command to its module) |
| ONIOM calculator, link atoms, MLIP backends | `mlmm/backends/mlmm_calc.py` |
| Bond changes and other chemistry helpers | `mlmm/domain/` |
| `summary.json`, energy diagrams, trajectories | `mlmm/io/` |
| Analytical MM Hessian | `hessian_ff/` |
| Optimizer, TS, and IRC internals | `pysisyphus/` |
| QRRHO thermochemistry | `thermoanalysis/` |
| MCP server | `mlmm/mcp/` ([mlmm-mcp](../mlmm-mcp/SKILL.md)) |
| Chemistry rules | search for `# CHEMISTRY-RULE:` markers |

The layer map, the import rules, and the invariants to keep are in [`docs/architecture.md`](../../docs/architecture.md); contributor recipes in [`CONTRIBUTING.md`](../../CONTRIBUTING.md). The import graph is checked by [`check_import_graph.py`](../../.github/scripts/check_import_graph.py), and the chemistry markers by [`check_engineering_markers.py`](../../.github/scripts/check_engineering_markers.py).

## Where to go next

- [ts-strategy.md](ts-strategy.md): studying a mechanism (hypothesis, TS precision, routes to a candidate, splitting the reaction, wrong n_imag, a TS that does not come out, multistep paths, comparisons, barriers).
- [outputs.md](outputs.md): `summary.json`, `result.json`, and the output tree.
- [mlmm-cli](../mlmm-cli/SKILL.md): running and judging each subcommand.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md): formats, residue and atom selection, layer encoding, charge and multiplicity, and building, trimming, and enlarging the ML region and the layers.
- [mlmm-install](../mlmm-install/SKILL.md): the core, backends, AmberTools, CUDA, and checking an unknown environment.
- [mlmm-hpc](../mlmm-hpc/SKILL.md): job scripts.
- [mlmm-mcp](../mlmm-mcp/SKILL.md): the MCP tools.
- [colab-local-gpu-runtime](../colab-local-gpu-runtime/SKILL.md): a Colab local runtime.
