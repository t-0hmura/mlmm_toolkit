# `all`

`mlmm all` runs the end-to-end ML/MM enzymatic-reaction workflow from full-system structures. Internally it coordinates active-site extraction, MM topology preparation, ML/MM layer assignment, an optional scan, MEP search (single-pass `path-opt` by default; recursive `path-search` with `--refine-path`), and optional post-processing (TS optimization, EulerPC IRC, thermochemistry, single-point DFT, and DFT//MLIP/MM diagrams). The default MLIP backend for the ML region is UMA; choose an alternative with `-b/--backend`.

That sequence describes `all`'s internally managed stages. For reusable files,
request a distinct `mm-parm` output prefix, then run `extract` and
`define-layer` on that exported PDB. It has the same atom identity and order as
the generated `parm7`; `mm-parm` fills its missing element columns.

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

`all` runs in one of three modes, chosen by what you pass:

- **Multi-structure MEP** — give at least two full structures in reaction order to drive a GSM (default) or DMF MEP search across the supplied structures.
- **Single-structure scan-defined workflow** — give one full structure plus `--scan-lists`. One literal defines one stage; several tuples within it are advanced concertedly. The relaxed stage endpoints become the MEP input series.
- **TSOPT-only** — give a single full structure and set `--tsopt` (no `--scan-lists`) to run TS optimization directly, with no MEP search.

Inputs may also be `.cif` / `.mmcif`; computation uses a temporary internal
PDB and public CIF companions restore the original identifiers.

```{important}
`--tsopt` produces **TS candidates** and reports numerical optimization and
terminal saddle order separately. `all` proceeds to IRC only when optimization
converged, terminal PHVA completed, and a negative reaction direction is
available. A converged higher-order stationary point may continue through
warning-labelled **diagnostic** IRC, but it is not a certified first-order TS.
Actual optimizer non-convergence, zero imaginary modes, failed/unavailable PHVA,
or no valid negative root stops after preserving TS artifacts and before IRC.
`--skip-final-freq` also stops before IRC because the reaction direction cannot
be validated. Always inspect the modes and endpoint connectivity.
```

## Examples

Command form:

```bash
mlmm all -i INPUT1 [INPUT2 ...] [-c SUBSTRATE] [--parm TOPOLOGY] [options]
```

`mlmm all --help` shows core options; `mlmm all --help-advanced` shows the full option list.

Multi-structure MEP with full post-processing:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft --out-dir ./result_all
```

Single-structure staged scan (two stages):

```bash
mlmm all -i A.pdb -c '308,309' --scan-lists '[(12,45,1.35)]' '[(10,55,2.20)]' \
    --multiplicity 1 --out-dir ./result_scan_all
# a single literal can drive several bonds at once: '[(10,55,2.20),(23,34,1.80)]'
```

TSOPT-only validation (single input, no MEP search):

```bash
mlmm all -i A.pdb -c 'GPP,MMT' -l 'GPP:-3,MMT:-1' \
    --tsopt --thermo --dft --out-dir result_tsopt_only
```

ORB backend:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --backend orb --out-dir ./result_all_orb
```

DMF with the CPU implementation (the default DMF backend is GPU):

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --mep-mode dmf --dmf-backend cpu --out-dir ./result_all_dmf
```

Converted trajectory/structure PDB companions are generated when reference
templates are available; control those conversions with `--convert-files` (on
by default). The inspectable ML-region PDB pair described below is a separate
artifact and is always written for PDB input.

## Workflow

1. **Active-site extraction and ML-region definition** (multi-structure union when multiple inputs)
   - Define the substrate via `-c/--center` (PDB path, residue IDs, or residue names) and optionally `--ligand-charge` as a total number (distributed) or a mapping such as `GPP:-3,MMT:-1`.
   - The extractor writes per-input pocket PDBs under `<out-dir>/_work/pockets/`. The first pocket is copied to `<out-dir>/ml_region.pdb` (a reusable deliverable you can pass back as `--model-pdb`) and defines the ML region for all subsequent ML/MM calculations.
   - `<out-dir>/ml_region_without_linkH.xyz` and `ml_region_with_linkH.xyz` expose the exact model system before and after link-H insertion. PDB inputs also produce matching `.pdb` companions. Automatic link pairs are parm7 bonds crossing the ML/MM selection, never distance-perceived bonds.
   - The **first-model net ML-region charge** becomes the net ML-region charge for later steps.
   - Omitting `-c/--center` skips extraction and uses the full input structures directly.
2. **ML/MM preparation (parm7 + layer assignment)**
   - `mm_parm` runs once on the first full input PDB and writes `<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7` (a reusable deliverable you can pass back as `--parm`), which are passed automatically as `--parm`.
   - `define-layer` runs on each full-system PDB and assigns 3-layer B-factors (ML = 0.0, Movable-MM = 10.0, Frozen = 20.0) based on the ML-region definition. The layered full-system PDBs are written under `<out-dir>/layered/`.
3. **Optional staged scan** (single-structure only)
   - When exactly one input PDB is provided and `--scan-lists` is given, the tool performs a staged, bond-length-driven scan on the layered full-system PDB using the ML/MM calculator.
   - Each stage's relaxed structure (`stage_XX/result.pdb`) is collected as an intermediate / product candidate. The ordered input series for the path search becomes `[initial layered PDB, stage_01/result.pdb, stage_02/result.pdb, ...]`.
4. **MEP search on full-system layered PDBs**
   - All MEP calculations run on full-system layered PDBs (with `--parm` and `--detect-layer`), not on pockets.
   - **`--refine-path`** runs recursive `path_search` with automatic refinement, detecting multistep reactions and building a detailed MEP per elementary step. Complex multistep mechanisms may need manual trial-and-error to obtain a converged pathway.
   - Select GSM (default) or DMF with `--mep-mode`. `--dmf-backend gpu` uses `dmf.torch`; use `--dmf-backend cpu` for the NumPy implementation or after a GPU out-of-memory error.
   - **`--no-refine-path` (default)** runs `path-opt` with the selected optimizer per adjacent pair, then concatenates trajectories, extracts the HEI per segment, detects bond changes, and writes `summary.json`. Both modes support Stage 5 post-processing.
   - For multi-input runs, the original full PDBs are supplied as merge references automatically. In the scan-derived series (single-structure case), the single original full PDB is reused as the reference template.
5. **Summary and optional post-processing**
   - The raw MEP-engine output (per-segment trajectories, the full MEP trajectory, and the engine `summary.json`) is written under `<out-dir>/_work/path_opt/` (or `<out-dir>/_work/path_search/` with `--refine-path`); the merged products (`mep.pdb`, optional `mep.cif`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`) are moved to `<out-dir>/` and `summary.{json,log}` copied there.
   - `--tsopt` runs TS optimization on each HEI. After the TS gate, `all` continues with EulerPC IRC and segment energy diagrams.
   - `--thermo` computes ML/MM thermochemistry on (R, TS, P) and adds a Gibbs diagram.
   - `--dft` runs model-region DFT single-points on (R, TS, P) and adds a model-DFT electronic diagram. With `--thermo`, the subtractive DFT//MLIP/MM total plus the ML/MM thermal correction produces the DFT//MLIP/MM Gibbs diagram.
   - TS optimization, IRC, frequency analysis, and flatten PHVA use the fixed constrained treatment, which removes only full-system rigid motions that leave frozen anchors fixed; realistic ML/MM boundaries normally have effective rank 0.
   - `--hessian-calc-mode` selects analytical or finite-difference Hessians where supported. Compare both on a target-system pilot because speed and memory depend on the backend and system.
6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists`)
   - Skips the MEP search and runs `tsopt` on the layered full-system PDB. After the TS gate, it performs EulerPC IRC, minimizes both ends, and optionally adds thermochemistry, DFT, and DFT//MLIP/MM diagrams.
   - When IRC runs, its ends are emitted as chemically unassigned `E1` and `E2` because no path/reference orientation is available. The summary reports the barrier from each endpoint to TS and does not emit R/P reaction energies. Inspect the structures before assigning chemical identities.

## Outputs

The tree has three zones: **deliverables at the root**, **per-segment deliverables under `segments/seg_NN/`**, and **pipeline scratch under `_work/`** (safe to remove once you have the results). The three you check first are `summary.log`, `summary.json`, and `mep.pdb` (the concatenated reaction path, moved to the root; raw engine output stays under `_work/path_opt/` by default, or `_work/path_search/` with `--refine-path`).

```text
<out-dir>/
  summary.json                   # mirrored top-level summary (when the MEP stage runs)
  summary.log
  mep.pdb · mep.cif             # path; CIF companion is emitted for bridged input
  mep_trj.xyz
  mep_plot.png                   # smooth MEP energy profile
  energy_diagram_MEP.png         # all-segment MEP barriers
  energy_diagram_MLIP_all.png           # aggregated post-processing diagrams (when enabled)
  energy_diagram_G_MLIP_all.png
  energy_diagram_DFT_all.png
  energy_diagram_G_DFT_plus_MLIP_all.png
  irc_plot_all.png
  ml_region.pdb                  # ML-region definition (reusable as --model-pdb for follow-up runs)
  ml_region_without_linkH.xyz    # exact ML model before link-H insertion
  ml_region_with_linkH.xyz       # exact ML model after parm7-derived link-H insertion
  ml_region_without_linkH.pdb    # topology-bearing companion for PDB input
  ml_region_with_linkH.pdb       # PDB companion with generated HL/LKH atoms
  mm_parm/<input1>.parm7,.rst7   # MM topology from the first full-enzyme input (reusable as --parm)
  layered/                       # Layered full-system PDBs (B-factor annotated; reusable inputs)
  segments/                      # per-reactive-segment deliverables
    seg_NN/                      # 1-based 2-digit index, e.g. seg_01, seg_02
      reactant.pdb · ts.pdb · product.pdb   # canonical R/TS/P for MEP runs
      e1.pdb · ts.pdb · e2.pdb              # unassigned endpoints for TSOPT-only
      *.cif                                 # bridged-input companions with original IDs
      ts/                        # TS optimization (--tsopt)
      irc/                       # EulerPC IRC after the TS gate
      freq/ (--thermo), dft/ (--dft)
      structures/{reactant,ts,product}.pdb  # MEP run nested copy
      structures/{endpoint_1,ts,endpoint_2}.pdb # TSOPT-only nested copy
      energy_diagram_{MLIP,G_MLIP,DFT,G_DFT_plus_MLIP}.png
  _work/                         # pipeline scratch (safe to delete)
    pockets/                     # Per-input pocket PDBs (multi-structure union)
    scan/                        # present only in single-structure + scan mode (stage_01/result.pdb …)
    path_opt/                    # raw MEP-engine output (path_search/ with --refine-path)
      summary.{json,log} · seg_NN_mep/    # raw per-segment MEP trajectories (merged products are moved to the root)
```

In **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) there is no MEP stage. `ts/` is written under `segments/seg_01/`; after the TS gate, the E1/TS/E2 structures and `irc/` are added there, followed by requested `freq/` and `dft/` outputs. `_work/path_opt/` is absent.

At `-v 2` the console summarises extraction, MM preparation, scan stages, MEP progress, and per-stage timing; see {ref}`verbosity-levels`.

### Reading `summary.log`

The header identifies the `all` entry route as `MEP`, `Scan`, or `TS-only` and
prints the absolute root output directory and, for `MEP` / `Scan`, the absolute
internal path module directory (`TS-only` reports `-`). Internal engine names
such as `path-opt` / `path-search`
remain machine-readable metadata rather than the user-facing pipeline mode.

The log is organized into numbered sections:

- **[1] Global MEP overview** — image / segment counts, MEP trajectory plot paths, aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (MLIP path)** — per-segment barriers, reaction energies, bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** — TS imaginary-frequency checks, gated IRC outputs, energy tables.
- **[4] Energy diagrams (overview)** — diagram tables for MEP / MLIP / Gibbs / DFT plus an optional cross-method summary.
- **[5] Output directory structure** — a compact tree of generated files with inline annotations.

### Reading `summary.json`

Top-level keys: `out_dir`, `n_images`, `n_segments` (run metadata and counts); `segments` (per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`); `energy_diagrams` (optional payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` paths).

When stage `result.json` files or `thermoanalysis.yaml` are written, their
`rigid_projection` block records the selected treatment, effective rank,
Hessian source, and Hessian shape. An all-frozen selection is rejected because
no active DOF remains.

## CLI options

Defaults shown are used when the option is not specified. The full flag list is in the generated [command reference](reference/commands/index.md); the tables below cover the options that need explanation.

### Input / output

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full structures in reaction order: PDB/mmCIF directly, or XYZ with `--ref-pdb` (single input allowed with `--scan-lists` or `--tsopt`). | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). Omit to skip extraction. | _None_ |
| `-l, --ligand-charge TEXT` | Total charge or residue-specific mapping (e.g. `GPP:-3,MMT:-1`). | _None_ |
| `-q, --charge INT` | Override the net charge of the ML region/model atoms (highest priority). | _None_ |
| `--freeze-atoms TEXT` | Comma-separated 1-based full-system atom indices frozen throughout scan/MEP/TSOPT/IRC/frequency stages. Merged with YAML `geom.freeze_atoms` and the automatically detected Frozen-MM layer. | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--parm FILE` | AMBER parm7 topology for the full (real) system. Auto-generated by `mm_parm` when omitted. | _None_ |
| `--model-pdb FILE` | Pre-built ML-region PDB. When provided, ML-region determination is skipped. | _None_ |
| `--ref-pdb FILE` | Reference PDB for XYZ input (required so PDB metadata can be recovered). | _None_ |
| `--convert-files / --no-convert-files` | Global toggle for XYZ / TRJ → PDB companions. | `True` |
| `--dump / --no-dump` | Save optional optimizer trajectories/restarts. An explicit parent toggle is forwarded to `path-search` / `path-opt` and `scan` / `tsopt`; when omitted, each child resolves its YAML/default. With `--thermo`, the required child `thermoanalysis.yaml` handoff is retained even under `--no-dump` so Gibbs assembly remains complete. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config / --no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run / --no-dry-run` | Run extraction/setup and charge/parity validation in a temporary directory, print the plan, and skip compute stages (shown in `--help-advanced`). | `False` |

### Extraction

| Option | Description | Default |
| --- | --- | --- |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). `0` is accepted and evaluated internally as `0.001 Å` (effectively off for ordinary radius neighbors). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff (Å). | `0.0` |
| `--include-h2o / --no-include-h2o` | Include water molecules (HOH / WAT / H2O / DOD / TIP / TIP3 / SOL). | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkh / --no-add-linkh` | Add link hydrogens for severed bonds. | `False` |
| `--selected-resn TEXT` | Force-include IDs/names such as `123`, `A:123A`, `SAM`, `A:SAM`, or `A:SAM:123` (comma/space separated). | `""` |
| `--modified-residue TEXT` | Comma-separated modified-residue names and integer charges for backbone truncation and charge assignment (e.g. `HD1:0,HD2:-1`). A known catalog residue may omit its charge (e.g. `SEP`). | `""` |

### MM preparation

| Option | Description | Default |
| --- | --- | --- |
| `--auto-mm-ff-set {ff19SB\|ff14SB}` | Force-field set for `mm_parm` (ff19SB → OPC3; ff14SB → TIP3P). | `ff19SB` |
| `--auto-mm-add-ter / --auto-mm-no-add-ter` | Control TER insertion around ligand / water / ion blocks. | `True` |
| `--auto-mm-disulfide / --auto-mm-no-disulfide` | Forwarded to mm_parm: detect disulfides from SG-SG geometry across CYS/CYM/CYX and bond them (renaming a bonded CYS to CYX). With `--auto-mm-no-disulfide` only residues already named CYX are bonded. | `True` |
| `--auto-mm-keep-temp` | Keep the `mm_parm` temporary working directory (for debugging). | `False` |
| `--auto-mm-ligand-mult TEXT` | Spin multiplicity mapping forwarded to `mm_parm` (e.g. `GPP:2,SAM:1`). If omitted, defaults to 1 for all ligands. | _None_ |

### MEP search

```{note}
`--max-cycles-gsm` and `--max-cycles-dmf` control only the selected MEP child
and each default to 300. Scan, TS optimization, IRC, and other stages keep
their dedicated cycle options and defaults.
```

| Option | Description | Default |
| --- | --- | --- |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--mep-mode [gsm\|dmf]` | MEP optimizer forwarded to both `path-opt` and recursive `path-search`. | `gsm` |
| `--dmf-backend [gpu\|cpu]` | DMF implementation. The parent forwards this only when explicitly set, so a child YAML `dmf.backend` remains effective otherwise. | `gpu` |
| `--max-nodes INT` | Internal nodes per GSM/DMF segment. | `20` |
| `--gsm-param [equi\|energy]` | GSM node parameterization after string growth. `energy` concentrates nodes in high-energy regions and may be tried when an equidistant path skips the reaction-coordinate region near the HEI; it does not identify a TS. | `equi` |
| `--max-cycles-gsm INT` | GSM string-optimizer cycle cap for the MEP child. | `300` |
| `--max-cycles-dmf INT` | DMF IPOPT iteration cap for the MEP child. | `3000` |
| `--climb / --no-climb` | Enable climbing-image TS refinement where supported by the selected optimizer. | `True` |
| `--opt-mode [grad\|hess]` | Fallback preset for TSOPT and post-IRC endpoint optimization (`grad` → Dimer / L-BFGS, `hess` → RS-P-RFO / RFO). `--opt-mode-post` takes precedence. | `grad` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset override for TSOPT / post-IRC endpoint optimizations (`grad` → Dimer / L-BFGS, `hess` → RS-P-RFO / RFO). | `hess` |
| `--thresh TEXT` | Convergence preset for single-structure optimizations and scan relaxations (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--thresh-gsm TEXT` | Convergence preset for the GSM string optimizer of the MEP stage (same presets as `--thresh`). | `gau_loose` |
| `--thresh-dmf TEXT` | IPOPT dual-infeasibility tolerance of the DMF MEP stage: `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive float. Not a Gaussian preset. | `tight` |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--preopt / --no-preopt` | Pre-optimize endpoints before segmentation. | `True` |
| `--refine-path / --no-refine-path` | `--no-refine-path` (default) → single-pass `path-opt`; `--refine-path` → recursive `path-search`, which discovers multistep mechanisms and also refines a single-step MEP, where it can improve a poor HEI or TS estimate. Both modes support Stage 5 (TSOPT / thermo / DFT). | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--precision [fp32\|fp64]` | Backend precision. Unset uses UMA/AIMNet2 fp32 and ORB/MACE fp64. AIMNet2 rejects fp64. | backend-specific |
| `--workers INT` | UMA predictor workers. Values greater than 1 require `fairchem-core[extras]` and are incompatible with an analytical Hessian. | `1` |
| `--workers-per-node INT` | Workers per node for the parallel UMA predictor. | _None_ |
| `--cmap / --no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
| `--hessian-calc-mode CHOICE` | ML/MM Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--detect-layer / --no-detect-layer` | Automatically read B-factor layers (B = 0 / 10 / 20). With explicit `--model-pdb`, retain only the MM sublayers; otherwise B-factors also define ML membership. | Enabled |

TSOPT optimizer selection order: `--opt-mode-post` (if set) → `--opt-mode` (only when explicitly provided) → TSOPT default (`hess` → RS-P-RFO).

### Scan (single-input runs)

| Option | Description | Default |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | Inline `(i, j, target_Å)` literals. One literal is one stage; several tuples within it form a concerted scan. Use standalone `scan` for YAML/JSON or bidirectional 4-tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | `<out-dir>/_work/scan` |
| `--scan-one-based / --scan-zero-based` | Interpret scan atom indices as 1-based or 0-based. | _None_ |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | `0.20` |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV / Å²). | `300.0` |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | `100000` |
| `--scan-preopt / --no-scan-preopt` | Override scan pre-optimization toggle. | _None_ |
| `--scan-endopt / --no-scan-endopt` | Override scan end-of-stage optimization. | _None_ |

### Post-processing + freq / DFT overrides

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt / --no-tsopt` | Run TS optimization and, after the TS gate, EulerPC IRC per reactive segment. | `False` |
| `--tsopt-from-mep-tan / --no-tsopt-from-mep-tan` | For Hessian TS optimizers, guide reaction-root identity with CPU/file-cached HEI tangent candidates. Turning it off disables cache creation/use and selects from initial Hessian modes. Not applicable to Dimer. | `True` |
| `--thermo / --no-thermo` | Run vibrational analysis (`freq`) on R/TS/P for MEP runs or E1/TS/E2 for TS-only runs. | `False` |
| `--dft / --no-dft` | Run single-point DFT on R/TS/P for MEP runs or E1/TS/E2 for TS-only runs. | `False` |
| `--flatten / --no-flatten` | Surplus-imaginary-mode flattening in `tsopt`. | `False` |
| `--reject-uphill / --no-reject-uphill` | Opt in to rejecting energy-raising RFO steps during post-IRC **endpoint re-optimization only**, using a `1e-4` Hartree tolerance (forwarded to the opt child); TS optimization forces rejection off, and path search is unaffected. At the emergency floor, the retained endpoint receives a final normal convergence check. | `False` |
| `--irc-step-size FLOAT` | Override the EulerPC maximum step (Bohr) for every post-TS IRC. If a branch stops after only a few frames, retry with a smaller value such as `0.05`. | IRC default `0.10` |
| `--irc-never-stop / --no-irc-never-stop` | Ignore IRC gradient and energy endpoint criteria and trace each branch to the cycle cap. Numerical/integration failures and external interruption still stop propagation. | `False` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | `100000` |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | `10` |
| `--freq-amplitude-ang FLOAT` | Mode-trajectory amplitude (Å). | `0.8` |
| `--freq-n-frames INT` | Frames per mode trajectory. | `20` |
| `--freq-sort TEXT` | Mode sorting behavior. | `value` |
| `--freq-temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--freq-pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Functional / basis pair. | `wb97m-v/def2-tzvpd` |
| `--dft-max-cycle INT` | SCF-iteration cap. | `100` |
| `--dft-conv-tol FLOAT` | SCF convergence tolerance. | `1e-9` |
| `--dft-grid-level INT` | PySCF grid level. | `3` |
| `--dft-engine [gpu\|cpu]` | DFT engine (GPU or CPU PySCF). | `gpu` |

## YAML configuration

`all` accepts `--config FILE` with the public precedence `defaults < config < explicit CLI`. The effective YAML is forwarded to downstream subcommands, and each tool reads the sections described in its own documentation:

| Subcommand | YAML sections |
|---|---|
| [`path-search`](path-search.md) | `geom`, `calc` / `mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc` / `mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc` / `mlmm`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc` / `mlmm`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

```yaml
# Minimal example
geom:
  tr_projection: constrained        # fixed internal PHVA treatment
calc:
  model_charge: 0
  model_mult: 1
  backend: uma                      # uma | orb | mace | aimnet2
  uma_model: uma-s-1p2              # uma-s-1p2 | uma-m-1p1
  hessian_calc_mode: Analytical     # compare with FiniteDifference on a pilot
gs:
  max_nodes: 20
  climb: true
dft:
  grid_level: 6
```

Full schema: [YAML Reference](yaml-reference.md).

The fixed constrained rigid-mode treatment is unrelated to `tsopt --ref-mode`,
which is an internal MEP-tangent handoff for TS root selection and overlap
tracking. A stale
non-constrained `geom.tr_projection` value fails explicitly.

## Notes

Input format depends on extraction:

- PDB and mmCIF inputs are accepted directly.
- XYZ inputs require `--ref-pdb`; XYZ supplies coordinates and the reference
  supplies residue, chain, and B-factor metadata for extraction and later stages.
- Multi-structure runs require ≥ 2 structures.

Charge priority is explicit `-q/--charge` → workflow-derived charge (extraction, or selected-model derivation using `--ligand-charge`) → YAML `calc.model_charge` → error. Multiplicity priority is explicit `--multiplicity` → YAML `calc.model_mult` → 1. Provide `--ligand-charge` for non-standard substrates. The first-model net ML-region charge is rounded to the nearest integer with a console note.

## See Also

[extract](extract.md) (called internally by `all`) · [mm-parm](mm-parm.md) (called internally by `all`) · [path-search](path-search.md) · [tsopt](tsopt.md) · [freq](freq.md) · [dft](dft.md) · [trj2fig](trj2fig.md) · [Common Error Recipes](recipes-common-errors.md) (symptom-first failure routing) · [Troubleshooting](troubleshooting.md) (common errors and fixes) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
