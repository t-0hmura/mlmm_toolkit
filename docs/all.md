# `all`

## Overview

> **Summary:** End-to-end enzymatic reaction workflow -- active-site extraction, ML/MM layer assignment, MM topology preparation, optional staged scan, MEP search (GSM) on full-system layered PDBs, with optional TS optimization, EulerPC IRC, thermochemistry, DFT, and DFT//MLIP diagrams.

`mlmm all` runs a one-shot pipeline that operates on full-system layered PDBs with ML/MM. It supports three modes:

- **Multi-structure ensemble** -- Provide two or more full PDBs in reaction order. The tool extracts the active-site region (for ML-region definition), builds MM topology, assigns ML/MM layers, runs GSM MEP search on the layered full-system PDBs, and optionally runs per-segment post-processing (TSOPT/freq/DFT).
- **Single-structure + staged scan** -- Provide one PDB plus `--scan-lists`. The scan generates intermediate/product candidates that become MEP endpoints.
  - One `--scan-lists` literal runs a single scan stage.
  - Multiple stages are passed as multiple values after a single `--scan-lists` flag.
- **TSOPT-only** -- Provide a single PDB, omit `--scan-lists`, and set `--tsopt`. The tool runs TS optimization on the layered full-system PDB, performs EulerPC IRC, minimizes both ends, and builds energy diagrams.

```{important}
`--tsopt` produces **TS candidates**. `all` automatically runs IRC and freq for validation, but always inspect the results (imaginary mode + endpoint connectivity) before mechanistic interpretation.
```

## Minimal example

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" --out-dir ./result_all
```

## Output checklist

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or `result_all/path_search/seg_*/`)

## Common examples

1. Run full post-processing in one command.

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

2. Single-structure staged scan route.

```bash
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
```

3. Validate parsing and plan only.

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" --dry-run
```

4. Use the ORB backend with xTB point-charge embedding.

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" -l "SAM:1,GPP:-3" \
 --backend orb --embedcharge --out-dir ./result_all_orb
```

PDB companions are generated when templates are available, controlled by `--convert-files/--no-convert-files` (enabled by default).

## Usage

```bash
mlmm all -i INPUT1 [INPUT2...] -c SUBSTRATE [options]
```

Run `mlmm all --help` for core options or `mlmm all --help-advanced` for the full option list.

### Examples

```bash
# Minimal end-to-end run with explicit substrate and ligand charges (multi-structure)
mlmm all -i reactant.pdb product.pdb -c "GPP,MMT" -l "GPP:-3,MMT:-1"

# Full ensemble with an intermediate, residue-ID substrate spec, and full post-processing
mlmm all -i A.pdb B.pdb C.pdb -c "308,309" -l "-1" \
 --multiplicity 1 --max-nodes 10 --max-cycles 100 --climb \
 --opt-mode grad --no-dump --config params.yaml --preopt \
 --out-dir result_all --tsopt --thermo --dft

# Single-structure + scan to build an ordered series
mlmm all -i A.pdb -c "308,309" --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
 --multiplicity 1 --out-dir result_scan_all --tsopt --thermo --dft

# Single-structure TSOPT-only mode (no path_search)
mlmm all -i A.pdb -c "GPP,MMT" -l "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_only
```

## Workflow

1. **Active-site extraction and ML-region definition** (multi-structure union when multiple inputs)
   - Define the substrate (`-c/--center`, by PDB, residue IDs, or residue names).
   - Optionally provide `--ligand-charge` as a total number (distributed) or a mapping (e.g., `GPP:-3,MMT:-1`).
   - The extractor writes per-input pocket PDBs under `<out-dir>/pockets/`. The first pocket is copied as `<out-dir>/ml_region.pdb` to define the ML region for all subsequent ML/MM calculations.
   - The extractor's **first-model total pocket charge** is used as the total charge in later steps, cast to the nearest integer with a console note if rounding occurs.
   - Additional extractor toggles: `--radius`, `--radius-het2het`, `--include-H2O/--no-include-H2O`, `--exclude-backbone/--no-exclude-backbone`, `--add-linkH/--no-add-linkH`, `--selected-resn`, `--verbose/--no-verbose`.
   - If `-c/--center` is omitted, extraction is skipped and full input structures are used directly.

2. **ML/MM preparation (parm7 + layer assignment)**
   - Run `mm_parm` once on the first full input PDB to build `<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7`, automatically passed as `--parm`.
   - Run `define-layer` on each full-system PDB to assign 3-layer B-factors (ML=0.0, MovableMM=10.0, FrozenMM=20.0) based on the ML-region definition. The layered full-system PDBs are written under `<out-dir>/layered/`.
   - Tune this stage with `--auto-mm-ff-set`, `--auto-mm-add-ter`, and `--auto-mm-keep-temp`.

3. **Optional staged scan (single-structure only)**
   - If exactly one full input PDB is provided and `--scan-lists` is given, the tool performs a staged, bond-length-driven scan on the layered full-system PDB using the ML/MM calculator.
   - For each stage, the final relaxed structure (`stage_XX/result.pdb`) is collected as an intermediate/product candidate.
   - The ordered input series for the path search becomes: `[initial layered PDB, stage_01/result.pdb, stage_02/result.pdb,...]`.

4. **MEP search on full-system layered PDBs**
   - All MEP calculations run on full-system layered PDBs (with `--parm` and `--detect-layer`), not on pockets.
   - **`--refine-path` (default):** Runs recursive `path_search` with kink-detection and refinement.
   - **`--no-refine-path`:** Runs single-pass `path-opt` GSM per adjacent pair, then concatenates trajectories, extracts HEI per segment, detects bond changes, and writes `summary.yaml` — enabling Stage 4 post-processing (TSOPT, thermo, DFT) on both modes.
   - For multi-input runs, the original full PDBs are supplied as merge references automatically. In the scan-derived series (single-structure case), the single original full PDB is reused (repeated) as the reference template.

5. **Summary and optional post-processing**
   - Per-segment trajectories, full MEP trajectory, and a `summary.yaml` are written under `<out-dir>/path_search/`.
   - `--tsopt`: run TS optimization on each HEI, follow with EulerPC IRC, and emit segment energy diagrams.
   - `--thermo`: Compute ML/MM thermochemistry on (R, TS, P) and add a Gibbs diagram.
   - `--dft`: Do DFT single-point on (R, TS, P) and add a DFT diagram. With `--thermo`, also generate a DFT//MLIP Gibbs diagram.
   - Shared overrides include `--opt-mode`, `--opt-mode-post` (overrides TSOPT and post-IRC endpoint optimization modes), `--flatten/--no-flatten`, `--hessian-calc-mode`, `--tsopt-max-cycles`, `--tsopt-out-dir`, `--freq-*`, `--dft-*`.
   - When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.

6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists`)
   - Skips steps (4)-(5) and runs `tsopt` on the layered full-system PDB, performs EulerPC IRC and minimizes both ends, builds ML/MM energy diagrams for R-TS-P, and optionally adds Gibbs, DFT, and DFT//MLIP diagrams.
   - In this mode only, the IRC endpoint with **higher energy** is adopted as the reactant (R).

### Charge and spin precedence

**Charge resolution (highest to lowest priority):**

| Priority | Source | When Used |
|----------|--------|-----------|
| 1 | `-q/--charge` | Explicit CLI override |
| 2 | Pocket extraction | When `-c` is provided (sums amino acids, ions, `--ligand-charge`) |
| 3 | `-l, --ligand-charge` | Fallback when extraction is skipped |
| 4 | Default | None (unresolved charge is an error) |

**Spin resolution:** `--multiplicity` (CLI) -> default (1)

> **Tip:** Always provide `--ligand-charge` for non-standard substrates to ensure correct charge propagation.

### Input expectations
- Extraction enabled (`-c/--center`): inputs must be **PDB** files so residues can be located.
- Extraction skipped: inputs may be **PDB/XYZ**.
- Multi-structure runs require at least 2 structures.

## CLI options

> **Note:** Default values shown are used when the option is not specified.

### Input/Output Options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order (single input allowed with `--scan-lists` or `--tsopt`). | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). Omit to skip extraction. | _None_ |
| `-l, --ligand-charge TEXT` | Total charge or residue-specific mapping (e.g., `GPP:-3,MMT:-1`). | _None_ |
| `-q, --charge INT` | Force total system charge (highest priority override). | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--parm FILE` | AMBER parm7 topology file for the full (real) system. When omitted, `all` auto-generates one via `mm_parm`. | _None_ |
| `--model-pdb FILE` | Pre-built ML-region PDB. When provided, pocket extraction is skipped and this file defines the ML region directly. | _None_ |
| `--ref-pdb FILE` | Reference PDB for XYZ input. Required when the input is XYZ so that PDB metadata (residues, chains, B-factors) can be recovered. | _None_ |
| `--convert-files/--no-convert-files` | Global toggle for XYZ/TRJ to PDB companions when templates are available. | `True` |
| `--dump/--no-dump` | Save optimizer dumps. Always forwarded to `path-search`/`path-opt`; forwarded to `scan`/`tsopt` only when explicitly set here. `freq` defaults to dump=True unless you pass `--no-dump`. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run/--no-dry-run` | Validate and print plan without running stages. Shown in `--help-advanced`. | `False` |

### Extraction Options

| Option | Description | Default |
| --- | --- | --- |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff (Å). | `0.0` |
| `--include-H2O/--no-include-H2O` | Include water molecules (HOH/WAT/TIP3/SOL). | `True` |
| `--exclude-backbone/--no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `True` |
| `--add-linkH/--no-add-linkH` | Add link hydrogens for severed bonds. | `False` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--verbose/--no-verbose` | Enable INFO-level extractor logging. | `True` |

### MM Preparation Options

| Option | Description | Default |
| --- | --- | --- |
| `--auto-mm-ff-set {ff19SB\|ff14SB}` | Force field set for `mm_parm` (ff19SB uses OPC3; ff14SB uses TIP3P). | `ff19SB` |
| `--auto-mm-add-ter/--auto-mm-no-add-ter` | Control TER insertion around ligand/water/ion blocks. | `True` |
| `--auto-mm-keep-temp` | Keep the `mm_parm` temporary working directory (for debugging). | `False` |
| `--auto-mm-ligand-mult TEXT` | Spin multiplicity mapping forwarded to `mm_parm` (e.g., `GPP:2,SAM:1`). If omitted, `mm_parm` defaults to 1 for all ligands. | _None_ |
### MEP Search Options

| Option | Description | Default |
| --- | --- | --- |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--max-nodes INT` | Internal nodes for segment GSM. | `20` |
| `--max-cycles INT` | Max GSM macro-cycles. | `300` |
| `--climb/--no-climb` | Enable TS refinement for segment GSM. | `True` |
| `--opt-mode [grad\|hess]` | Optimizer preset for scan/path-search and single optimizations (`grad` -> LBFGS/Dimer, `hess` -> RFO/RSIRFO). | `grad` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset override for TSOPT/post-IRC endpoint optimizations (`grad` -> Dimer/LBFGS, `hess` -> RS-I-RFO/RFO). | `hess` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--preopt/--no-preopt` | Pre-optimize endpoints before segmentation. | `False` |
| `--refine-path/--no-refine-path` | If True, run recursive `path-search`; if False, chain `path-opt` segments (single-pass GSM per pair, with trajectory concatenation, HEI extraction, bond-change detection, and summary.yaml). Both modes support Stage 4 (TSOPT/thermo/DFT). | `True` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap/--no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in model parm7. Default: disabled (consistent with Gaussian ONIOM). | `--no-cmap` |
| `--hessian-calc-mode CHOICE` | ML/MM Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--detect-layer/--no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20) in downstream tools. If disabled, downstream tools require `--model-pdb` or `--model-indices`. | `True` |

TSOPT optimizer selection order: `--opt-mode-post` (if set) -> `--opt-mode` (only when explicitly provided) -> TSOPT default (`hess` -> `heavy`).

### Scan Options (Single-Input Runs)

| Option | Description | Default |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | Staged scans: `(i, j, target_Å)` tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based/--no-scan-one-based` | Override scan indexing interpretation (True = 1-based, False = 0-based). | _None_ |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | _Default_ |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV/Å²). | _Default_ |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | _Default_ |
| `--scan-preopt/--no-scan-preopt` | Override scan pre-optimization toggle. | _None_ |
| `--scan-endopt/--no-scan-endopt` | Override scan end-of-stage optimization. | _None_ |

### Post-Processing Options

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | Run TS optimization + EulerPC IRC per reactive segment. | `False` |
| `--thermo/--no-thermo` | Run vibrational analysis (`freq`) on R/TS/P. | `False` |
| `--dft/--no-dft` | Run single-point DFT on R/TS/P. | `False` |
| `--flatten/--no-flatten` | Enable extra-imaginary-mode flattening in `tsopt`. | `False` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | _Default_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |

### Freq Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | _Default_ |
| `--freq-amplitude-ang FLOAT` | Mode animation amplitude (Å). | _Default_ |
| `--freq-n-frames INT` | Frames per mode animation. | _Default_ |
| `--freq-sort TEXT` | Mode sorting behavior. | _Default_ |
| `--freq-temperature FLOAT` | Thermochemistry temperature (K). | _Default_ |
| `--freq-pressure FLOAT` | Thermochemistry pressure (atm). | _Default_ |

### DFT Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Functional/basis pair. | _Default_ |
| `--dft-max-cycle INT` | Maximum SCF iterations. | _Default_ |
| `--dft-conv-tol FLOAT` | SCF convergence tolerance. | _Default_ |
| `--dft-grid-level INT` | PySCF grid level. | _Default_ |

## Outputs

```text
<out-dir>/
 ml_region.pdb                         # ML-region definition (copy of the first pocket)
 pockets/
  pocket_<input1_basename>.pdb
  pocket_<input2_basename>.pdb
  ...
 mm_parm/
  <input1_basename>.parm7              # Generated from the first full-enzyme input PDB
  <input1_basename>.rst7
 layered/                              # Layered full-system PDBs (B-factor annotated)
 scan/                                 # present only in single-structure+scan mode
  stage_01/result.pdb
  stage_02/result.pdb
  ...
 summary.yaml                          # mirrored top-level summary (when path_search runs)
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 energy_diagram_UMA_all.png            # aggregated post-processing diagrams (when enabled)
 energy_diagram_G_UMA_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_UMA_all.png
 irc_plot_all.png
 path_search/                          # present when path_search is executed
  mep_trj.xyz
  mep.pdb
  summary.yaml
  summary.log
  mep_plot.png
  energy_diagram_MEP.png
  post_seg_XX/                         # when post-processing is enabled
   ts/...
   irc/...
   freq/...                            # with --thermo
   dft/...                             # with --dft
   energy_diagram_UMA.png
   energy_diagram_G_UMA.png
   energy_diagram_DFT.png
   energy_diagram_G_DFT_plus_UMA.png
 tsopt_single/                         # present only in single-structure TSOPT-only mode
  ts/...
  irc/...
  structures/
   reactant.pdb
   ts.pdb
   product.pdb
  freq/...                             # with --thermo
  dft/...                              # with --dft
  energy_diagram_UMA.png
  energy_diagram_G_UMA.png
  energy_diagram_DFT.png
  energy_diagram_G_DFT_plus_UMA.png
```

### Reading `summary.log`
The log is organized into numbered sections:
- **[1] Global MEP overview** -- image/segment counts, MEP trajectory plot paths, and the aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (MLIP path)** -- per-segment barriers, reaction energies, and bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** -- per-segment TS imaginary frequency checks, IRC outputs, and energy tables.
- **[4] Energy diagrams (overview)** -- diagram tables for MEP/MLIP/Gibbs/DFT series plus an optional cross-method summary table.
- **[5] Output directory structure** -- a compact tree of generated files with inline annotations.

### Reading `summary.yaml`
The YAML is a compact, machine-readable summary. Common top-level keys include:
- `out_dir`, `n_images`, `n_segments` -- run metadata and total counts.
- `segments` -- list of per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, and `bond_changes`.
- `energy_diagrams` (optional) -- diagram payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, and `image` paths.

## YAML configuration

`all` supports layered YAML:

- `--config FILE`: base settings.

`defaults < config < CLI < override-yaml`

The resulting effective YAML is forwarded to downstream subcommands. Each tool reads the sections described in its own documentation:

| Subcommand | YAML Sections |
|------------|---------------|
| [`path-search`](path_search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`/`mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc`/`mlmm`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc`/`mlmm`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

> **Note:** Applied after CLI values.

**Minimal example:**
```yaml
calc:
 charge: 0
 spin: 1
mlmm:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
 backend: uma                    # MLIP backend: uma | orb | mace | aimnet2
 embedcharge: false              # xTB point-charge embedding correction
 uma_model: uma-s-1p1            # uma-s-1p1 | uma-m-1p1
 hessian_calc_mode: Analytical     # recommended when VRAM permits
gs:
 max_nodes: 12
 climb: true
dft:
 grid_level: 6
```

For a complete reference of all YAML options, see **[YAML Configuration Reference](yaml_reference.md)**.

## See Also

- [extract](extract.md) -- Standalone pocket extraction (called internally by `all`)
- [mm_parm](mm_parm.md) -- Build AMBER topology (called internally by `all`)
- [path-search](path_search.md) -- Standalone recursive MEP search
- [tsopt](tsopt.md) -- Standalone TS optimization
- [freq](freq.md) -- Vibrational analysis and thermochemistry
- [dft](dft.md) -- Single-point DFT calculations
- [trj2fig](trj2fig.md) -- Plot energy profiles from trajectories
- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Common errors and fixes
- [YAML Reference](yaml_reference.md) -- Complete YAML configuration options
- [Glossary](glossary.md) -- Definitions of MEP, TS, IRC, GSM
