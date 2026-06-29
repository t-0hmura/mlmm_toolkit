# `all`

`mlmm all` runs the end-to-end ML/MM enzymatic-reaction workflow on full-system layered PDBs in one command, instead of chaining `extract` → `mm-parm` → `define-layer` → `scan` / `path-search` → `tsopt` → `irc` / `freq` / `dft` by hand. It chains active-site extraction, MM topology preparation, ML/MM layer assignment, an optional staged scan, MEP search (single-pass `path-opt` by default; recursive `path-search` with `--refine-path`), and optional post-processing (TS optimization, EulerPC IRC, thermochemistry, single-point DFT, and DFT//MLIP diagrams). The default MLIP backend for the ML region is UMA; choose an alternative with `-b/--backend`.

`all` runs in one of three modes, chosen by what you pass:

- **Multi-structure ensemble** — give ≥ 2 full PDBs in reaction order to drive a GSM MEP search across the supplied structures.
- **Single-structure staged scan** — give one PDB plus `--scan-lists`; each literal is a scan stage and the relaxed endpoints become the MEP endpoints.
- **TSOPT-only** — give a single PDB and set `--tsopt` (no `--scan-lists`) to run TS optimization directly, with no MEP search.

```{important}
`--tsopt` produces **TS candidates**. `all` runs IRC and freq automatically for validation, but always inspect the results (imaginary mode count + endpoint connectivity) before mechanistic interpretation.
```

## Examples

Command form:

```bash
mlmm all -i INPUT1 [INPUT2 ...] -c SUBSTRATE [options]
```

`mlmm all --help` shows core options; `mlmm all --help-advanced` shows the full option list.

Multi-structure MEP with full post-processing:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --tsopt True --thermo True --dft True --out-dir ./result_all
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
    --tsopt True --thermo True --dft True --out-dir result_tsopt_only
```

ORB backend with xTB point-charge embedding:

```bash
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' \
    --backend orb --embedcharge --out-dir ./result_all_orb
```

PDB companion files are generated when reference templates are available; control with `--convert-files` (on by default).

## Workflow

1. **Active-site extraction and ML-region definition** (multi-structure union when multiple inputs)
   - Define the substrate via `-c/--center` (PDB path, residue IDs, or residue names) and optionally `--ligand-charge` as a total number (distributed) or a mapping such as `GPP:-3,MMT:-1`.
   - The extractor writes per-input pocket PDBs under `<out-dir>/_work/pockets/`. The first pocket is copied to `<out-dir>/ml_region.pdb` (a reusable deliverable you can pass back as `--model-pdb`) and defines the ML region for all subsequent ML/MM calculations.
   - The **first-model net ML-region charge** becomes the net ML-region charge for later steps.
   - Omitting `-c/--center` skips extraction and uses the full input structures directly.
2. **ML/MM preparation (parm7 + layer assignment)**
   - `mm_parm` runs once on the first full input PDB and writes `<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7` (a reusable deliverable you can pass back as `--parm`), which are passed automatically as `--parm`.
   - `define-layer` runs on each full-system PDB and assigns 3-layer B-factors (ML = 0.0, Movable-MM = 10.0, Frozen-MM = 20.0) based on the ML-region definition. The layered full-system PDBs are written under `<out-dir>/layered/`.
3. **Optional staged scan** (single-structure only)
   - When exactly one input PDB is provided and `--scan-lists` is given, the tool performs a staged, bond-length-driven scan on the layered full-system PDB using the ML/MM calculator.
   - Each stage's relaxed structure (`stage_XX/result.pdb`) is collected as an intermediate / product candidate. The ordered input series for the path search becomes `[initial layered PDB, stage_01/result.pdb, stage_02/result.pdb, ...]`.
4. **MEP search on full-system layered PDBs**
   - All MEP calculations run on full-system layered PDBs (with `--parm` and `--detect-layer`), not on pockets.
   - **`--refine-path`** runs recursive `path_search` with automatic refinement, detecting multistep reactions and building a detailed MEP per elementary step. Complex multistep mechanisms may need manual trial-and-error to obtain a converged pathway.
   - **`--no-refine-path` (default)** runs `path-opt` GSM per adjacent pair, then concatenates trajectories, extracts the HEI per segment, detects bond changes, and writes `summary.json`. Both modes support Stage 5 post-processing.
   - For multi-input runs, the original full PDBs are supplied as merge references automatically. In the scan-derived series (single-structure case), the single original full PDB is reused as the reference template.
5. **Summary and optional post-processing**
   - The raw MEP-engine output (per-segment trajectories, the full MEP trajectory, and the engine `summary.json`) is written under `<out-dir>/_work/path_opt/` (or `<out-dir>/_work/path_search/` with `--refine-path`); the merged products (`mep.pdb`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`) are moved to `<out-dir>/` and `summary.{json,log}` copied there.
   - `--tsopt` runs TS optimization on each HEI, follows with EulerPC IRC, and emits segment energy diagrams.
   - `--thermo` computes ML/MM thermochemistry on (R, TS, P) and adds a Gibbs diagram.
   - `--dft` runs DFT single-point on (R, TS, P) and adds a DFT diagram. With `--thermo`, a DFT//MLIP Gibbs diagram is also produced.
   - When VRAM allows, set `--hessian-calc-mode Analytical` (strongly recommended over the FiniteDifference default).
6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists`)
   - Skips steps 4–5 and runs `tsopt` on the layered full-system PDB, performs EulerPC IRC, minimizes both ends, builds ML/MM energy diagrams for R-TS-P, and optionally adds Gibbs, DFT, and DFT//MLIP diagrams.
   - In this mode only, the IRC endpoint with **higher energy** is adopted as the reactant (R).

## Outputs

The tree has three zones: **deliverables at the root**, **per-segment deliverables under `segments/seg_NN/`**, and **pipeline scratch under `_work/`** (safe to remove once you have the results). The three you check first are `summary.log`, `summary.json`, and `mep.pdb` (the concatenated reaction path, moved to the root; raw engine output stays under `_work/path_opt/` by default, or `_work/path_search/` with `--refine-path`).

```text
<out-dir>/
  summary.json                   # mirrored top-level summary (when the MEP stage runs)
  summary.log
  mep.pdb                        # concatenated MEP path (copied to the root)
  mep_trj.xyz
  mep_plot.png                   # smooth MEP energy profile
  energy_diagram_MEP.png         # all-segment MEP barriers
  energy_diagram_UMA_all.png            # aggregated post-processing diagrams (when enabled)
  energy_diagram_G_UMA_all.png
  energy_diagram_DFT_all.png
  energy_diagram_G_DFT_plus_UMA_all.png
  irc_plot_all.png
  ml_region.pdb                  # ML-region definition (reusable as --model-pdb for follow-up runs)
  mm_parm/<input1>.parm7,.rst7   # MM topology from the first full-enzyme input (reusable as --parm)
  layered/                       # Layered full-system PDBs (B-factor annotated; reusable inputs)
  segments/                      # per-reactive-segment deliverables
    seg_NN/                      # 1-based 2-digit index, e.g. seg_01, seg_02
      reactant.pdb · ts.pdb · product.pdb   # canonical R/TS/P
      ts/, irc/                  # TS optimisation + EulerPC IRC (--tsopt)
      freq/ (--thermo), dft/ (--dft)
      structures/{reactant,ts,product}.pdb  # nested copy + raw IRC endpoints
      energy_diagram_{UMA,G_UMA,DFT,G_DFT_plus_UMA}.png
  _work/                         # pipeline scratch (safe to delete)
    pockets/                     # Per-input pocket PDBs (multi-structure union)
    scan/                        # present only in single-structure + scan mode (stage_01/result.pdb …)
    path_opt/                    # raw MEP-engine output (path_search/ with --refine-path)
      summary.{json,log} · seg_NN_mep/    # raw per-segment GSM trajectories (merged MEP products are moved to the root)
```

In **TSOPT-only mode** (single input + `--tsopt`, no `--scan-lists`) there is no MEP stage: the optimized R/TS/P plus `ts/`, `irc/`, `freq/`, and `dft/` land under `segments/seg_01/`, and `_work/path_opt/` is absent.

At `-v 2` the console summarises extraction, MM preparation, scan stages, MEP progress (GSM), and per-stage timing; see {ref}`verbosity-levels`.

### Reading `summary.log`

The log is organised into numbered sections:

- **[1] Global MEP overview** — image / segment counts, MEP trajectory plot paths, aggregate MEP energy diagram.
- **[2] Segment-level MEP summary (MLIP path)** — per-segment barriers, reaction energies, bond-change summaries.
- **[3] Per-segment post-processing (TSOPT / Thermo / DFT)** — TS imaginary-frequency checks, IRC outputs, energy tables.
- **[4] Energy diagrams (overview)** — diagram tables for MEP / MLIP / Gibbs / DFT plus an optional cross-method summary.
- **[5] Output directory structure** — a compact tree of generated files with inline annotations.

### Reading `summary.json`

Top-level keys: `out_dir`, `n_images`, `n_segments` (run metadata and counts); `segments` (per-segment entries with `index`, `tag`, `kind`, `barrier_kcal`, `delta_kcal`, `bond_changes`); `energy_diagrams` (optional payloads with `labels`, `energies_kcal`, `energies_au`, `ylabel`, `image` paths).

## CLI options

Defaults shown are used when the option is not specified. The full flag list is in the generated [command reference](reference/commands/index.md); the tables below cover the options that need explanation.

### Input / output

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order (single input allowed with `--scan-lists` or `--tsopt`). | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). Omit to skip extraction. | _None_ |
| `-l, --ligand-charge TEXT` | Total charge or residue-specific mapping (e.g. `GPP:-3,MMT:-1`). | _None_ |
| `-q, --charge INT` | Force net system charge (highest-priority override). | _None_ |
| `-o, --out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--parm FILE` | AMBER parm7 topology for the full (real) system. Auto-generated by `mm_parm` when omitted. | _None_ |
| `--model-pdb FILE` | Pre-built ML-region PDB. When provided, ML-region determination is skipped. | _None_ |
| `--ref-pdb FILE` | Reference PDB for XYZ input (required so PDB metadata can be recovered). | _None_ |
| `--convert-files / --no-convert-files` | Global toggle for XYZ / TRJ → PDB companions. | `True` |
| `--dump / --no-dump` | Save optimizer dumps. Always forwarded to `path-search` / `path-opt`; forwarded to `scan` / `tsopt` only when explicitly set. `freq` defaults to `dump=True` unless you pass `--no-dump`. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config / --no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run / --no-dry-run` | Validate and print plan without running stages (shown in `--help-advanced`). | `False` |

### Extraction

| Option | Description | Default |
| --- | --- | --- |
| `-r, --radius FLOAT` | Pocket inclusion cutoff (Å). | `2.6` |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff (Å). | `0.0` |
| `--include-h2o / --no-include-h2o` | Include water molecules (HOH / WAT / H2O / DOD / TIP / TIP3 / SOL). | `True` |
| `--exclude-backbone / --no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | `False` |
| `--add-linkh / --no-add-linkh` | Add link hydrogens for severed bonds. | `False` |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--modified-residue TEXT` | Comma-separated residue names (with optional charge) to treat as amino acids for backbone truncation and charge assignment (e.g. `HD1,HD2,HD3` or `HD1:0,SEP:-2`). | `""` |

### MM preparation

| Option | Description | Default |
| --- | --- | --- |
| `--auto-mm-ff-set {ff19SB\|ff14SB}` | Force-field set for `mm_parm` (ff19SB → OPC3; ff14SB → TIP3P). | `ff19SB` |
| `--auto-mm-add-ter / --auto-mm-no-add-ter` | Control TER insertion around ligand / water / ion blocks. | `True` |
| `--auto-mm-keep-temp` | Keep the `mm_parm` temporary working directory (for debugging). | `False` |
| `--auto-mm-ligand-mult TEXT` | Spin multiplicity mapping forwarded to `mm_parm` (e.g. `GPP:2,SAM:1`). If omitted, defaults to 1 for all ligands. | _None_ |

### MEP search

| Option | Description | Default |
| --- | --- | --- |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--max-nodes INT` | Internal nodes for segment GSM. | `20` |
| `--max-cycles INT` | Maximum GSM macro-cycles. | `300` |
| `--climb / --no-climb` | Enable TS refinement for segment GSM. | `True` |
| `--opt-mode [grad\|hess]` | Optimizer preset for scan / path-search and single optimizations (`grad` → LBFGS / Dimer, `hess` → RFO / RSIRFO). | `grad` |
| `--opt-mode-post [grad\|hess]` | Optimizer preset override for TSOPT / post-IRC endpoint optimizations (`grad` → Dimer / LBFGS, `hess` → RS-I-RFO / RFO). | `hess` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). Effective default: `gau_loose` for path-opt, `gau` for scan. | _None_ |
| `--thresh-post TEXT` | Convergence preset for post-IRC endpoint optimizations. | `baker` |
| `--preopt / --no-preopt` | Pre-optimize endpoints before segmentation. | `True` |
| `--refine-path / --no-refine-path` | `--no-refine-path` (default) → single-pass `path-opt`; `--refine-path` → recursive `path-search`. Both modes support Stage 5 (TSOPT / thermo / DFT). | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding correction for MM-to-ML environmental effects (experimental). | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap / --no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in the model parm7. Disabled by default, consistent with Gaussian ONIOM. | `--no-cmap` |
| `--hessian-calc-mode CHOICE` | ML/MM Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B = 0 / 10 / 20). If disabled, downstream tools require `--model-pdb` or `--model-indices`. | `True` |

TSOPT optimizer selection order: `--opt-mode-post` (if set) → `--opt-mode` (only when explicitly provided) → TSOPT default (`hess` → RS-I-RFO).

### Scan (single-input runs)

| Option | Description | Default |
| --- | --- | --- |
| `-s, --scan-lists TEXT...` | Staged scans: `(i, j, target_Å)` tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based / --no-scan-one-based` | Override scan indexing (True = 1-based, False = 0-based). | _None_ |
| `--scan-max-step-size FLOAT` | Maximum step size (Å). | _Default_ |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV / Å²). | _Default_ |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | _Default_ |
| `--scan-preopt / --no-scan-preopt` | Override scan pre-optimization toggle. | _None_ |
| `--scan-endopt / --no-scan-endopt` | Override scan end-of-stage optimization. | _None_ |

### Post-processing + freq / DFT overrides

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt / --no-tsopt` | Run TS optimization + EulerPC IRC per reactive segment. | `False` |
| `--thermo / --no-thermo` | Run vibrational analysis (`freq`) on R / TS / P. | `False` |
| `--dft / --no-dft` | Run single-point DFT on R / TS / P. | `False` |
| `--flatten / --no-flatten` | Surplus-imaginary-mode flattening in `tsopt`. | `False` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | _Default_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | _Default_ |
| `--freq-amplitude-ang FLOAT` | Mode animation amplitude (Å). | _Default_ |
| `--freq-n-frames INT` | Frames per mode animation. | _Default_ |
| `--freq-sort TEXT` | Mode sorting behavior. | _Default_ |
| `--freq-temperature FLOAT` | Thermochemistry temperature (K). | _Default_ |
| `--freq-pressure FLOAT` | Thermochemistry pressure (atm). | _Default_ |
| `--dft-out-dir PATH` | Base directory override for DFT outputs. | _None_ |
| `--dft-func-basis TEXT` | Functional / basis pair. | _Default_ |
| `--dft-max-cycle INT` | Maximum SCF iterations. | _Default_ |
| `--dft-conv-tol FLOAT` | SCF convergence tolerance. | _Default_ |
| `--dft-grid-level INT` | PySCF grid level. | _Default_ |
| `--dft-engine [gpu\|cpu]` | DFT engine (GPU or CPU PySCF). | _None_ |

## YAML configuration

`all` supports layered YAML — `--config FILE` for base settings, with the precedence `defaults < config < CLI < override-yaml`. The effective YAML is forwarded to downstream subcommands, and each tool reads the sections described in its own documentation:

| Subcommand | YAML sections |
|---|---|
| [`path-search`](path-search.md) | `geom`, `calc` / `mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc` / `mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc` / `mlmm`, `opt`, `hessian_dimer`, `rsirfo` |
| [`freq`](freq.md) | `geom`, `calc` / `mlmm`, `freq`, `thermo` |
| [`dft`](dft.md) | `dft` |

```yaml
# Minimal example
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  backend: uma                      # uma | orb | mace | aimnet2
  embedcharge: false                # xTB point-charge embedding correction
  uma_model: uma-s-1p1              # uma-s-1p1 | uma-m-1p1
  hessian_calc_mode: Analytical     # recommended when VRAM permits
gs:
  max_nodes: 12
  climb: true
dft:
  grid_level: 6
```

Full schema: [YAML Reference](yaml-reference.md).

## Notes

Input format depends on extraction:

- Extraction enabled (`-c/--center`): inputs must be **PDB** so residues can be located.
- Extraction skipped: inputs may be **PDB / XYZ**.
- Multi-structure runs require ≥ 2 structures.

Charge is resolved in order of priority — `-q/--charge` (explicit CLI override) → pocket extraction (when `-c` is provided, summing amino acids + ions + `--ligand-charge`) → `-l, --ligand-charge` fallback (when extraction is skipped) → default (unresolved charge is an error). Spin resolution: `--multiplicity` (CLI) → default (1). Always provide `--ligand-charge` for non-standard substrates so the correct net charge propagates downstream. The first-model net ML-region charge is cast to the nearest integer, with a console note if rounding occurs.

## See Also

[extract](extract.md) (called internally by `all`) · [mm_parm](mm-parm.md) (called internally by `all`) · [path-search](path-search.md) · [tsopt](tsopt.md) · [freq](freq.md) · [dft](dft.md) · [trj2fig](trj2fig.md) · [Common Error Recipes](recipes-common-errors.md) (symptom-first failure routing) · [Troubleshooting](troubleshooting.md) (common errors and fixes) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
