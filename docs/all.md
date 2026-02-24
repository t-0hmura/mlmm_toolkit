# `all`

## Overview

> **Summary:** End-to-end enzymatic reaction workflow -- pocket extraction, optional staged scan, recursive MEP search (GSM), merge to full systems, with optional TS optimization, pseudo-IRC, thermochemistry, DFT, and DFT//UMA diagrams.

`mlmm all` runs a one-shot pipeline centered on pocket models. It supports three modes:

- **Multi-structure ensemble** -- Provide two or more full PDBs in reaction order. The tool extracts pockets, builds MM topology, runs recursive GSM MEP search, merges back into the full system, and optionally runs per-segment post-processing (TSOPT/freq/DFT).
- **Single-structure + staged scan** -- Provide one PDB plus `--scan-lists`. The scan generates intermediate/product candidates that become MEP endpoints.
- **TSOPT-only** -- Provide a single PDB, omit `--scan-lists`, and set `--tsopt`. The tool runs TS optimization on the pocket, performs pseudo-IRC, minimizes both ends, and builds energy diagrams.

## Minimal example

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --out-dir ./result_all
```

## Output checklist

- `result_all/summary.log`
- `result_all/summary.yaml`
- `result_all/path_search/mep.pdb` (or `result_all/path_search/seg_*/`)

## Common examples

1. Run full post-processing in one command.

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" \
 --tsopt --thermo --dft --out-dir ./result_all
```

2. Single-structure staged scan route.

```bash
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --out-dir ./result_scan_all
```

3. Validate parsing and plan only.

```bash
mlmm all -i R.pdb P.pdb -c "SAM,GPP" --ligand-charge "SAM:1,GPP:-3" --dry-run
```

## Usage

```bash
# Standard (multi-structure ensemble in reaction order)
mlmm all -i R.pdb [I1.pdb...] P.pdb -c <substrate-spec> [--ligand-charge <map-or-number>]
 [--multiplicity <2S+1>] [--max-nodes N] [--max-cycles N]
 [--climb/--no-climb] [--opt-mode light|heavy] [--opt-mode-post light|heavy]
 [--preopt/--no-preopt] [--hessian-calc-mode Analytical|FiniteDifference] [--out-dir DIR]
 [--tsopt/--no-tsopt] [--thermo/--no-thermo] [--dft/--no-dft]
 [--tsopt-max-cycles N] [--freq-* overrides] [--dft-* overrides]

# Single-structure + staged scan
mlmm all -i A.pdb -c "308,309" --scan-lists "[(12,45,1.35)]" "[(10,55,2.20)]" \
 --multiplicity 1 --opt-mode light --preopt \
 --out-dir result_all --tsopt --thermo --dft

# Single-structure TSOPT-only mode (no path search)
mlmm all -i single.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_single
```

For help output, `mlmm all --help` shows core options and `mlmm all --help-advanced` shows the full option list.

### Examples

```bash
# Minimal end-to-end run with explicit substrate and ligand charges (multi-structure)
mlmm all -i reactant.pdb product.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1"

# Full ensemble with an intermediate, residue-ID substrate spec, and full post-processing
mlmm all -i A.pdb B.pdb C.pdb -c "308,309" --ligand-charge "-1" \
 --multiplicity 1 --max-nodes 10 --max-cycles 100 --climb \
 --opt-mode light --no-dump --config params.yaml --preopt \
 --out-dir result_all --tsopt --thermo --dft

# Single-structure + scan to build an ordered series
mlmm all -i A.pdb -c "308,309" --scan-lists "[(10,55,2.20),(23,34,1.80)]" \
 --multiplicity 1 --out-dir result_scan_all --tsopt --thermo --dft

# Single-structure TSOPT-only mode (no path_search)
mlmm all -i A.pdb -c "GPP,MMT" --ligand-charge "GPP:-3,MMT:-1" \
 --tsopt --thermo --dft --out-dir result_tsopt_only
```

## Workflow

1. **Active-site pocket extraction** (multi-structure union when multiple inputs)
 - Define the substrate (`-c/--center`, by PDB, residue IDs, or residue names).
 - Optionally provide `--ligand-charge` as a total number (distributed) or a mapping (e.g., `GPP:-3,MMT:-1`).
 - The extractor writes per-input pocket PDBs under `<out-dir>/pockets/`.
 - The extractor's **first-model total pocket charge** is used as the total charge in later steps, cast to the nearest integer with a console note if rounding occurs.
 - Additional extractor toggles: `--radius`, `--radius-het2het`, `--include-H2O/--no-include-H2O`, `--exclude-backbone/--no-exclude-backbone`, `--add-linkH/--no-add-linkH`, `--selected-resn`, `--verbose/--no-verbose`.
 - If `-c/--center` is omitted, extraction is skipped and full input structures are used directly.

### Charge resolution order
1. `-q/--charge` (explicit CLI override) — highest priority.
2. Pocket extraction charge summary (sum of amino acids, ions, and `--ligand-charge`).
3. `--ligand-charge` fallback when extraction is skipped.
4. Default: none (abort if unresolved).

2. **ML/MM preparation**
 - Use the first pocket as `<out-dir>/ml_region.pdb` for `--model-pdb`. Keep `--no-add-linkH` (the default) during extraction if you prefer to omit link hydrogens from this definition.
 - Run `mm_parm` once on the first full input PDB to build `<out-dir>/mm_parm/<input_basename>.parm7` / `.rst7`, automatically passed as `--real-parm7`.
 - Tune this stage with `--auto-mm-ff-set`, `--auto-mm-add-ter`, and `--auto-mm-keep-temp`.

3. **Optional staged scan (single-structure only)**
 - If exactly one full input PDB is provided and `--scan-lists` is given, the tool performs a staged, bond-length-driven scan on the extracted pocket PDB using the UMA calculator.
 - For each stage, the final relaxed structure (`stage_XX/result.pdb`) is collected as an intermediate/product candidate.
 - The ordered input series for the path search becomes: `[initial pocket, stage_01/result.pdb, stage_02/result.pdb,...]`.

4. **MEP search (recursive GSM) on pocket inputs**
 - Runs `path_search` with options forwarded from this command.
 - For multi-input runs, the original full PDBs are supplied as merge references automatically. In the scan-derived series (single-structure case), the single original full PDB is reused (repeated) as the reference template for all pocket inputs.

5. **Merge to full systems and optional post-processing**
 - The pocket MEP is merged back into the original full-system template(s) within `<out-dir>/path_search/`.
 - Pocket-only and full-system trajectories, per-segment merged PDBs, and a summary are written.
 - `--thermo`: Compute UMA thermochemistry on (R, TS, P) and add a Gibbs diagram.
 - `--dft`: Do DFT single-point on (R, TS, P) and add a DFT diagram. With `--thermo`, also generate a DFT//UMA Gibbs diagram.

6. **TSOPT-only mode** (single input, `--tsopt`, no `--scan-lists`)
 - Skips steps (4)-(5) and runs `tsopt` on the pocket, does a pseudo-IRC and minimizes both ends, builds UMA energy diagrams for R-TS-P, and optionally adds UMA Gibbs, DFT, and DFT//UMA diagrams.
 - In this mode only, the IRC endpoint with **higher energy** is adopted as the reactant (R).

## CLI options

### Input/Output Options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full PDBs in reaction order (single input allowed with `--scan-lists` or `--tsopt`). | Required |
| `-c, --center TEXT` | Substrate specification (PDB path, residue IDs, or residue names). Omit to skip extraction. | _None_ |
| `--ligand-charge TEXT` | Total charge or residue-specific mapping (e.g., `GPP:-3,MMT:-1`). | _None_ |
| `-q, --charge INT` | Force total system charge (highest priority override). | _None_ |
| `--out-dir PATH` | Top-level output directory. | `./result_all/` |
| `--dump/--no-dump` | Save optimizer dumps. | `False` |
| `--config FILE` | Base YAML applied first. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration before execution. | `False` |
| `--dry-run/--no-dry-run` | Validate and print plan without running stages. | `False` |

### Extraction Options

| Option | Description | Default |
| --- | --- | --- |
| `--radius FLOAT` | Pocket inclusion cutoff (Angstrom). | Extractor default |
| `--radius-het2het FLOAT` | Independent hetero-hetero cutoff. | Extractor default |
| `--include-H2O/--no-include-H2O` | Include water molecules. | Extractor default |
| `--exclude-backbone/--no-exclude-backbone` | Remove backbone atoms on non-substrate amino acids. | Extractor default |
| `--add-linkH/--no-add-linkH` | Add link hydrogens for severed bonds. | Extractor default |
| `--selected-resn TEXT` | Residues to force include. | `""` |
| `--verbose/--no-verbose` | Enable INFO-level extractor logging. | Extractor default |

### MM Preparation Options

| Option | Description | Default |
| --- | --- | --- |
| `--auto-mm-ff-set TEXT` | Force field set for `mm_parm`. | `mm_parm` default |
| `--auto-mm-add-ter {True\|False}` | Add TER records. | `mm_parm` default |
| `--auto-mm-keep-temp {True\|False}` | Keep temporary files. | `mm_parm` default |

### MEP Search Options

| Option | Description | Default |
| --- | --- | --- |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--max-nodes INT` | Internal nodes for segment GSM. | `10` |
| `--max-cycles INT` | Max GSM macro-cycles. | `300` |
| `--climb/--no-climb` | Enable TS refinement for segment GSM. | `True` |
| `--opt-mode [light\|heavy]` | Optimizer preset for scan/path-search and single optimizations (`light` → LBFGS/Dimer, `heavy` → RFO/RSIRFO). | `light` |
| `--opt-mode-post [light\|heavy]` | Optimizer preset override for TSOPT/post-IRC endpoint optimizations. | _None_ |
| `--preopt/--no-preopt` | Pre-optimize endpoints before segmentation. | `True` |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode (`Analytical` or `FiniteDifference`). | _Default_ |

TSOPT optimizer selection order: `--opt-mode-post` (if set) -> `--opt-mode` (only when explicitly provided) -> TSOPT default (`heavy`).

### Scan Options (Single-Input Runs)

| Option | Description | Default |
| --- | --- | --- |
| `--scan-lists TEXT...` | Staged scans: `(i,j,target_A)` tuples. | _None_ |
| `--scan-out-dir PATH` | Override the scan output directory. | _None_ |
| `--scan-one-based/--no-scan-one-based` | Force 1-based indexing for atom selectors. | `True` |
| `--scan-max-step-size FLOAT` | Maximum step size (Angstrom). | _Default_ |
| `--scan-bias-k FLOAT` | Harmonic bias strength (eV/Angstrom^2). | _Default_ |
| `--scan-relax-max-cycles INT` | Relaxation max cycles per step. | _Default_ |
| `--scan-preopt/--no-scan-preopt` | Override scan pre-optimization toggle. | _Default_ |
| `--scan-endopt/--no-scan-endopt` | Override scan end-of-stage optimization. | _Default_ |

### Post-Processing Options

| Option | Description | Default |
| --- | --- | --- |
| `--tsopt/--no-tsopt` | Run TS optimization + pseudo-IRC per reactive segment. | `False` |
| `--thermo/--no-thermo` | Run vibrational analysis (`freq`) on R/TS/P. | `False` |
| `--dft/--no-dft` | Run single-point DFT on R/TS/P. | `False` |
| `--tsopt-max-cycles INT` | Override `tsopt --max-cycles`. | _Default_ |
| `--tsopt-out-dir PATH` | Custom tsopt subdirectory. | _None_ |

### Freq Overrides

| Option | Description | Default |
| --- | --- | --- |
| `--freq-out-dir PATH` | Base directory override for freq outputs. | _None_ |
| `--freq-max-write INT` | Maximum modes to write. | _Default_ |
| `--freq-amplitude-ang FLOAT` | Mode animation amplitude (Angstrom). | _Default_ |
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
 ml_region.pdb # ML-region definition (copy of the first pocket)
 pockets/
 pocket_<input1_basename>.pdb
 pocket_<input2_basename>.pdb
...
 mm_parm/
 <input1_basename>.parm7 # Generated from the first full-enzyme input PDB
 <input1_basename>.rst7
 scan/ # present only in single-structure+scan mode
 stage_01/result.pdb
 stage_02/result.pdb
...
 summary.yaml # mirrored top-level summary (when path_search runs)
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 energy_diagram_UMA_all.png # aggregated post-processing diagrams (when enabled)
 energy_diagram_G_UMA_all.png
 energy_diagram_DFT_all.png
 energy_diagram_G_DFT_plus_UMA_all.png
 irc_plot_all.png
 path_search/ # present when path_search is executed
 mep_trj.xyz
 mep.pdb
 mep_w_ref.pdb
 mep_w_ref_seg_XX.pdb
 summary.yaml
 summary.log
 mep_plot.png
 energy_diagram_MEP.png
 post_seg_XX/ # when post-processing is enabled
 ts/...
 irc/...
 freq/... # with --thermo
 dft/... # with --dft
 energy_diagram_UMA.png
 energy_diagram_G_UMA.png
 energy_diagram_DFT.png
 energy_diagram_G_DFT_plus_UMA.png
 tsopt_single/ # present only in single-structure TSOPT-only mode
 ts/...
 irc/...
 structures/
 reactant.pdb
 ts.pdb
 product.pdb
 freq/... # with --thermo
 dft/... # with --dft
 energy_diagram_UMA.png
 energy_diagram_G_UMA.png
 energy_diagram_DFT.png
 energy_diagram_G_DFT_plus_UMA.png
```

## YAML configuration

`all` now supports two YAML layers:

- `--config FILE`: base settings.

`defaults < config < CLI < override-yaml`

The resulting effective YAML is forwarded to downstream subcommands. Each tool reads the sections described in its own documentation:

| Subcommand | YAML Sections |
|------------|---------------|
| [`path-search`](path_search.md) | `geom`, `calc`/`mlmm`, `gs`, `opt`, `lbfgs`, `bond`, `search` |
| [`scan`](scan.md) | `geom`, `calc`/`mlmm`, `opt`, `lbfgs` |
| [`tsopt`](tsopt.md) | `geom`, `calc`/`mlmm`, `opt` |
| [`freq`](freq.md) | `geom`, `calc`/`mlmm`, `freq` |
| [`dft`](dft.md) | `dft` |

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- **Python >= 3.11** is required.
- **Substrate (`-c/--center`) and (when needed) `--ligand-charge` are practically required.**
- Single-structure mode requires either `--scan-lists` or `--tsopt`; otherwise at least two structures are needed.
- Preflight checks now validate repeated `-i/--input` values as files and verify AmberTools commands (`tleap`, `antechamber`, `parmchk2`) before auto `mm_parm`.
- Energies in diagrams are plotted relative to the first state in kcal/mol (converted from Hartree).
- Charge handling: the extractor's first-model total pocket charge is used as the path/scan/TSOPT total charge (rounded to int).

---

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
