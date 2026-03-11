# `path-search`

## Overview

> **Summary:** Build a continuous MEP from two or more structures with recursive GSM segmentation. Automatically refines only regions with bond changes and exports the highest-energy image (HEI) as a TS candidate.

### At a glance
- **Use when:** You have R -> ... -> P structures (2+ inputs) and want a single stitched MEP with automatic refinement.
- **Method:** Chains GSM segments and recursively refines only sub-intervals that still contain covalent changes.
- **Outputs:** `mep_trj.xyz` (main trajectory), `summary.yaml` (segment-by-segment results), and optional plots/merged PDBs when enabled.
- **Defaults:** `--opt-mode grad` (LBFGS), `--preopt`, `--align`, `--thresh gau_loose` (GSM) / `gau` (single-structure).
- **Next step:** HEI output alone does **not** validate a TS. Follow with [tsopt](tsopt.md), [freq](freq.md), and [irc](irc.md).

`mlmm path-search` builds a continuous minimum-energy path (MEP) across two or more structures using GSM. It selectively refines only those regions where covalent bond changes are detected, then stitches the resolved subpaths into a single trajectory.

If you only have **two** endpoints and do not need recursive refinement, [path-opt](path_opt.md) is the simpler option.

## Minimal example

```bash
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --out-dir ./result_path_search
```

## Output checklist

- `result_path_search/mep_trj.xyz`
- `result_path_search/summary.yaml`
- `result_path_search/summary.log`
- `result_path_search/mep_plot.png` (when plotting succeeds)

## Common examples

1. Build a multistep path with explicit intermediates.

```bash
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q -1 --out-dir ./result_path_search_multi
```

2. Merge pocket trajectories back into a full template.

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --ref-pdb holo_template.pdb \
 --out-dir ./result_path_search_merge
```

3. Run a lighter pass without pre-optimization or alignment.

```bash
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --no-preopt --no-align --max-nodes 8 \
 --out-dir ./result_path_search_fast
```

## Usage

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb \
 --parm real.parm7 --model-pdb ml_region.pdb -q CHARGE [-m MULT]
 [--mep-mode gsm|dmf] [--refine-mode peak|minima]
 [--freeze-atoms "1,3,5"] [--max-nodes N] [--max-cycles N] [--climb/--no-climb]
 [--opt-mode grad|hess]
 [--thresh PRESET] [--dump/--no-dump] [--out-dir DIR]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

### Examples

```bash
# Minimal pocket-only MEP between two states
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0

# Multistep path with YAML overrides, frozen atoms, and merged full-system output
mlmm path-search -i R.pdb IM1.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q -1 --freeze-atoms "1,3,5" \
 --ref-pdb holo_template.pdb --out-dir ./run_ps
```

## Workflow

1. **Initial segment per pair (GSM/DMF)** -- Run the selected MEP engine (`--mep-mode`) between each adjacent input (A->B) to obtain a coarse MEP and identify the highest-energy image (HEI).
2. **Local relaxation around HEI** -- Seed refinement from `--refine-mode` (`peak`: HEI+/-1, `minima`: nearest local minima), then optimize with the chosen single-structure optimizer (`opt-mode`) to recover nearby minima (`End1`, `End2`).
3. **Decide between kink vs. refinement**:
 - If no covalent bond change is detected between `End1` and `End2`, treat the region as a *kink*: insert `search.kink_max_nodes` linear nodes and optimize each individually.
 - Otherwise, launch a **refinement segment (GSM)** between `End1` and `End2` to sharpen the barrier.
4. **Selective recursion** -- Compare bond changes for `(A->End1)` and `(End2->B)` using the `bond` thresholds. Recurse only on sub-intervals that still contain covalent updates. Recursion depth is capped by `search.max_depth`.
5. **Stitching & bridging** -- Concatenate resolved subpaths, dropping duplicate endpoints when RMSD <= `search.stitch_rmsd_thresh`. If the RMSD gap between two stitched pieces exceeds `search.bridge_rmsd_thresh`, insert a bridge MEP segment using the selected `--mep-mode`. When the interface itself shows a bond change, a brand-new recursive segment replaces the bridge.
6. **Optional alignment and merge** -- After pre-opt, `--align` rigidly co-aligns inputs and refines freezes. With `--ref-pdb`, pocket trajectories merge into full templates and segments are annotated for plotting/analysis.

Bond-change detection relies on `bond_changes.compare_structures` with thresholds surfaced under the `bond` YAML section.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more full-enzyme PDBs in reaction order. Repeat `-i` or pass multiple paths after one flag. | Required |
| `--parm PATH` | Amber parm7 topology for the full enzyme complex. | Required |
| `--model-pdb PATH` | PDB defining the ML (high-level) region atoms for ML/MM. Optional when `--detect-layer` or `--model-indices` is used. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed like `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20). If disabled, you must provide `--model-pdb` or `--model-indices`. | `True` |
| `-q, --charge INT` | Charge of the ML region (integer). | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--mep-mode [gsm\|dmf]` | MEP backend for segment/bridge searches. | `gsm` |
| `--refine-mode [peak\|minima]` | HEI refinement seed rule. | `peak` for `gsm`, `minima` for `dmf` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. Applied to movable MM atoms. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. Providing `--movable-cutoff` disables `--detect-layer`. | _None_ |
| `--max-nodes INT` | Internal nodes for segment GSM. | `10` |
| `--max-cycles INT` | Max GSM macro-cycles. | `300` |
| `--climb/--no-climb` | Enable TS refinement for segment GSM. | `True` |
| `--opt-mode [grad\|hess]` | Single-structure optimizer preset (`grad` = LBFGS, `hess` = RFO). | `grad` |
| `--preopt/--no-preopt` | Pre-optimize endpoints with LBFGS before segmentation. | `True` |
| `--align / --no-align` | Rigidly align inputs after pre-opt. | `True` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (effective: `gau_loose`) |
| `--dump/--no-dump` | Save optimizer dumps. | `False` |
| `-o, --out-dir PATH` | Output directory. | `./result_path_search/` |
| `--ref-pdb PATH...` | Full template PDB(s) for final merge. | _None_ |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layer metadata) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running path search. Shown in `--help-advanced`. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## Outputs

```text
out_dir/ (default: ./result_path_search/)
 summary.yaml # MEP-level run summary (no full settings dump)
 summary.log # Human-readable summary
 mep_trj.xyz # Final MEP (always written)
 mep.pdb # Final MEP (PDB when ref template available)
 mep_w_ref.pdb # Full-system merged MEP (requires --ref-pdb)
 mep_w_ref_seg_XX.pdb # Per-segment merged MEPs (bond-change segments; requires --ref-pdb)
 mep_seg_XX_trj.xyz / mep_seg_XX.pdb # Pocket-only per-segment paths
 hei_seg_XX.xyz / hei_seg_XX.pdb # Pocket HEI and optional PDB per bond-change segment
 hei_w_ref_seg_XX.pdb # Merged HEI per bond-change segment (requires --ref-pdb)
 mep_plot.png # Delta-E profile vs image index (from trj2fig)
 energy_diagram_MEP.png # State-level energy diagram relative to the reactant (kcal/mol)
 seg_000_*/ # Segment-level GSM and refinement artifacts
```

## YAML configuration

Merge order is **defaults < config < explicit CLI < override**.
The YAML root must be a mapping. Accepted sections:

- **`geom`** -- `coord_type` (`"cart"` default), `freeze_atoms` (0-based indices).
- **`calc` / `mlmm`** -- ML/MM calculator settings: `input_pdb`, `real_parm7`, `model_pdb`, `model_charge`, `model_mult`, backend selection (`backend`, `embedcharge`), UMA controls (`uma_model`, `uma_task_name`, `ml_hessian_mode`), device selection, freeze atoms.
- **`gs`** -- Growing String settings: `max_nodes`, `climb`, `climb_rms`, `climb_fixed`, `reparam_every_full`, `reparam_check`.
- **`opt`** -- StringOptimizer controls: `max_cycles`, `print_every`, `dump`, `dump_restart`, `out_dir`.
- **`lbfgs`** -- Single-structure optimizer controls for HEI+/-1 refinement: `keep_last`, `beta`, `gamma_mult`, `max_step`, `control_step`, `double_damp`, `mu_reg`, `max_mu_reg_adaptions`.
- **`bond`** -- Bond-change detection: `bond_factor`, `margin_fraction`, `delta_fraction`.
- **`search`** -- Recursion logic: `max_depth`, `stitch_rmsd_thresh`, `bridge_rmsd_thresh`, `max_nodes_segment`, `max_nodes_bridge`, `kink_max_nodes`, `max_seq_kink`, `refine_mode`.

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Inputs: provide at least two full-enzyme PDBs to `-i/--input` in reaction order.
- Preflight checks validate `-i/--input` and `--ref-pdb` paths before starting GSM.
- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli_conventions.md).
- Freeze atoms: `--freeze-atoms "1,3,5"` stores zero-based indices and merges with YAML `geom.freeze_atoms`.
- Nodes and recursion: segment vs bridge nodes differ via `search.max_nodes_segment` and `search.max_nodes_bridge`. Kinks use `search.kink_max_nodes` (default 3) linear nodes. Recursion depth is capped by `search.max_depth` (default 10).
- Optimizers: `--mep-mode gsm` uses pysisyphus `GrowingString` + `StringOptimizer`; `--mep-mode dmf` uses Direct Max Flux. Single-structure refinements always use LBFGS.
- Final merge rule with `--align`: when `--ref-pdb` is provided, the first reference PDB is used for all pairs.
- Console output prints the state sequence (e.g., `R --> TS1 --> IM1 -->... --> P`) plus the labels/energies used to build the energy diagram.
- `summary.log` rendering is resilient to missing payload fields. Internally, defaults are applied for:
 `root_out_dir`, `path_module_dir`, `pipeline_mode`, `segments`, `energy_diagrams`.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path-opt](path_opt.md) -- Single-pass MEP optimization (no recursive refinement)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow that calls path-search internally
- [trj2fig](trj2fig.md) -- Plot energy profiles from MEP trajectories
