# `path-search`

`mlmm path-search` builds a continuous minimum-energy path (MEP) across two or more structures using the selected MEP engine (GSM by default, or DMF). It selectively refines only those regions where covalent bond changes are detected, then stitches the resolved subpaths into a single trajectory. Use it to drive a multistep mechanism from R + (optional intermediates) + P, where the recursive segmentation auto-detects elementary steps. Complex multistep mechanisms may require manual trial-and-error—adjusting input intermediates, MEP-engine settings, or convergence thresholds—to obtain a satisfactory pathway.

## Examples

```bash
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --out-dir ./result_path_search
```

Build a multistep path with explicit intermediates:

```bash
# Build a multistep path with explicit intermediates
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q -1 --out-dir ./result_path_search_multi
```

Lighter pass without pre-optimization or alignment:

```bash
# Lighter pass without pre-optimization or alignment
mlmm path-search -i reactant.pdb product.pdb --parm real.parm7 \
 --model-pdb ml_region.pdb -q 0 --no-preopt --no-align --max-nodes 8 \
 --out-dir ./result_path_search_fast
```

General command form:

```bash
mlmm path-search -i R.pdb IM1.pdb P.pdb \
 --parm real.parm7 --model-pdb ml_region.pdb -q CHARGE [-m MULT]
 [--mep-mode gsm|dmf] [--refine-mode peak|minima]
 [--freeze-atoms "1,3,5"] [--max-nodes N] [--max-cycles-gsm N] [--max-cycles-dmf N] [--climb/--no-climb]
 [--thresh PRESET] [--dump/--no-dump] [--out-dir DIR]
 [--show-config/--no-show-config] [--dry-run/--no-dry-run]
```

## Workflow

1. **Initial segment per pair (GSM/DMF)** -- Run the selected MEP engine (`--mep-mode`) between each adjacent input (A->B) to obtain a coarse MEP and identify the highest-energy image (HEI).
2. **Local relaxation around HEI** -- Seed refinement from `--refine-mode` (`peak`: HEI+/-1, `minima`: nearest local minima), then use L-BFGS to recover nearby minima (`End1`, `End2`).
3. **Decide between kink vs. refinement**:
 - If no covalent bond change is detected between `End1` and `End2`, treat the region as a *kink*: insert `search.kink_max_nodes` linear nodes and optimize each individually.
 - Otherwise, launch a **refinement segment with the selected MEP engine** between `End1` and `End2` to sharpen the barrier.
4. **Selective recursion** -- Compare bond changes for `(A->End1)` and `(End2->B)` using the `bond` thresholds. Recurse only on sub-intervals that still contain covalent bond changes. Recursion depth is capped by `search.max_depth`.
5. **Stitching & bridging** -- Concatenate resolved subpaths, dropping duplicate endpoints when RMSD <= `search.stitch_rmsd_thresh`. If the RMSD gap between two stitched pieces exceeds `search.bridge_rmsd_thresh`, insert a bridge MEP segment using the selected `--mep-mode`. When the interface itself shows a bond change, a new recursive segment replaces the bridge.
6. **Optional alignment/refinement** -- After optional preoptimization, `--align` rigidly aligns inputs to the first input. With frozen anchors, the shared owner also performs a freeze-guided scan and L-BFGS relaxation toward the reference, then re-matches the freeze-atom selection. Segments are annotated for plotting/analysis.

Bond-change detection relies on `bond_changes.compare_structures` with thresholds surfaced under the `bond` YAML section.

## Outputs

```text
out_dir/ (default: ./result_path_search/)
 summary.json # MEP-level run summary (no full settings dump)
 summary.log # Human-readable summary
 mep_trj.xyz # Final MEP (always written)
 mep.pdb # Final MEP (PDB when ref template available)
 mep_seg_XX_trj.xyz / mep_seg_XX.pdb # Per-segment paths
 hei_seg_XX.xyz / hei_seg_XX.pdb # HEI per bond-change segment
 mep_plot.png # Delta-E profile vs image index (from trj2fig)
 energy_diagram_MEP.png # State-level energy diagram relative to the reactant (kcal/mol)
 seg_000_*/ # Segment-level GSM and refinement artifacts
```

## CLI options

`mlmm path-search --help` shows core options; `mlmm path-search --help-advanced` shows the full option list. The full flag list is also in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH...` | Two or more PDB/mmCIF structures, or XYZ files with corresponding `--ref-pdb` entries, in reaction order. Repeat `-i` or pass multiple paths after one flag. | Required |
| `--parm PATH` | Amber parm7 topology for the full enzyme complex. | Required |
| `--model-pdb PATH` | PDB defining the ML (high-level) region atoms for ML/MM. Optional when `--detect-layer` or `--model-indices` is used. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed like `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer` | Automatically read B-factor layers (B=0/10/20). With explicit ML membership, only the MM sublayers are retained; otherwise B-factors also define ML membership. | Enabled |
| `-q, --charge INT` | Net charge of the ML region (integer). Required unless `--ligand-charge` is provided. | _None_ |
| `-l, --ligand-charge TEXT` | Per-residue charge map, e.g. `SAM:1,PHN:-1`. Derives total charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--mep-mode [gsm\|dmf]` | MEP backend for segment/bridge searches. | `gsm` |
| `--dmf-backend [cpu\|gpu]` | DMF compute backend (`--mep-mode dmf` only): `gpu` (`dmf.torch`/CUDA) or `cpu` (`dmf`/NumPy). Retry `cpu` on a GPU out-of-memory error. Requires `pydmf>=1.2`. | `gpu` |
| `--refine-mode [peak\|minima]` | HEI refinement seed rule. | `peak` for `gsm`, `minima` for `dmf` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. Providing `--movable-cutoff` disables `--detect-layer`. | _None_ |
| `--max-nodes INT` | Movable internal images per GSM or DMF segment (`max_nodes + 2` total images). | `20` |
| `--max-cycles-gsm INT` | Cycle budget for the GSM string optimizer. | `300` |
| `--max-cycles-dmf INT` | IPOPT iteration budget for DMF. | `300` |
| `--climb/--no-climb` | Enable TS refinement for segment GSM. | `True` |
| `--preopt/--no-preopt` | Pre-optimize endpoints with L-BFGS before segmentation. | `True` |
| `--align/--no-align` | After preoptimization, align inputs and, with frozen anchors, run freeze-guided scan/relaxation before re-matching freeze atoms. | `True` |
| `--thresh TEXT` | Convergence preset for single-structure L-BFGS runs only (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (effective: `gau`) |
| `--thresh-gsm TEXT` | Convergence preset for the GSM string optimizer (`stopt.thresh`; same presets as `--thresh`). | _None_ (effective: `gau_loose`) |
| `--thresh-dmf TEXT` | IPOPT dual-infeasibility tolerance of the DMF optimizer (`dmf.tol`): `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive float. Gaussian presets are rejected. | _None_ (effective: `tight`) |
| `--mm-backend [hessian_ff\|openmm]` | MM backend. Hessians use finite differences by default; set `calc.mm_fd: false` for the `hessian_ff` analytical path. | `hessian_ff` |
| `--dump/--no-dump` | Save optimizer dumps. | `False` |
| `-o, --out-dir PATH` | Output directory. | `./result_path_search/` |
| `--ref-pdb PATH...` | Full template PDB(s) for XYZ→PDB conversion and topology reference. | _None_ |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layer metadata) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running path search. Shown in `--help-advanced`. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--cmap/--no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## YAML configuration

Merge order is **defaults < config < explicit CLI**. The YAML root must be a mapping. The relevant sections are `geom`/`calc`(alias `mlmm`)/`gs`/`opt` (shared with `path-opt`) plus `lbfgs` (HEI+/-1 single-structure refinement), `bond` (bond-change detection), and `search` (recursive segmentation logic, path-search only).

```yaml
# Minimal path-search YAML (every key and default: see YAML Reference)
calc:
  backend: uma
search:
  max_depth: 10            # recursion depth cap
  refine_mode: null        # peak | minima | null (auto)
bond:
  bond_factor: 1.2         # covalent-radius scaling for bond-change cutoff
```

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Notes

- If you only have **two** endpoints and do not need recursive refinement, prefer [path-opt](path-opt.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [path-opt](path-opt.md) — Single-pass MEP optimization (no recursive refinement)
- [opt](opt.md) — Single-structure geometry optimization
- [all](all.md) — End-to-end workflow (uses single-pass path-opt by default; `--refine-path` for recursive path-search)
- [trj2fig](trj2fig.md) — Plot energy profiles from MEP trajectories
- [YAML Reference](yaml-reference.md) — Full `gs`, `bond`, `search` configuration options
