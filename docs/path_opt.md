# `path-opt`

## Overview

> **Summary:** Find an MEP between **exactly two** enzyme structures using GSM or DMF with the ML/MM calculator. Writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate.

`mlmm path-opt` optimizes a minimum-energy path between two enzyme states using PySisyphus `GrowingString` with the ML/MM calculator. The ML/MM calculator keeps the full enzyme complex without link atoms: the ML region is defined by `--model-pdb`, the Amber topology comes from `--parm`, and both endpoints are supplied as PDBs containing the full system coordinates.

For workflows that start from **two or more** structures and automatically refine only the reactive region, use [path-search](path_search.md).

## Minimal example

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_path_opt
```

## Output checklist

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb` (when PDB conversion is available)

## Common examples

1. Pre-optimize both endpoints before path growth.

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. Disable climbing-image refinement for a quick first pass.

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-climb --max-nodes 8 --out-dir ./result_path_opt_fast
```

3. Freeze selected atoms and keep optimizer dumps.

```bash
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_path_opt_dump
```

## Usage
```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --parm real.parm7 --model-pdb model.pdb \
 -q CHARGE [-m MULT] [--mep-mode gsm|dmf] [--fix-ends/--no-fix-ends] [options]
```

### Examples
```bash
# Minimal invocation
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0

# With frozen atoms, more nodes, and layered YAML overrides
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
 --freeze-atoms "1,3,5,7" --max-nodes 10 --max-cycles 200 --dump --out-dir ./result_path_opt/

# With endpoint pre-optimization
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
 --preopt --preopt-max-cycles 20000
```

## Workflow
1. **Load endpoints** -- Read both PDB structures and resolve charge/spin from CLI or defaults.
 Set up the ML/MM calculator with `--parm`, `--model-pdb`, and charge/spin.
2. **Pre-alignment** -- All endpoints after the first are Kabsch-aligned to the first
 structure. If `freeze_atoms` is defined, only those atoms participate in the RMSD
 fit; the resulting transform is applied to all atoms.
3. **Optional pre-optimization** -- With `--preopt`, each endpoint is pre-optimized
 by LBFGS (using the same ML/MM calculator) before alignment and string growth.
 The number of LBFGS cycles is controlled by `--preopt-max-cycles` (default: 10000).
4. **Path optimization** -- `--mep-mode gsm` uses PySisyphus `GrowingString` with `(max_nodes + 2)` images including endpoints; `--mep-mode dmf` uses Direct Max Flux.
5. **Climbing image (GSM only)** -- With `--climb`, a climbing-image refinement is applied after string growth, and the highest-energy image (HEI) is reported.
6. **Output** -- Final path trajectory and HEI are written as XYZ and PDB files.
 PDB conversion is performed when the inputs are PDBs.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product PDB structures. | Required |
| `--parm PATH` | Amber prmtop for the full REAL system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs). Optional when `--detect-layer` or `--model-indices` is used. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed like `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20). If disabled, you must provide `--model-pdb` or `--model-indices`. | `True` |
| `-q, --charge INT` | Total ML-region charge. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-residue charge map, e.g. `SAM:1,PHN:-1`. Derives total charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--mep-mode [gsm\|dmf]` | MEP backend. | `gsm` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze (converted to 0-based; merged with YAML `geom.freeze_atoms`). | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (Å) from ML region for MM atoms to include in Hessian calculation. Applied to movable MM atoms. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. Providing `--movable-cutoff` disables `--detect-layer`. | _None_ |
| `--fix-ends/--no-fix-ends` | Fix endpoint structures during GSM growth (`gs.fix_first/fix_last`). The effective default is `True` because `gs.fix_first` and `gs.fix_last` default to `True` in `defaults.py`. | `False` (CLI); `True` (effective via GS_KW) |
| `--max-nodes INT` | Number of internal string nodes (total images = `max_nodes + 2`). | `20` |
| `--max-cycles INT` | Optimizer macro-iteration cap (growth + refinement). Also sets `opt.stop_in_when_full`. | `300` |
| `--climb/--no-climb` | Enable climbing-image refinement after full string growth. | `True` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with LBFGS before alignment/string growth. | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint pre-optimization cycles. | `10000` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (effective: `gau_loose`) |
| `--dump/--no-dump` | Dump optimizer trajectories and restarts inside `out_dir`. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layers) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running optimization. Shown in `--help-advanced`. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## Outputs
```
out_dir/ (default: ./result_path_opt/)
├─ final_geometries_trj.xyz # XYZ trajectory with per-image energies in the comment line
├─ final_geometries.pdb # Same path as final_geometries_trj.xyz, mapped back to the reference PDB ordering
├─ hei.xyz # Highest-energy image (XYZ, always written)
├─ hei.pdb # HEI in PDB format (when reference PDB is available)
├─ align_refine/ # External alignment/refinement artifacts
├─ preopt/ # Endpoint pre-optimization outputs (present when --preopt)
└─ <optimizer dumps> # Present when --dump or opt.dump_restart > 0
```

## YAML configuration

Merge order is **defaults < config < explicit CLI < override**.

### Section `geom`
- `coord_type`: Coordinate type (cartesian vs dlc internals).
- `freeze_atoms`: 0-based frozen atoms merged with CLI `--freeze-atoms`.

### Section `calc` / `mlmm`
- ML/MM calculator setup: `charge`, `spin`, `backend`, `embedcharge`, UMA-specific `model`/`task_name`, `device`, neighbor radii, Hessian options, etc.

### Section `gs`
- Growing String controls: `max_nodes`, `perp_thresh`, reparameterization cadence, `max_micro_cycles`, DLC resets, climb toggles/thresholds.

### Section `opt`
- StringOptimizer settings: `stop_in_when_full`, `scale_step`, `max_cycles`, dumping flags, `reparam_thresh`, `coord_diff_thresh`, `out_dir`, `print_every`.

## Exit codes

| Code | Meaning |
| --- | --- |
| `0` | Success |
| `3` | Optimization failure |
| `4` | Final trajectory write error |
| `5` | HEI dump error |
| `130` | Keyboard interrupt |
| `1` | Unhandled exception |

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path-search](path_search.md) -- Recursive MEP search with automatic refinement (for 2+ structures)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow (uses path-search by default)
- [YAML Reference](yaml_reference.md) -- Full `gs`, `opt` configuration options
