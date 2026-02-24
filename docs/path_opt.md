# `path-opt`

## Overview

> **Summary:** Find an MEP between **exactly two** enzyme structures using the Growing String Method (GSM) with the ML/MM calculator. Writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate.

### Quick reference
- **Use when:** You have reactant and product endpoints (R -> P) and want a first-pass MEP.
- **Method:** PySisyphus `GrowingString` with the ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`).
- **Outputs:** `final_geometries_trj.xyz` (path) and `hei.xyz` (HEI), plus optional `.pdb` companions.
- **Defaults:** `--climb`, `--max-nodes 10`, `--max-cycles 300`.
- **Next step:** Validate the HEI with `tsopt` -> `freq` (expect one imaginary mode) -> `irc`.

`mlmm path-opt` optimizes a minimum-energy path between two enzyme states using PySisyphus `GrowingString` with the ML/MM calculator. The ML/MM calculator keeps the full enzyme complex without link atoms: the ML region is defined by `--model-pdb`, the Amber topology comes from `--real-parm7`, and both endpoints are supplied as PDBs containing the full system coordinates.


## Minimal example

```bash
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_path_opt
```

## Output checklist

- `result_path_opt/final_geometries_trj.xyz`
- `result_path_opt/hei.xyz`
- `result_path_opt/hei.pdb` (when PDB conversion is available)

## Common examples

1. Pre-optimize both endpoints before path growth.

```bash
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

2. Disable climbing-image refinement for a quick first pass.

```bash
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-climb --max-nodes 8 --out-dir ./result_path_opt_fast
```

3. Freeze selected atoms and keep optimizer dumps.

```bash
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_path_opt_dump
```

## Usage
```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
 -q CHARGE [-m MULT] [options]
```

### Examples
```bash
# Minimal invocation
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# With frozen atoms, more nodes, and layered YAML overrides
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
 --freeze-atoms "1,3,5,7" --max-nodes 10 --max-cycles 200 --dump --out-dir ./result_path_opt/ \

# With endpoint pre-optimization
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
 --preopt --preopt-max-cycles 20000
```

## Workflow
1. **Load endpoints** -- Read both PDB structures and resolve charge/spin from CLI or defaults.
 Set up the ML/MM calculator with `--real-parm7`, `--model-pdb`, and charge/spin.
2. **Pre-alignment** -- All endpoints after the first are Kabsch-aligned to the first
 structure. If `freeze_atoms` is defined, only those atoms participate in the RMSD
 fit; the resulting transform is applied to all atoms.
3. **Optional pre-optimization** -- With `--preopt`, each endpoint is pre-optimized
 by LBFGS (using the same ML/MM calculator) before alignment and string growth.
 The number of LBFGS cycles is controlled by `--preopt-max-cycles` (default: 10000).
4. **String growth** -- PySisyphus `GrowingString` grows the path between the endpoints
 using `(max_nodes + 2)` images including endpoints.
5. **Climbing image** -- With `--climb`, a climbing-image refinement is applied after
 the string grows fully, and the highest-energy image (HEI) is reported.
6. **Output** -- Final path trajectory and HEI are written as XYZ and PDB files.
 PDB conversion is performed when the inputs are PDBs.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product PDB structures. | Required |
| `--real-parm7 PATH` | Amber prmtop for the full REAL system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs). Optional when `--detect-layer` or `--model-indices` is used. | _None_ |
| `-q, --charge INT` | Total ML-region charge. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze (converted to 0-based; merged with YAML `geom.freeze_atoms`). | _None_ |
| `--max-nodes INT` | Number of internal string nodes (total images = `max_nodes + 2`). | `10` |
| `--max-cycles INT` | Optimizer macro-iteration cap (growth + refinement). Also sets `opt.stop_in_when_full`. | `300` |
| `--climb/--no-climb` | Enable climbing-image refinement after full string growth. | `True` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with LBFGS before alignment/string growth. | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint pre-optimization cycles. | `10000` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--dump/--no-dump` | Dump optimizer trajectories and restarts inside `out_dir`. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layers) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running optimization. | `False` |

## Outputs
```
out_dir/ (default:./result_path_opt/)
├─ final_geometries_trj.xyz # XYZ trajectory with per-image energies in the comment line
├─ final_geometries.pdb # Same as_trj.xyz but mapped back to the reference PDB ordering
├─ hei.xyz # Highest-energy image (XYZ, always written)
├─ hei.pdb # HEI in PDB format (when reference PDB is available)
├─ align_refine/ # External alignment/refinement artifacts
├─ preopt/ # Endpoint pre-optimization outputs (present when --preopt)
└─ <optimizer dumps> # Present when --dump or opt.dump_restart > 0
```

Merge order is **defaults < config < explicit CLI < override**.
### Section `geom`
- `coord_type`: Coordinate type (cartesian vs dlc internals).
- `freeze_atoms`: 0-based frozen atoms merged with CLI `--freeze-atoms`.

### Section `calc` / `mlmm`
- ML/MM calculator setup: `charge`, `spin`, UMA `model`, `task_name`, `device`, neighbor radii, Hessian options, etc.

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

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path-search](path_search.md) -- Recursive MEP search with automatic refinement (for 2+ structures)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow (uses path-search by default)
- [YAML Reference](yaml_reference.md) -- Full `gs`, `opt` configuration options
