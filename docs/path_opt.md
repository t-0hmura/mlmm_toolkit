# `path-opt`

## Overview

> **Summary:** Find an MEP between **exactly two** enzyme structures using the Growing String Method (GSM) with the ML/MM calculator. Writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate.

### Quick reference
- **Use when:** You have reactant and product endpoints (R -> P) and want a first-pass MEP.
- **Method:** PySisyphus `GrowingString` with the ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`).
- **Outputs:** `final_geometries.trj` (path) and `gsm_hei.xyz` (HEI), plus optional `.pdb` companions.
- **Defaults:** `--climb True`, `--max-nodes 10`, `--max-cycles 300`.
- **Next step:** Validate the HEI with `tsopt` -> `freq` (expect one imaginary mode) -> `irc`.

`mlmm path-opt` optimizes a minimum-energy path between two enzyme states using PySisyphus `GrowingString` with the ML/MM calculator. The ML/MM calculator keeps the full enzyme complex without link atoms: the ML region is defined by `--model-pdb`, the Amber topology comes from `--real-parm7`, and both endpoints are supplied as PDBs containing the full system coordinates.

Configuration comes from YAML (`geom`, `calc`/`mlmm`, `gs`, `opt`) with precedence **CLI > YAML > defaults**. The CLI always overrides the ML/MM input files, charge, and multiplicity. `StringOptimizer.align` remains disabled; instead, an external Kabsch-based alignment/refinement routine aligns all images prior to path growth (respecting any frozen atoms from YAML/CLI).

## Usage
```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
              -q CHARGE [-m MULT] [options]
```

### Examples
```bash
# Minimal invocation
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# With frozen atoms, more nodes, and YAML overrides
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --freeze-atoms "1,3,5,7" --max-nodes 10 --max-cycles 200 --dump True --out-dir ./result_path_opt/ \
    --args-yaml ./args.yaml

# With endpoint pre-optimization
mlmm path-opt -i reac.pdb prod.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 \
    --preopt True --preopt-max-cycles 20000
```

## Workflow
1. **Load endpoints** -- Read both PDB structures and resolve charge/spin from CLI or defaults.
   Set up the ML/MM calculator with `--real-parm7`, `--model-pdb`, and charge/spin.
2. **Pre-alignment** -- All endpoints after the first are Kabsch-aligned to the first
   structure. If `freeze_atoms` is defined, only those atoms participate in the RMSD
   fit; the resulting transform is applied to all atoms.
3. **Optional pre-optimization** -- With `--preopt True`, each endpoint is pre-optimized
   by LBFGS (using the same ML/MM calculator) before alignment and string growth.
   The number of LBFGS cycles is controlled by `--preopt-max-cycles` (default: 10000).
4. **String growth** -- PySisyphus `GrowingString` grows the path between the endpoints
   using `(max_nodes + 2)` images including endpoints.
5. **Climbing image** -- With `--climb True`, a climbing-image refinement is applied after
   the string grows fully, and the highest-energy image (HEI) is reported.
6. **Output** -- Final path trajectory and HEI are written as XYZ and PDB files.
   PDB conversion is performed when the inputs are PDBs.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product PDB structures. | Required |
| `--real-parm7 PATH` | Amber prmtop for the full REAL system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs). | Required |
| `-q, --charge INT` | Total ML-region charge. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze (converted to 0-based; merged with YAML `geom.freeze_atoms`). | _None_ |
| `--max-nodes INT` | Number of internal string nodes (total images = `max_nodes + 2`). | `10` |
| `--max-cycles INT` | Optimizer macro-iteration cap (growth + refinement). Also sets `opt.stop_in_when_full`. | `300` |
| `--climb {True\|False}` | Enable climbing-image refinement after full string growth. | `True` |
| `--preopt {True\|False}` | Pre-optimize each endpoint with LBFGS before alignment/string growth. | `False` |
| `--preopt-max-cycles INT` | Cap for endpoint pre-optimization cycles. | `10000` |
| `--thresh TEXT` | Convergence preset override (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--dump {True\|False}` | Dump optimizer trajectories and restarts inside `out_dir`. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--args-yaml FILE` | YAML overrides (sections `geom`, `calc`/`mlmm`, `gs`, `opt`). | _None_ |

## Outputs
```
out_dir/ (default: ./result_path_opt/)
├─ final_geometries.trj        # XYZ trajectory with per-image energies in the comment line
├─ final_geometries.pdb        # Same as .trj but mapped back to the reference PDB ordering
├─ gsm_hei.xyz                 # Highest-energy image (XYZ, always written)
├─ gsm_hei.pdb                 # HEI in PDB format (when reference PDB is available)
├─ align_refine/               # External alignment/refinement artifacts
├─ preopt/                     # Endpoint pre-optimization outputs (present when --preopt True)
└─ <optimizer dumps>           # Present when --dump True or opt.dump_restart > 0
```

## YAML configuration (`--args-yaml`)
YAML parameters override internal defaults; CLI overrides YAML.

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

- [path-search](path_search.md) -- Recursive MEP search with automatic refinement (for 2+ structures)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow (uses path-search by default)
- [YAML Reference](yaml-reference.md) -- Full `gs`, `opt` configuration options
