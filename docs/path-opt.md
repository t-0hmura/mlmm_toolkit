# `path-opt`

`mlmm path-opt` finds a minimum-energy path (MEP) between **exactly two** layered enzyme structures with GSM (default) or DMF (`--mep-mode dmf`), using the ML/MM calculator on the full enzyme complex. It writes the path trajectory and exports the highest-energy image (HEI) as a TS candidate. Use it when two layered endpoints are well defined and no intermediates are expected. It is the simpler MEP-only sibling of `path-search` (no recursive segmentation, no bond-change-driven decomposition). For workflows that start from **two or more** structures and automatically refine only the reactive region, use [path-search](path-search.md) instead.

## Examples

```bash
# Minimal invocation
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_path_opt
```

```bash
# Pre-optimize both endpoints before path growth
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --preopt --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

```bash
# Disable climbing-image refinement for a quick first pass
mlmm path-opt -i reac.pdb prod.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-climb --max-nodes 8 --out-dir ./result_path_opt_fast
# freeze selected atoms and keep optimizer dumps: --freeze-atoms "1,3,5,7" --dump
```

General command form:

```bash
mlmm path-opt -i REACTANT.pdb PRODUCT.pdb --parm real.parm7 --model-pdb model.pdb \
 -q CHARGE [-m MULT] [--mep-mode gsm|dmf] [--fix-ends/--no-fix-ends] [options]
```

`mlmm path-opt --help` shows core options; `mlmm path-opt --help-advanced` shows the full option list.

## Workflow
1. **Load endpoints** -- Read PDB/mmCIF structures, or XYZ coordinates with matching `--ref-pdb` topology, and resolve charge/spin.
    Set up the ML/MM calculator with `--parm`, `--model-pdb`, and charge/spin.
2. **Optional pre-optimization** -- With `--preopt`, each endpoint is pre-optimized
    by L-BFGS (using the same ML/MM calculator) before alignment and string growth.
    `--preopt-max-cycles` sets the L-BFGS cycle cap (default: 100000).
3. **Alignment and freeze-guided refinement** -- Endpoints after the first are rigidly
    aligned to the first. With `freeze_atoms`, the shared owner then performs its
    freeze-guided scan and L-BFGS relaxation toward the reference before string growth.
4. **Path optimization** -- `--mep-mode gsm` uses pysisyphus `GrowingString` with `(max_nodes + 2)` images including endpoints; `--mep-mode dmf` uses Direct Max Flux.
5. **Climbing image (GSM only)** -- With `--climb`, a climbing-image refinement is applied after string growth, and the highest-energy image (HEI) is reported.
6. **Output** -- Final path trajectory and HEI are written as XYZ files.
    PDB/CIF companions are written when conversion is enabled and a reference
    topology is available.

## Outputs

```text
out_dir/ (default: ./result_path_opt/)
├─ final_geometries_trj.xyz # XYZ trajectory with per-image energies in the comment line
├─ final_geometries.pdb # Same path as final_geometries_trj.xyz, mapped back to the reference PDB ordering
├─ hei.xyz # Highest-energy image (XYZ, always written)
├─ hei.pdb # HEI in PDB format (when reference PDB is available)
├─ align_refine/ # External alignment/refinement artifacts
├─ preopt/ # Endpoint pre-optimization outputs (present when --preopt)
└─ <optimizer dumps> # Present when --dump or opt.dump_restart > 0
```

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH PATH` | Reactant and product PDB/mmCIF structures, or XYZ coordinates with corresponding `--ref-pdb` entries. | Required |
| `--parm PATH` | Amber prmtop for the full REAL system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs). Optional when `--detect-layer` or `--model-indices` is used. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed like `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Automatically read B-factor layers (B=0/10/20). With explicit ML membership, only the MM sublayers are retained; otherwise B-factors also define ML membership. | Enabled |
| `-q, --charge INT` | Net ML-region charge. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-residue charge map, e.g. `SAM:1,PHN:-1`. Derives total charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--mep-mode [gsm\|dmf]` | MEP backend. | `gsm` |
| `--dmf-backend [cpu\|gpu]` | DMF compute backend (`--mep-mode dmf` only): `gpu` (`dmf.torch`/CUDA) or `cpu` (`dmf`/NumPy). Retry `cpu` on a GPU out-of-memory error. Requires `pydmf>=1.2`. | `gpu` |
| `--freeze-atoms TEXT` | Comma-separated 1-based atom indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. MM atoms beyond this are frozen. Providing `--movable-cutoff` disables `--detect-layer`. | _None_ |
| `--fix-ends/--no-fix-ends` | Fix endpoint structures during GSM growth (`gs.fix_first/fix_last`). | `True` |
| `--max-nodes INT` | Number of internal string nodes (total images = `max_nodes + 2`). | `20` |
| `--gsm-param [equi\|energy]` | GSM node parameterization after string growth. `energy` concentrates nodes in high-energy regions and may be tried when an equidistant path skips the reaction-coordinate region near the HEI; it does not identify a TS. | `equi` |
| `--max-cycles-gsm INT` | GSM string-optimizer cycle cap; also sets `stopt.stop_in_when_full`. | `300` |
| `--max-cycles-dmf INT` | DMF IPOPT iteration cap. | `3000` |
| `--climb/--no-climb` | Enable climbing-image refinement after full string growth. | `True` |
| `--preopt/--no-preopt` | Pre-optimize each endpoint with L-BFGS before alignment/string growth. | `True` |
| `--preopt-max-cycles INT` | Endpoint pre-optimization cycle cap. | `100000` |
| `--thresh TEXT` | Convergence preset override for endpoint pre-optimization only (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--thresh-gsm TEXT` | Convergence preset for the GSM string optimizer (`stopt.thresh`; same presets as `--thresh`). | `gau_loose` |
| `--thresh-dmf TEXT` | IPOPT dual-infeasibility tolerance of the DMF optimizer (`dmf.tol`): `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive float. Gaussian presets are rejected. | `tight` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend. Hessians use finite differences by default; set `calc.mm_fd: false` for the `hessian_ff` analytical path. | `hessian_ff` |
| `--dump/--no-dump` | Dump optimizer trajectories and restarts inside `out_dir`. | `False` |
| `-o, --out-dir TEXT` | Output directory. | `./result_path_opt/` |
| `--config FILE` | Base YAML configuration layer applied before explicit CLI values. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration (including YAML layers) and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print the execution plan without running optimization. Shown in `--help-advanced`. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--cmap/--no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## YAML configuration

Merge order is **defaults < config < explicit CLI**. The relevant sections are `geom` (`coord_type`, `freeze_atoms`), `calc` / `mlmm` (ML/MM calculator setup), `gs` (Growing String controls), and `opt` (StringOptimizer settings).

Full schema (every key and default): [YAML Reference](yaml-reference.md).

## Exit codes

| Code | Meaning |
| --- | --- |
| `0` | Success |
| `2` | CLI usage or configuration failure |
| `3` | Optimization failure |
| `4` | Final trajectory write error |
| `5` | HEI dump error |
| `130` | Keyboard interrupt |
| `1` | Unhandled exception |

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [path-search](path-search.md) — Recursive MEP search with automatic refinement (for 2+ structures)
- [opt](opt.md) — Single-structure geometry optimization
- [all](all.md) — End-to-end workflow (uses single-pass path-opt by default; add `--refine-path` for recursive path-search)
- [YAML Reference](yaml-reference.md) — Full `gs`, `opt` configuration options
