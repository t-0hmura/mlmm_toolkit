# `mlmm path-opt`

## Purpose

MEP optimization for **one** segment between two endpoints. The
building block of `path-search` (which runs `path-opt` internally
once per segment, then bond-segments any multi-step paths). Use it
standalone to refine one segment without re-running the whole
recursive search.

## Synopsis

```bash
mlmm path-opt -i reactant.pdb product.pdb \
    [--mep-mode gsm|dmf] [--max-nodes 20] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_opt/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

Beyond the cluster-style flags below (inherited from `pdb2reaction`),
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); overrides `--model-pdb` |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (= 2) | Reactant and product, identical atom ordering |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--max-nodes` | int | 20 | Max internal nodes (final string ≤ `max-nodes + 2`) |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent name |
| `-o, --out-dir` | path | `./result_path_opt/` | Output directory |

## Examples

### Default GSM single segment

```bash
mlmm path-opt -i R.xyz P.xyz -q 0 -m 1 -b uma -o result_path_opt
```

### DMF for tough strings

```bash
mlmm path-opt -i R.pdb P.pdb -l 'GPP:-3' --mep-mode dmf -b mace \
    -o result_path_opt_dmf
```

## Output

```
result_path_opt/
├── result.json
├── final_string_geoms.xyz          # converged string trajectory
├── nodes/                          # per-node geometries
└── hei.{xyz,pdb}                   # highest-energy image (TS candidate)
```

`result.json` reports converged string energies, gradient norm, and
the path-opt status (`converged` / `not_converged`).

## When to use vs path-search

- **`path-search`** if you want recursive bond-change segmentation
  for a possibly multi-step mechanism.
- **`path-opt`** if you already know the segment is a single-step
  reaction and just want the MEP between two endpoints, without the
  segmentation overhead.

## Caveats

- Convergence is sensitive to initial endpoint geometries. If
  `not_converged`, try running `mlmm opt` on each endpoint
  first to pre-relax to local minima.
- For systems with > 200 atoms, GSM's per-node Hessian-free curvature
  estimation may stall — try `--mep-mode dmf` or relax `--max-nodes`.

## See also

- `path-search.md` — the recursive driver around this command.
- `opt.md` — pre-relax endpoints before path-opt.
- Defaults: `import mlmm.defaults as d; print(d.GS_KW, d.DMF_KW, d.STOPT_KW)`