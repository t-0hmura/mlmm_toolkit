# `mlmm path-opt`

## Purpose

MEP optimization for **one** segment between two endpoints. The
building block of `path-search` (which runs `path-opt` internally
once per segment, then bond-segments any multi-step paths). Use it
standalone to refine one segment without re-running the whole
recursive search.

## Synopsis

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm real.parm7 \
    [--mep-mode gsm|dmf] [--max-nodes 20] \
    [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_opt/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below, **`mlmm-toolkit` requires an
Amber topology** and supports layer-aware selection. Most subcommands
accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | **ML-region** charge (not whole-system); stored as `model_charge` |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (= 2) | Reactant and product, identical atom ordering |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--max-nodes` | int | 20 | Max internal nodes (final string ≤ `max-nodes + 2`) |
| `--thresh` | str | `gau` | Endpoint preoptimization convergence preset |
| `--thresh-gsm` | str | `gau_loose` | GSM string-optimizer convergence preset |
| `--thresh-dmf` | str/float | `tight` | DMF IPOPT dual-infeasibility tolerance: `tight`, `middle`, `loose`, or a positive float |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_path_opt/` | Output directory |

## Examples

### Default GSM single segment

```bash
mlmm path-opt -i R.pdb P.pdb --parm real.parm7 -q 0 -m 1 -b uma -o result_path_opt
```

### DMF for hard-to-converge strings

```bash
mlmm path-opt -i R.pdb P.pdb --parm real.parm7 -l 'GPP:-3' --mep-mode dmf -b mace \
    -o result_path_opt_dmf
```

## Output

```
result_path_opt/
├── result.json                     # written when --out-json
├── final_geometries_trj.xyz        # converged string trajectory
├── final_geometries.pdb            # PDB conversion when input is PDB
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
- If one MEP optimizer stalls, compare GSM and DMF on the actual
  backend/model/system and inspect the path before changing `--max-nodes`.
  System size alone is not a portable optimizer-selection rule.

## See also

- `path-search.md` — the recursive driver around this command.
- `opt.md` — pre-relax endpoints before path-opt.
- Defaults: `import mlmm.core.defaults as d; print(d.GS_KW, d.DMF_KW, d.STOPT_KW)`
