# `mlmm path-search`

## Purpose

Recursive minimum-energy-path (MEP) search across two or more
endpoints. Detects bond changes along the candidate MEP and
**recursively re-segments** the path until each segment crosses
exactly one transition state. Output: one `seg_NN/` per elementary
step, plus a stitched `mep.pdb` and energy diagrams.

This is the engine behind `mlmm all` in endpoint-MEP mode.

## Synopsis

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb [-i 1.R.pdb 2.IM.pdb 3.P.pdb] --parm real.parm7 \
    [--mep-mode gsm|dmf] [--refine-mode peak|minima] \
    [--max-nodes 20] [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2] \
    [-o ./result_path_search/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below, **`mlmm-toolkit` requires an
Amber topology** and supports layer-aware selection. Most subcommands
accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); used only when `--model-pdb` is omitted (`--model-pdb` takes precedence) |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system); stored as `model_charge` |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (≥ 2) | Two or more endpoints in reaction order |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--refine-mode` | str | mode-dep | `peak` (HEI±1) or `minima` (nearest local minima) |
| `--max-nodes` | int | 20 | Max internal nodes per segment string |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / multiplicity (see common conventions) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_path_search/` | Output directory |
| `--config` / `--show-config` / `--dry-run` | — | — | YAML config + preview |

## Examples

### 2-endpoint MEP, GSM, default refinement

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma \
    -o result_path_search
```

### 3-endpoint with explicit intermediate

```bash
mlmm path-search -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -l 'SAM:1,GPP:-3' -b uma --max-nodes 30 \
    -o result_path_search
```

### DMF mode (sometimes better for ill-conditioned strings)

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb \
    --mep-mode dmf --refine-mode minima \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search
```

## Output

```
result_path_search/
├── summary.json                       # full result, see below
├── summary.log                        # human-readable
├── mep_seg_NN_trj.xyz                 # stitched string nodes per elementary step (2-digit tag)
├── mep_seg_NN.pdb                     # PDB conversion when input is PDB or --ref-pdb supplied
├── hei_seg_NN.{xyz,pdb}               # highest-energy image (TS candidate) per segment
├── mep_trj.xyz                        # full stitched MEP across all segments
└── energy_diagram_*.png
```

The per-segment TS refinement seeds (`segments/seg_NN/structures/`) are
produced only by `mlmm all`, not by standalone `path-search`.

`summary.json["segments"]` lists each elementary step with:

```python
{
  "index": 1,
  "tag": "seg_000_refine",
  "kind": "seg",
  "barrier_kcal": 21.5,
  "delta_kcal": -0.7,
  "bond_changes": "...summary text..."
}
```

## Caveats

- All `-i` inputs must have identical atom counts and ordering.
- Recursive segmentation can produce **more** segments than `len(-i) - 1`
  — that's the whole point: it finds intermediates the user didn't
  supply.
- `--max-nodes` bigger than 30 rarely helps; if a segment doesn't
  converge with 20 nodes, the chemistry is usually the problem.
- Output **does not** include refined TSs; those are `segments/seg_NN/tsopt/`
  in the `all` pipeline.

## See also

- `path-opt.md` — single-segment MEP optimization (the building block).
- `tsopt.md` — runs after path-search on each `segments/seg_NN/structures/ts.pdb`.
- `bond-summary.md` — same bond-change algorithm used here, standalone.
- Defaults: `import mlmm.core.defaults as d; print(d.SEARCH_KW, d.GS_KW, d.DMF_KW)`