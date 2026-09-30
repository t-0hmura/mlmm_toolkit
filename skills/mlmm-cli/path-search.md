# `mlmm path-search`

## Purpose

Recursive minimum-energy-path (MEP) search across two or more
endpoints. Detects bond changes along the candidate MEP and
**recursively re-segments** the path into candidate reaction intervals.
Output: flat per-segment files
(`mep_seg_NN_trj.xyz`, `hei_seg_NN.{xyz,pdb}`),
plus a stitched `mep_trj.pdb`/`mep_trj.xyz`, a `mep_trj.cif` companion for bridged
input, and energy diagrams.

`mlmm all --refine-path` selects this engine. Validate each HEI with TS/IRC
before treating its segment as an elementary step.

## Synopsis

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb [-i 1.R.pdb 2.IM.pdb 3.P.pdb] --parm7 real.parm7 \
    [--mep-mode gsm|dmf] [--refine-mode peak|minima] \
    [--max-nodes 20] [-l 'RES:Q,...'] [-b uma|orb|mace|aimnet2|dft] \
    [-o ./result_path_search/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below, **`mlmm-toolkit` requires an
Amber topology** and supports layer-aware selection. Most subcommands
accept:

| flag | purpose |
|---|---|
| `--parm7 FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | **ML-region** charge (not whole-system); stored as `model_charge` |
| `-l, --ligand-charge` | Per-residue ML-region charge mapping |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required (≥ 2) | Two or more endpoints in reaction order |
| `--mep-mode` | str | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--refine-mode` | str | mode-dep | `peak` (HEI±1) or `minima` (nearest local minima) |
| `--max-nodes` | int | 20 | Max internal nodes per segment string |
| `--max-depth` | int | 10 | Recursive subdivision levels; `0` disables it. A capped interval is tagged `seg_NNN_maxdepth` |
| `--thresh` | str | `gau` | Single-structure optimization and input-alignment convergence preset |
| `--thresh-gsm` | str | `gau_loose` | GSM string-optimizer convergence preset |
| `--dmf-tol` | str/float | `tight` | DMF IPOPT dual-infeasibility tolerance: `tight`, `middle`, `loose`, or a positive float |
| `-q, --charge` / `-l` / `-m` | — | — | Charge / multiplicity (see common conventions) |
| `--opt-mode` | str | `grad` | Single-structure optimizer: L-BFGS (`grad`) or RFO (`hess`) |
| `-b, --backend` | str | `uma` | High-level backend (MLIP or optional DFT) |
| `-o, --out-dir` | path | `./result_path_search/` | Output directory |
| `--config` / `--show-config` / `--dry-run` | — | — | YAML config + preview |

## Examples

### 2-endpoint MEP, GSM, default refinement

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb --parm7 real.parm7 \
    -l 'SAM:1,GPP:-3' -b uma \
    -o result_path_search
```

### 3-endpoint with explicit intermediate

```bash
mlmm path-search -i 1.R.pdb 2.IM.pdb 3.P.pdb --parm7 real.parm7 \
    -l 'SAM:1,GPP:-3' -b uma --max-nodes 30 \
    -o result_path_search
```

### DMF mode (sometimes better for ill-conditioned strings)

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb --parm7 real.parm7 \
    --mep-mode dmf --refine-mode minima \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search
```

## Output

```
result_path_search/
├── summary.json                       # full result, see below
├── summary.log                        # human-readable
├── mep_seg_NN_trj.xyz                 # stitched string nodes per elementary step (NN = summary.json `index`, 1-based, 2 digits; not the 3-digit `tag`)
├── mep_seg_NN.pdb                     # PDB conversion when input is PDB or --ref-pdb supplied
├── hei_seg_NN.{xyz,pdb}               # highest-energy image (TS candidate) per segment
├── mep_trj.xyz                        # full stitched MEP across all segments
└── energy_diagram_*.png
```

The nested copies of the post-tsopt/IRC structures (`segments/seg_NN/structures/`)
are produced only by `mlmm all`, not by standalone `path-search`.

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
- Increasing `--max-nodes` trades cost for path resolution but does not repair
  chemically inconsistent endpoints. Benchmark convergence on the actual
  system and inspect the trajectory and bond changes.
- Output **does not** include refined TSs; `all` writes the optimized TS to
  `segments/seg_NN/ts.pdb` (canonical), and `ts/` holds the tsopt working output.

## See also

- `path-opt.md` — single-segment MEP optimization (the building block).
- `tsopt.md` — starts from each `hei_seg_NN.{xyz,pdb}` (TS candidate).
- `bond-summary.md` — same bond-change algorithm used here, standalone.
- Defaults: `import mlmm.core.defaults as d; print(d.SEARCH_KW, d.GS_KW, d.DMF_KW)`
