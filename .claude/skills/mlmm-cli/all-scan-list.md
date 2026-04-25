# `mlmm all` — scan-list mode

## When to use

You have **only the reactant** (no product structure) and you can
articulate the chemistry as a sequence of staged distance scans —
e.g. "first push the methyl from S of SAM to C7 of GPP, then snap H11
to OE2 of GLU186". `mlmm all` runs each stage in order, ties
the resulting trajectories into an MEP, and the recursive bond-change
segmentation slots in any intermediates it finds.

This is how the bezA case study in the published benchmark was driven.

## Synopsis

```bash
mlmm all -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists \
        '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
        '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each `--scan-lists '...'` argument is **one stage**. Stages run
sequentially; the final geometry of stage *k* is the input geometry of
stage *k+1*.

## `--scan-lists` syntax

Each argument is a Python literal-eval expression: a list of bond
tuples, where each tuple is `(atom_a, atom_b, target_distance_Å)`.

```
[ ("<atom-spec>", "<atom-spec>", <float>) , ... ]
```

`<atom-spec>` formats:

| Form | Meaning |
|---|---|
| `"NAME RESNAME RESID"` | Atom by PDB name + residue name + residue index, separated by single spaces |
| `"RESNAME\`RESID/NAME"` | Compact "tilde-form" with backticks and slash; same fields, different separators |
| `"chainID:RESID:NAME"` | Chain-aware lookup |

Multiple bonds per stage drive simultaneously. If you want them done
**sequentially**, split them into separate `--scan-lists` arguments.

Examples:

```bash
# One stage, two bonds driven together (concerted SN2):
--scan-lists '[("CS1 SAM 320","GPP 321 C7",1.60),("GPP 321 C7","S SAM 320",3.0)]'

# Two stages, one bond each (stepwise mechanism):
--scan-lists '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
             '[("GPP`321/H11","GLU`186/OE2",0.90)]'
```

## Mode-specific flags

| Flag | Default | Meaning |
|---|---|---|
| `--scan-lists` | required | One or more stages of distance-restraint scans |
| `--scan-mode` | `staged` | (advanced) governs how nodes are seeded between stages |
| `--mep-mode` | `gsm` | After scans complete, MEP refinement uses GSM unless `dmf` |

Unlike endpoint-MEP mode, `-i` is **a single PDB** (the reactant). The
toolkit synthesizes intermediate / product geometries from the scan
trajectories.

## Output

Same overall tree as in `all.md`, plus per-stage scan output:

```
result_scan/
├── path_search/
│   ├── scan/
│   │   ├── stage_01/  scan_*.xyz   # raw distance-restraint scan trajectory
│   │   ├── stage_02/  scan_*.xyz
│   │   └── ...
│   ├── mep.pdb / mep_trj.xyz       # stitched MEP
│   ├── seg_NN/                     # segments after bond-change splitting
│   └── post_seg_NN/                # per-segment refinements
├── seg_NN/                         # canonical R/TS/P per segment (top-level)
└── summary.json
```

`summary.json["scan"]` (or under `path_search.scan`) carries the
stage-by-stage record:

```python
import json
d = json.load(open("result_scan/summary.json"))
for stage in d["path_search"]["scan"]["stages"]:
    print(stage["index"], stage["bonds"], stage["status"], stage["final_distance_A"])
```

## Distinctive failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Stage k goes to a different geometry than expected | Distance restraint not strong enough; SCF found a side product | Tighten the target distance, or split a complex stage into two simpler ones |
| `--scan-lists` triggers a Python literal-eval error | Quoting mistake | Wrap each stage in single quotes outside, double quotes inside; backticks survive bash without escaping |
| Path search reports more segments than expected | Bond-change detector found a "free" intermediate | This is usually correct; check the IM geometry in `seg_01/product.pdb` (= `seg_03/reactant.pdb`) |

## Caveats

- The atom specs must match the **exact** atom names in the input PDB
  (case sensitive). PyMOL/Maestro sometimes rename `CB` ↔ `CB1`.
- `--scan-lists` is incompatible with multiple `-i` inputs (the latter
  triggers `all-endpoint-mep.md`).
- Each stage can take longer than path-search itself; budget walltime
  accordingly.

## See also

- `all.md` — base orientation.
- `scan.md`, `scan2d.md`, `scan3d.md` — standalone distance scan
  subcommands (without the surrounding pipeline).
- `path-search.md` — what happens after all scans complete.
- Defaults: `import mlmm.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)`.
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
