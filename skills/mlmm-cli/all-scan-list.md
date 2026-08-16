# `mlmm all` — scan-list mode

## When to use

You have **only the reactant** (no product structure) and you can
articulate the chemistry as a sequence of staged distance scans —
e.g. "first push the methyl from S of SAM to C7 of GPP, then snap H11
to OE2 of GLU186". `mlmm all` runs each stage in order, then ties
the resulting trajectories into an MEP with single-pass `path-opt`. With
`--refine-path`, the recursive bond-change segmentation inserts any
intermediates it finds.

Typical use: multistep methyltransferase mechanisms where the user
encodes successive distance scans (methyl transfer → proton abstraction,
etc.) as separate stages.

## Synopsis

```bash
mlmm all --parm enzyme.parm7 -i 1.R.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --scan-lists \
        '[("CS1 SAM 320","GPP 321 C7",1.60)]' \
        '[("GPP`321/H11","GLU`186/OE2",0.90)]' \
    --tsopt --thermo \
    -o result_scan
```

Each literal after the single `--scan-lists` flag is **one stage**; do not
repeat the flag. Stages run
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
| `"RESNAME RESID NAME"` | Atom by residue name + residue index + PDB name, separated by single spaces |
| `"RESNAME\`RESID/NAME"` | Compact form with backticks and slash; same three fields, different separators |
| `"CHAIN:RESNAME:RESID[ICODE]:NAME"` | Exact chain-qualified form for repeated or mmCIF identifiers |

All bonds in a stage are driven simultaneously. If you want them done
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

After scans complete, `mlmm all` stitches the scan trajectories with
single-pass `path-opt` (GSM) by default; pass `--refine-path` to run the
recursive `path-search` instead. `--mep-mode dmf` selects DMF for either
route.

Unlike endpoint-MEP mode, `-i` is **a single PDB** (the reactant). The
toolkit synthesizes intermediate / product geometries from the scan
trajectories.

## Output

Same overall tree as in `all.md`, plus per-stage scan output:

```
result_scan/
├── mep.pdb / mep.cif / mep_trj.xyz # CIF companion for bridged input
├── segments/
│   └── seg_NN/                     # canonical R/TS/P + post-processing per segment
└── _work/                          # pipeline scratch
    ├── scan/
    │   ├── stage_01/  scan_*.xyz   # raw distance-restraint scan trajectory
    │   ├── stage_02/  scan_*.xyz
    │   └── ...
    └── path_opt/                   # raw MEP-engine output (path_search/ with --refine-path)
        └── seg_NN_mep/             # one MEP per stitched pair (recursive bond-change splitting only with --refine-path)
```

The `all` pipeline runs the scan in `_work/scan/` and does **not** emit a JSON
record (the top-level `summary.json` is the `all` envelope and carries no scan
stages). To get the stage-by-stage record as JSON, run the scan standalone with
`--out-json`; its `summary.json` then holds the record under the top-level
`stages` key:

```python
import json
d = json.load(open("result_scan/summary.json"))  # from `mlmm scan ... --out-json`
for stage in d["stages"]:
    print(stage["index"], stage["converged"], stage["bond_changes"], stage["target_distances_angstrom"])
```

## Distinctive failure modes

| Symptom | Cause | Fix |
|---|---|---|
| Stage k goes to a different geometry than expected | Distance restraint not strong enough; SCF found a side product | Tighten the target distance, or split a complex stage into two simpler ones |
| `--scan-lists` triggers a Python literal-eval error | Quoting mistake | Wrap each stage in single quotes outside, double quotes inside; backticks survive bash without escaping |
| Path search reports more segments than expected | Bond-change detector found a "free" intermediate | This is usually correct; check the IM geometry in `seg_01/product.pdb` (= `seg_02/reactant.pdb`) |

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
- Defaults: `import mlmm.core.defaults as d; print(d.SEARCH_KW, d.STOPT_KW)`.
## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — optional; when omitted, `mm_parm` generates a parm7 from the input PDB |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over extraction- or B-factor-derived ML membership |
| `--detect-layer` | Automatically read valid B-factor MM sublayers; without explicit or extraction-derived ML membership, B-factors also define ML membership. Enabled by default. |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | Override the net ML-region/model charge (highest priority) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.
