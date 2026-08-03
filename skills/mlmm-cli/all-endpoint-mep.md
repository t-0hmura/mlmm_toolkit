# `mlmm all` — endpoint-MEP mode

## When to use

You have **two or more reaction-ordered structures** (reactant, optional
intermediate(s), product), all with the **same atom count and atom
ordering**. The pipeline interpolates an MEP between adjacent
structures and segments multi-step paths automatically.

This is the most common mode for a published-mechanism reproduction
where you have R and P (and sometimes IM) coordinates from a prior QM
or QM/MM study.

## Synopsis

```bash
mlmm all --parm enzyme.parm7 -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    [--dft --dft-func-basis 'wb97m-v/def2-svp'] \
    -o result_mep
```

For a known multistep mechanism, supply each intermediate explicitly:

```bash
mlmm all --parm enzyme.parm7 -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep_3pt
```

By default each adjacent pair is connected with a single-pass `path-opt`,
so the endpoints you pass are taken as the elementary steps. Add
`--refine-path` to enable the recursive bond-change segmentation in
`path-search`, which splits a pair further when it detects intermediate
bond changes — then you don't have to provide every elementary step,
just the "obvious" ones from the literature.

## Mode-specific flags

| Flag | Default | Meaning |
|---|---|---|
| `--max-nodes` | 20 | Maximum string nodes per segment (final string ≤ `max-nodes + 2`) |
| `--mep-mode gsm\|dmf` | `gsm` | MEP optimizer for both the default and recursive path routes |
| `--dmf-backend gpu\|cpu` | `gpu` | DMF implementation; set `cpu` after a GPU out-of-memory error |

By default `mlmm all` runs single-pass `path-opt` between adjacent pairs;
`--refine-path` selects recursive `path-search`. `--mep-mode` controls the
optimizer in either route. The finer-grained `--refine-mode` remains a
standalone `path-search` option.

`--scan-lists` is **not** allowed in this mode — it triggers
`all-scan-list.md` instead.

## Atom-count consistency requirement

All `-i` inputs must have:

- the same number of atoms,
- the same element sequence (atom ordering),
- the same residue assignments.

If the inputs come from different programs or were re-numbered, run
them through `extract` once to canonicalize ordering:

```bash
mlmm extract -i 1.R_raw.pdb 3.P_raw.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o "1.R.pdb" "3.P.pdb"
```

## Output

Same as the base `all.md`. Specifically for endpoint-MEP mode:

- `mep.pdb` (and `mep.cif` for bridged input) at the output root — the full MEP across all segments
  (raw engine copy under `_work/path_opt/`, or `_work/path_search/` with `--refine-path`)
- `segments/seg_01/ … seg_NN/` — per-segment string of nodes
- `segments/seg_NN/{reactant,ts,product}.pdb` (plus CIF companions for bridged input) — canonical R/TS/P per
  segment after IRC + RFO endpoint optimization (`--opt-mode-post hess` default; `grad` selects L-BFGS)
- `summary.json["segments"]` — list of `{index, barrier_kcal,
  delta_kcal, bond_changes, ...}` entries

## Distinctive failure modes

| Symptom in `summary.json` | Likely cause | Fix |
|---|---|---|
| `status == "partial"`, or `bond-summary` reports extra changes vs the optimized MEP | Bond-change detector found extra changes; the reaction in the inputs and the reaction the optimizer found don't match. | Check which bonds changed via `bond-summary -i 1.R.pdb 3.P.pdb`; rerun standalone `path-search` with `--refine-mode minima`, or supply IM explicitly. |
| `tsopt.n_imaginary_modes > 1` for a segment | Higher-order saddle or unresolved soft modes | Compare a Hessian-based mode and Dimer on the same seed/backend, then rerun frequency analysis and IRC connectivity checks. |
| Different atom counts across `-i` inputs | Inconsistent extractions | Re-extract per the snippet above, verify with `wc -l 1.R.pdb 3.P.pdb`. |

## Caveats

- GSM is the default. Use `--mep-mode dmf`; choose `--dmf-backend cpu`
  when the GPU implementation runs out of memory.
- Under `--refine-path`, path search may discover **more** segments than
  you have inputs: if `summary.json["n_segments"] > len(inputs) - 1`,
  that's the recursive bond-change segmentation finding intermediates the
  inputs didn't contain — often the *correct* answer. (Default single-pass
  `path-opt` yields one segment per adjacent input pair.)

## See also

- `all.md` — base orientation (output tree, summary.json schema).
- `path-search.md` — recursive MEP search internals.
- `bond-summary.md` — what bond-change detection looks like.
- `mlmm-workflows-output/SKILL.md` — interpreting multi-segment
  results.
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
| `--embedcharge / --no-embedcharge` | Unavailable in v0.3.3; use `--no-embedcharge` |
| `-q, --charge` | Override the net ML-region/model charge (highest priority) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.
