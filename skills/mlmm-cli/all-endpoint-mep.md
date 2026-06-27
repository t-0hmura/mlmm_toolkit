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

`--mep-mode` / `--refine-mode` are exposed on standalone `path-search`,
not on `mlmm all`; by default `mlmm all` runs single-pass `path-opt`
(GSM) between adjacent pairs, and recursive `path-search` runs only with
`--refine-path`.

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

- `path_opt/mep.pdb` — the full MEP across all segments
  (`path_search/...` under `--refine-path`)
- `path_opt/seg_01/ … seg_NN/` — per-segment string of nodes
- `seg_NN/{reactant,ts,product}.pdb` — canonical R/TS/P per segment
  after IRC + LBFGS endpoint optimization
- `summary.json["segments"]` — list of `{index, barrier_kcal,
  delta_kcal, bond_changes, ...}` entries

## Distinctive failure modes

| Symptom in `summary.json` | Likely cause | Fix |
|---|---|---|
| `path_search.status == "wrong_reaction"` | Bond-change detector found extra changes; the reaction in the inputs and the reaction the optimizer found don't match. | Check which bonds changed via `bond-summary -i 1.R.pdb 3.P.pdb`; rerun standalone `path-search` with `--refine-mode minima`, or supply IM explicitly. |
| `tsopt.n_imaginary > 1` for a segment | Multi-imaginary-mode TS — common with Orb on tricky systems | Re-run that segment with `mlmm tsopt --opt-mode rsirfo -b uma` (RS-I-RFO is more robust than Dimer for ill-conditioned Hessians). |
| Different atom counts across `-i` inputs | Inconsistent extractions | Re-extract per the snippet above, verify with `wc -l 1.R.pdb 3.P.pdb`. |

## Caveats

- The MEP method used internally by `mlmm all` is GSM; to switch to
  DMF, drive `mlmm path-search` directly with `--mep-mode dmf` and
  then feed the result into `mlmm tsopt`/`freq`/`irc` manually.
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
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | Force total system charge (highest priority over derived charges) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.
