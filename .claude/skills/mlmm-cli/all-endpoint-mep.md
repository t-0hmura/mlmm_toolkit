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
mlmm all -i 1.R.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    [--dft --func-basis 'wb97m-v/def2-svp'] \
    -o result_mep
```

For a known multistep mechanism, supply each intermediate explicitly:

```bash
mlmm all -i 1.R.pdb 2.IM.pdb 3.P.pdb \
    -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo \
    -o result_mep_3pt
```

The recursive bond-change segmentation in `path-search` will still run
inside each adjacent pair, so you don't have to provide every
elementary step — just the "obvious" ones from the literature.

## Mode-specific flags

| Flag | Default | Meaning |
|---|---|---|
| `--mep-mode` | `gsm` | `gsm` (Growing String) or `dmf` (Direct Max Flux) |
| `--refine-mode` | mode-dependent | `peak` (HEI±1 around max) or `minima` (nearest local minima) |
| `--max-nodes` | 20 | Maximum string nodes per segment (final string ≤ `max-nodes + 2`) |

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

- `path_search/mep.pdb` — the full MEP across all segments
- `path_search/seg_01/ … seg_NN/` — per-segment string of nodes
- `seg_NN/{reactant,ts,product}.pdb` — canonical R/TS/P per segment
  after IRC + LBFGS endpoint optimization
- `summary.json["segments"]` — list of `{index, barrier_kcal,
  delta_kcal, bond_changes, ...}` entries

## Distinctive failure modes

| Symptom in `summary.json` | Likely cause | Fix |
|---|---|---|
| `path_search.status == "wrong_reaction"` | Bond-change detector found extra changes; the reaction in the inputs and the reaction the optimizer found don't match. | Check which bonds changed via `bond-summary -i 1.R.pdb 3.P.pdb`; consider `--refine-mode minima` or supply IM explicitly. |
| `tsopt.n_imaginary > 1` for a segment | Multi-imaginary-mode TS — common with Orb on tricky systems | Re-run that segment with `mlmm tsopt --opt-mode rsirfo -b uma` (RS-I-RFO is more robust than Dimer for ill-conditioned Hessians). |
| Different atom counts across `-i` inputs | Inconsistent extractions | Re-extract per the snippet above, verify with `wc -l 1.R.pdb 3.P.pdb`. |

## Caveats

- The "right" choice between `--mep-mode gsm` and `--mep-mode dmf`
  depends on system size and how much you trust the initial MEP guess;
  GSM is the safer default.
- Path search may discover **more** segments than you have inputs:
  if `summary.json["n_segments"] > len(inputs) - 1`, that's the
  recursive bond-change segmentation finding intermediates the inputs
  didn't contain. Often this is the *correct* answer.

## See also

- `all.md` — base orientation (output tree, summary.json schema).
- `path-search.md` — recursive MEP search internals.
- `bond-summary.md` — what bond-change detection looks like.
- `mlmm-workflows-output/SKILL.md` — interpreting multi-segment
  results.
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
