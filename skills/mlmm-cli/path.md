# `mlmm path-opt` and `mlmm path-search`

Find a minimum-energy path (MEP) between layered full-system structures with
GSM (default) or DMF, and write its highest-energy image (HEI) as a TS
candidate. Run `mlmm path-opt -i R.pdb P.pdb --parm7 real.parm7 -q 0 -m 1` or
`mlmm path-search -i R.pdb IM.pdb P.pdb --parm7 real.parm7 -q 0 -m 1`.
Success is `hei.xyz` (path-opt) or one `hei_seg_NN.xyz` per reactive segment
in `summary.json` (path-search), each confirmed later by `tsopt` and `irc`.

## When to use

- `path-opt`: exactly two endpoints, one MEP in one pass, no segmentation.
  Use it when the step is known to be single, to redo one segment without
  rerunning the whole search, or to compare `--mep-mode gsm` and `dmf` on
  the same pair.
- `path-search`: two or more structures in reaction order. It finds where
  covalent bonds change and recursively splits the path into candidate
  reaction segments, so it suits a possibly multistep mechanism. It writes
  per-segment files, one stitched MEP, and energy diagrams.
- `all` runs `path-opt` for its MEP step; `all --refine-path` runs
  `path-search` instead.

## Minimal run

```bash
mlmm path-opt -i R.pdb P.pdb --parm7 real.parm7 -q 0 -m 1 -b uma -o result_path_opt
mlmm path-opt -i R.pdb P.pdb --parm7 real.parm7 -l 'GPP:-3' --mep-mode dmf -b mace \
    -o result_path_opt_dmf
```

```bash
mlmm path-search -i 1.R.pdb 3.P.pdb --parm7 real.parm7 \
    -l 'SAM:1,GPP:-3' -b uma -o result_path_search
mlmm path-search -i 1.R.pdb 2.IM.pdb 3.P.pdb --parm7 real.parm7 \
    -l 'SAM:1,GPP:-3' -b uma --max-nodes 30 -o result_path_search
mlmm path-search -i 1.R.pdb 3.P.pdb --parm7 real.parm7 \
    --mep-mode dmf --refine-mode minima -l 'SAM:1,GPP:-3' -b uma -o result_path_search
```

Give all structures after one `-i`. `--max-nodes` (default 20) sets the
movable images of each string, so a string has at most `max-nodes + 2`
images. `-b` selects `uma`, `orb`, `mace`, `aimnet2`, or `dft`.
`--refine-mode` is `peak` (HEI ± 1) for GSM and `minima` (nearest local
minima) for DMF unless given; `--max-depth` (default 10) limits the
recursion, and `0` turns it off.

## Judge success

**path-opt.** The console line `[write] Wrote '…/hei.xyz'` shows that the
TS candidate was written. With `--out-json`, `result.json` has `converged`,
`scientific_status` (`success` when the endpoint pre-optimization and the
MEP converged, otherwise `partial` or `failed`), `image_energies_hartree`,
`barrier_kcal` (HEI relative to the first image), `hei_index`, and
`n_images`. Files: `final_geometries_trj.xyz` (final string, energies on the
comment lines), `final_geometries.pdb`, `hei.xyz`, and `hei.pdb`. A HEI with
`hei_index` between 1 and `n_images − 2` is a TS candidate; at an endpoint,
no image lies above the higher endpoint and there is no candidate.

**path-search.** `summary.json` and `summary.log` are always written (no
`--out-json`); read `[2] Segment-level MEP summary` in the log.
`scientific_status` is `success` when the pre-optimizations and every path
run converged, otherwise `partial` or `failed`. Each entry of
`summary.json["segments"]` looks like:

```
{"index": 1, "tag": "seg_000_refine", "kind": "seg", "converged": true,
 "barrier_kcal": 21.5, "delta_kcal": -0.7, "bond_changes": "...summary text..."}
```

Files: `mep_trj.xyz` and `mep_trj.pdb` (the stitched MEP), `mep_plot.png`,
`energy_diagram_MEP.png`, and, for each segment with bond changes,
`mep_seg_NN_trj.xyz` (`mep_seg_NN.pdb`) and `hei_seg_NN.xyz`
(`hei_seg_NN.pdb`). NN is the segment `index`, counted from 01; NNN in a
`seg_NNN` tag counts the GSM/DMF runs from 000, so the two differ.

- A segment with bond changes and its `hei_seg_NN.xyz`: a TS candidate.
- A tag `seg_NNN_maxdepth` or `seg_NNN_kinklimit`: splitting stopped there,
  and the segment may hold more than one step. Raise `--max-depth` or give
  intermediates.
- Only `_kink` tags, or the warning `HEI is at an endpoint`: no bond change
  or no peak. Check the inputs or give intermediates.

A segment is a candidate, not a proven elementary step. Accept a HEI as a TS
only after `tsopt` gives n_imag = 1 and IRC reaches the intended R and P.

## Pitfalls and recovery

- Every input must have the same atoms in the same order as `--parm7`. For
  XYZ inputs, `path-opt` takes one `--ref-pdb` for both endpoints, and
  `path-search` takes one per input, in `-i` order.
- Convergence depends on the endpoints. Both are pre-optimized by default
  (`--preopt`, capped by `--preopt-max-cycles`); if the MEP still does not
  converge, relax each endpoint with [opt.md](opt.md) first.
- If one optimizer stalls, compare GSM and DMF on the real system and look
  at the path before changing `--max-nodes`; system size alone does not pick
  the optimizer. More nodes buy resolution, not repair: they cannot fix
  chemically inconsistent endpoints, so check the trajectory and the bond
  changes.
- `path-search` can return more segments than inputs minus one; it finds
  intermediates you did not give.
- DMF needs `cyipopt` and `pydmf`, which are not installed with
  mlmm-toolkit: `conda install -c conda-forge cyipopt -y`, then
  `pip install 'pydmf[torch]>=1.2'` (GPU) or `pip install 'pydmf>=1.2'`
  (CPU). The default `--dmf-backend gpu` stops when CUDA is missing; use
  `--dmf-backend cpu` then, or after a GPU out-of-memory error.
- DMF holds frozen atoms with a harmonic restraint, so they can drift a
  little; GSM keeps them fixed. DMF ignores `--climb` and `--fix-ends`.
- Neither command optimizes the TS. `all --tsopt` writes the optimized R,
  TS, and P under `segments/seg_NN/`.

## Next step

- [tsopt.md](tsopt.md): optimize `hei.pdb` or `hei_seg_NN.pdb` with the same
  `--parm7` and `--model-pdb`; for an `.xyz`, add `--ref-pdb`.
- [irc.md](irc.md): check that the TS connects the intended R and P.
- [opt.md](opt.md): pre-relax the endpoints.
- [bond-summary](utilities.md#bond-summary): the same bond-change rules on
  any two structures.
- [all.md](all.md): the full workflow around these commands.
