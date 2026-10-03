# mlmm ts-strategy

How to study a reaction mechanism with an ML/MM ONIOM model and get a correct barrier: precision, routes to a TS candidate, wrong n_imag, the barrier direction, splitting the reaction, a TS that does not come out, and controlled comparisons. Before any run, write down the bonds that form, the bonds that break, and every H atom that moves; this list becomes the `-s` coordinates and the points to check at the IRC endpoints, and every atom on it belongs in the ML region.

## 1. Precision: keep backend defaults

- Without `--precision`, UMA runs in fp32, ORB and MACE in fp64, and AIMNet2 in fp32 only, where an explicit fp64 is rejected. When precision matters, compare the supported precisions on the target system.
- `--precision fp32|fp64` is case-insensitive and is accepted by `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `all`.
- fp64 changes numerical precision only. `--deterministic` requests deterministic algorithms, but exact repeatability still has to be checked on the installed backend, model, and hardware ([`docs/backends.md`](../../docs/backends.md#precision)). AIMNet2 rejects `--deterministic` as well as fp64.

## 2. Two routes to a TS candidate

- **MEP search**, when you have R and P, optionally with intermediates in order: `path-search` splits the path recursively with GSM or DMF and returns an HEI per segment; `path-opt` optimizes one given segment and returns its HEI.
- **Restrained build-up**, when there is no usable second endpoint or TS guess: `scan` (or `scan2d` / `scan3d`) drives the reacting coordinates with a harmonic restraint E = ½ k (r − target)², k = 300 by default (`--restraint-k`), and relaxes the rest at each step.
- `opt` takes restraints only through `--distance-restraint`, with `--restraint-k` (also 300 by default). `scan` additionally drives staged targets up to a TS candidate.
- Pass a candidate from either route to `tsopt` and then `irc`, or use `all --tsopt`. The final PHVA (partial Hessian vibrational analysis) of `tsopt` gives the saddle order; add `freq` for the full modes or thermochemistry.

## 3. Wrong n_imag after TS optimization

A first-order saddle has exactly one imaginary mode, along the reaction; check its displacement and the IRC ends. Two or more fail regardless of their size.

- **Extra imaginary modes**: look at every mode's displacement, the MEP guess, the optimizer's stop reason, and the backend's numerical behavior. Then retry with another coordinate, flattening, or precision setting and read the new final PHVA.
- **n_imag = 0**: a failure, not a TS. Improve the MEP or the starting structure; `--flatten` only removes surplus modes and cannot create a missing reaction direction.
- **Poor MEP or HEI**: in `all`, try `--refine-path` before the TS optimization. It can split a poor path into several stages and costs more, so it is off by default.
- **Still no clean saddle**: revise the endpoints or the scan coordinates, and check that the single imaginary mode moves the reacting atoms.

Settings to retry with:

- `--flatten` runs the loop that removes extra imaginary modes, on `opt`, `tsopt`, and `all`; the TS commands run it when n_imag > 1. `--flatten` uses the iteration cap (50 by default) and `--no-flatten` sets it to zero.
- `--coord-type` takes `cart`, `redund`, `dlc`, or `tric` in `opt` and `tsopt`, and `cart` or `dlc` in `all`; the default is `cart`. Keep `cart` for ML/MM: `dlc` (delocalized internal coordinates) is slow to build for these models and its convergence depends on the system, so compare it with `cart` from the same start. `opt` accepts `dlc` with L-BFGS (`--opt-mode grad`) or RFO (`--opt-mode hess`). `path-opt` and `path-search` have no `--coord-type`; they read `geom.coord_type` from the YAML, `cart` or `dlc` only.
- Check any change of coordinate system with a frequency analysis and the IRC ends.
- `--ref-mode` is an advanced input, not a routine fix. It takes Cartesian 3N vectors in the input atom order (`.npz`, `.npy`, or text) and guides which negative Hessian root the TS optimizer follows; it does not replace the Hessian, and Dimer ignores it. `all` supplies the MEP tangent this way by default; with `all --no-tsopt-from-mep-tan`, `tsopt` picks its starting root from the Hessian of the starting structure.

Keep convergence and saddle order apart. A converged higher-order saddle is not a TS, although `all` still runs a diagnostic IRC, with a warning, when a valid negative root exists. No convergence, no imaginary mode, a failed or skipped final Hessian, or no valid negative root stops `all` before the IRC, with the TS files kept. The full list of remedies is in [Wrong imaginary-mode count after optimization](../../docs/tsopt.md#wrong-imaginary-mode-count-after-optimization).

## 4. Reading the barrier when the scan started from P

If the scan or path starts from P, the raw barrier is the reverse one.

- Forward barrier: E(TS) − E(R). Reverse barrier, the raw number of a P-start run: E(TS) − E(P). In `summary.json`, the barrier from the other end is `barrier_kcal − delta_kcal`.
- This is how you read the numbers, not a CLI flag. Before quoting a barrier, compare `segments/seg_NN/reactant.*` and `product.*` with the intended R and P instead of trusting the scan direction. In TS-only mode, R is the higher-energy IRC end ([outputs.md](outputs.md#oriented-rtsp-paths)).
- The top of a scan is only a TS candidate. Read the barrier from the energy after the TS optimization (`post_segments[].mlip.barrier_kcal`), counted from the minimum just before that stage.

## 5. Staged vs concerted scan

Each literal after `-s` is one stage. `all` and `scan` also accept `-s` repeated, but write `-s` once and list every literal after it.

- **Concerted**: one literal with several tuples; all coordinates move together in one stage. The mechanism need not be split up front, and `path-search` can split the path afterwards.
- **Staged**: several literals; each is one restrained relaxation in sequence, written to `stage_NN/`. Each stage needs its part of the mechanism defined up front, and each prescribed change appears as its own stage.

```bash
# Concerted (one stage, two coordinates driven together):
mlmm scan -i r.pdb --parm7 e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40),(7,9,1.60)]' -o result_concerted

# Staged (two sequential stages):
mlmm scan -i r.pdb --parm7 e.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' '[(7,9,0.95)]' -o result_staged
```

- Put every moving coordinate in its stage: the bond that breaks and each H that moves, not only the bond that forms.
- In `scan`, a 4-tuple `(i,j,low,high)` scans one distance both ways from the input and becomes two stages; in `all`, a 4-tuple is an angle target.

## 6. When the TS does not come out

A TS needs n_imag = 1 and IRC ends that are the intended R and P; an exit code of 0 is not that evidence.

- **Extra modes**: `--flatten` and `--refine-path` are in [Wrong n_imag after TS optimization](#3-wrong-n_imag-after-ts-optimization). After `--flatten`, check both IRC ends again, since n_imag can reach 1 on a TS of another reaction.
- **An extra mode that moves the ML region together with movable MM atoms**, such as a water moving with the substrate: with microiteration, the TS steps move only the ML atoms and the MM atoms bonded to them, and the other movable MM atoms only relax. Rerun `tsopt --no-microiter --flatten` so that RS-P-RFO steps over every movable atom. `all` has no microiteration switch, so pass the `--parm7` and `--model-pdb` of the `all` run to `tsopt`.
- **A narrowed Hessian**: if `tsopt` ran with a small `--hessian-cutoff`, widen it or leave it out.
- **Start from the extra mode**: `vib/imag_*_trj.xyz` and `vib/imag_*.pdb` hold 20 frames, and frames 6 and 16 are the largest displacements either way. Save each as a PDB and start TS-only mode or `tsopt` from each, with the same charge and multiplicity.
- **Start from another structure**: the HEI of the unsplit MEP (`hei.xyz` of `path-opt`), the HEI of each segment (`hei_seg_NN.xyz` of `path-search`; under `_work/path_opt/` in `all`, or `_work/path_search/` with `--refine-path`), or a scan frame near the top. If the forming and breaking bonds are both long at the HEI, or the HEI lies after the bond broke, the TS optimization tends to lose the reaction mode.
- **Path settings**: the number of internal images (`--max-nodes`, 20 by default), `--mep-mode dmf`, and `--gsm-param energy` all move the HEI.
- **An extra mode that also moves the reacting bonds**: two stages may overlap in one candidate; run them as separate stages.
- **Split the reaction differently**: run each plausible stage order and compare the energy diagrams. Run it stepwise and check, with `--scan-endopt`, whether the intermediate survives once the restraints are released; if it relaxes back, the concerted path fits better. When R and P can be prepared, give both to an MEP search instead of a scan.
- **ML region**: when a residue, water, or cofactor that takes part lies outside the ML region, enlarge it with a larger `-r` or a hand-built `--model-pdb` ([Enlarge when the model is too small](../mlmm-model-setup/SKILL.md#enlarge-when-the-model-is-too-small)), and recheck the charge.

The full guide is [`docs/mechanism-tips.md`](../../docs/mechanism-tips.md).

## 7. Controlled mutant-vs-WT comparison

- Form each barrier within one system (TS − R or TS − P), then compare the barriers: ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT. Do not subtract total energies of systems with different compositions.
- Give R and P of each system in MEP mode, so that R is the chemical reactant; G_TS − G_R is `post_segments[].gibbs_mlip.barrier_kcal`.
- Keep the backend and model, precision, force field, convergence criteria, restraints, thermochemistry settings, and temperature the same.
- Check each TS on its own: one imaginary mode along the reaction, its displacement, and both IRC ends before naming R and P.
- The same ML region and layers for both systems, the added or deleted atoms, and the charge and multiplicity of each system are in [Same atoms across states and variants](../mlmm-model-setup/SKILL.md#same-atoms-across-states-and-variants). Two radius-based selections can differ at the boundary, so compare the two `ml_region.pdb` files.

## See also

- [cli/tsopt.md](../mlmm-cli/tsopt.md), [cli/scan.md](../mlmm-cli/scan.md), [cli/path.md](../mlmm-cli/path.md), [cli/define-layer.md](../mlmm-cli/define-layer.md): running and judging each command.
- [cli/all-ts-only.md](../mlmm-cli/all-ts-only.md): TS-only mode.
- [outputs.md](outputs.md): R/TS/P paths and bond changes.
- [mlmm-hpc](../mlmm-hpc/SKILL.md): job templates and CPU or GPU resources.
- [`docs/backends.md`](../../docs/backends.md#determinism-and-reproducibility): fp64 and `--deterministic`.
