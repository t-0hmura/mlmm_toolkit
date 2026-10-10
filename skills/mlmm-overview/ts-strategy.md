# mlmm ts-strategy

How to study a reaction mechanism with an ML/MM ONIOM model and get a correct barrier: precision, routes to a TS candidate, wrong n_imag, the barrier direction, multistep paths, splitting the reaction, a TS that does not come out, and controlled comparisons. Before any run, write down the bonds that form, the bonds that break, and every H atom that moves; this list becomes the `-s` coordinates and the points to check at the IRC endpoints, and every atom on it belongs in the ML region.

## 1. Precision: keep backend defaults

- Without `--precision`, UMA runs in fp32, ORB and MACE in fp64, and AIMNet2 in fp32 only, where an explicit fp64 is rejected. When precision matters, compare the supported precisions on the target system before production, and time one Hessian in each precision on the production GPU. fp64 is slow on consumer GPUs, where UMA in fp32 gives the best balance; on HPC GPUs, ORB in fp64 is cost-effective. ORB and MACE in fp32 leave extra imaginary modes more often.
- `--precision fp32|fp64` is case-insensitive and is accepted by `sp`, `opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `all`.
- fp64 changes numerical precision only. `--deterministic` requests deterministic algorithms, but exact repeatability still has to be checked on the installed backend, model, and hardware ([`docs/backends.md`](../../docs/backends.md#precision)). AIMNet2 rejects `--deterministic` as well as fp64.

## 2. Two routes to a TS candidate

- **MEP search**, when you have R and P, optionally with intermediates in order: `path-search` splits the path recursively with GSM or DMF and returns an HEI per segment; `path-opt` optimizes one given segment and returns its HEI.
- **Restrained build-up**, when there is no usable second endpoint or TS guess: `scan` (or `scan2d` / `scan3d`) drives the reacting coordinates with a harmonic restraint E = ½ k (r − target)², k = 300 by default (`--restraint-k`), and relaxes the rest at each step.
- Pick the scan atoms by measuring the structure, not from atom numbers in a paper: a forming bond is long in R and close to a bond length in P. Before scanning a proton transfer, measure the donor–acceptor heavy-atom distance in R; if it is well beyond a hydrogen-bond contact (roughly 3.2 Å) and no water or residue in the model bridges it, revise the hypothesis instead of driving it. If the TS or product puts the moving H on an atom you did not name, the coordinate was wrong; rebuild it rather than accepting the new acceptor.
- Scan relaxations and endpoint pre-optimization use `gau` unless you set `--thresh`, while TS and post-IRC endpoint optimizations in `all` use `--thresh-post baker`. When the scan seeds the TS candidate or its profile looks wrong (for example a barrier near zero), rerun it with `--thresh baker` and a smaller step (`--scan-max-step-size` in `all`, `--max-step-size` in `scan`) before changing the ML region.
- `opt` takes restraints only through `--distance-restraint`, with `--restraint-k` (also 300 by default). `scan` additionally drives staged targets up to a TS candidate.
- Pass a candidate from either route to `tsopt` and then `irc`, or use `all --tsopt`. The final PHVA (partial Hessian vibrational analysis) of `tsopt` gives the saddle order; add `freq` for the full modes or thermochemistry.

## 3. Wrong n_imag after TS optimization

A first-order saddle has exactly one imaginary mode, along the reaction; check its displacement and the IRC ends. Two or more fail, however small the extra imaginary frequencies are.

- **A small extra mode in a model with a Frozen-MM layer** (a few to about 20 cm⁻¹): first re-optimize with a tighter preset (`tsopt --thresh gau_tight` or `gau_vtight`, or `all --thresh-post`) and count n_imag again. If the small mode vanishes or changes sign, it was residual curvature; only a persistent one calls for `--flatten`.
- **Extra imaginary modes**: look at every mode's displacement, the MEP guess, the optimizer's stop reason, and the backend's numerical behavior. Then retry with another coordinate, flattening, or precision setting and read the new final PHVA.
- **n_imag = 0**: a failure, not a TS. Improve the MEP or the starting structure; `--flatten` only removes surplus modes and cannot create a missing reaction direction.
- **Poor MEP or HEI**: in `all`, try `--refine-path` before the TS optimization. It can split a poor path into several stages and costs more, so it is off by default.
- **Still no clean saddle**: revise the endpoints or the scan coordinates, and check that the single imaginary mode moves the reacting atoms.

Settings to retry with:

- `--flatten` runs the loop that removes extra imaginary modes, on `opt`, `tsopt`, and `all`, and is off by default (`--no-flatten`); with it, the TS commands run the loop when n_imag > 1. `--flatten` uses the iteration cap (50 by default) and `--no-flatten` sets it to zero.
- If n_imag is still ≥ 2 after `--flatten`, read `flatten_skip_reason` in the tsopt `result.json`. `target mode is not negative` (or `target mode sign never determined`) means that the reaction mode followed along a reference direction (the MEP tangent in `all`, or `--ref-mode`) is no longer imaginary, so flattening was refused to avoid making a saddle of another reaction: start again from the HEI or another candidate instead of retrying `--flatten`. `max-cycles budget exhausted before flattening` (or `during flattening`) means the flatten rounds share `--max-cycles`.
- `--coord-type` takes `cart`, `redund`, `dlc`, or `tric` in `opt` and `tsopt`, and `cart` or `dlc` in `all`; the default is `cart`. Keep `cart` for ML/MM: `dlc` (delocalized internal coordinates) is slow to build for these models and its convergence depends on the system, so compare it with `cart` from the same start. `opt` accepts `dlc` with L-BFGS (`--opt-mode grad`) or RFO (`--opt-mode hess`). `path-opt` and `path-search` have no `--coord-type`; they read `geom.coord_type` from the YAML, `cart` or `dlc` only.
- Check any change of coordinate system with a frequency analysis and the IRC ends.
- `--ref-mode` is an advanced input, not a routine fix. It takes Cartesian 3N vectors in the input atom order (`.npz`, `.npy`, or text) and guides which negative Hessian root the TS optimizer follows; it does not replace the Hessian, and Dimer ignores it. `all` supplies the MEP tangent this way by default; with `all --no-tsopt-from-mep-tan`, `tsopt` picks its starting root from the Hessian of the starting structure.

Keep convergence and saddle order apart. A converged higher-order saddle is not a TS, although `all` still runs a diagnostic IRC, with a warning, when a valid negative root exists. No convergence, no imaginary mode, a failed or skipped final Hessian, or no valid negative root stops `all` before the IRC, with the TS files kept. The full list of remedies is in [Wrong imaginary-mode count after optimization](../../docs/tsopt.md#wrong-imaginary-mode-count-after-optimization).

## 4. Reading the barrier when the scan started from P

If the scan or path starts from P, the raw barrier is the reverse one.

- Forward barrier: E(TS) − E(R). Reverse barrier, the raw number of a P-start run: E(TS) − E(P). In `summary.json`, the barrier from the other end is `barrier_kcal − delta_kcal`.
- This is how you read the numbers, not a CLI flag. Before quoting a barrier, compare `segments/seg_NN/reactant.*` and `product.*` with the intended R and P instead of trusting the scan direction. In TS-only mode, R is the higher-energy IRC end ([outputs.md](outputs.md#oriented-rtsp-paths)).
- The top of a scan is only a TS candidate. Read the barrier from the energy after the TS optimization (`post_segments[].mlip.barrier_kcal`), counted from the minimum just before that stage.

## 5. Multistep paths

For each TS of a multistep path, report the local barrier (counted as above) and the TS height above one common reference: the lowest optimized R of the same parm7, ML region, and protonation state, which after IRC and endpoint optimization is often several kcal/mol below the structure you started from. Do not add local barriers, do not count a drop that comes only from switching the reference R, and do not join steps from different topologies, ML regions, snapshots, or protonation states into one profile.

When only one step gives a validated TS, or the optimized IRC ends show only part of the intended change (for example the heavy-atom transfer without the proton transfer), treat that optimized end as an intermediate, not a failure. Start the next step from it rather than from a scan stage end: a new Scan-list run from `segments/seg_NN/product.*` (or `reactant.*`) without `-c`, with the same `--parm7` and `--model-pdb` (the `mm_parm/` parm7 and `ml_region.pdb` of the earlier run), charge, and multiplicity, driving only the remaining coordinates. Before joining steps computed separately, compare the optimized P of one step with the optimized R of the next: covalent bonds, the owner of every H, each residue's protonation, then all-atom RMSD and energy. Identical bonding can still differ by a few tenths of an Å and several kcal/mol; connect that conformational gap with a path between the two minima, or report it as unverified. Then draw all steps on one diagram from the common R with `energy-diagram` ([outputs.md](outputs.md#energy-diagrams)).

## 6. Staged vs concerted scan

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

- When the TS does not come out, try putting every moving coordinate in its stage: the bond that breaks and each H that moves, not only the bond that forms.
- **Hold what should not move yet**: a stage restrains only the coordinates it lists, so an X–H that should react later can move in an earlier stage. Add that distance to each earlier stage with its target set to the value measured on the structure the stage starts from; a tuple whose target equals its start does not move, so its restraint holds the distance until the next stage starts. For the first stage, `all` pre-optimizes the input inside the scan, so measure on that result (`_work/scan/preopt/result.*`, or your own `opt` result) and start `all --no-scan-preopt` from it. The hold is harmonic: check the held distance in `scan_trj.xyz` and again at the MEP ends.
- **A stage end is not an intermediate**: it is a restrained structure, and the MEP pre-optimizes each stage end without restraints (`--preopt`, on by default), which often breaks a bond the stage formed or re-forms one it broke. Compare the driven distances in `_work/scan/stage_NN/result.*` with the same distances at the MEP ends, or rerun with `--scan-endopt` (off by default) and check whether the intermediate survives once the restraints are released; if it relaxes back, the concerted path fits better.
- **Count the steps from the TSs**, not from the `-s` stages: one literal with several tuples often refines into separate TSs, and two stages can collapse into one segment. For each TS, check which bonds and H atoms its mode moves and where they cross in the IRC frames. A proton still on its donor at the TS and at the raw IRC end (`segments/seg_NN/structures/*_irc.*`) that moves only during endpoint optimization completes downhill after that TS; it shows neither a concerted event nor a separate step.
- In `scan`, a 4-tuple `(i,j,low,high)` scans one distance both ways from the input and becomes two stages; in `all`, a 4-tuple is an angle target.

## 7. When the TS does not come out

A TS needs n_imag = 1 and IRC ends that are the intended R and P; an exit code of 0 is not that evidence.

- **Extra modes**: `--flatten` and `--refine-path` are in [Wrong n_imag after TS optimization](#3-wrong-n_imag-after-ts-optimization). After `--flatten`, check both IRC ends again, since n_imag can reach 1 on a TS of another reaction.
- **An extra mode that moves the ML region together with movable MM atoms**, such as a water moving with the substrate: with microiteration, the TS steps move only the ML atoms and the MM atoms bonded to them, and the other movable MM atoms only relax. Rerun `tsopt --no-microiter --flatten` so that RS-P-RFO steps over every movable atom. `all` has no microiteration switch, so pass the `--parm7` and `--model-pdb` of the `all` run to `tsopt`.
- **A narrowed Hessian**: if `tsopt` ran with a small `--hessian-cutoff`, widen it or leave it out.
- **Start from the extra mode**: `vib/imag_*_trj.xyz` and `vib/imag_*.pdb` hold 20 frames, and frames 6 and 16 are the largest displacements either way. Save each as a PDB and start TS-only mode or `tsopt` from each, with the same charge and multiplicity.
- **No imaginary mode, the candidate slid toward R or P, or the step's MEP barrier is low**: the HEI may not yet be a TS. Check its imaginary mode and the forming and breaking bond lengths, then read the MEP energy profile. Two maxima mean the pair spans two steps: give the intermediate or use `--refine-path` rather than raising `--max-nodes`. For a single HEI that does not look like a TS, rerun with `--refine-path` or change the path settings, which all move the HEI: the number of internal images (`--max-nodes`, 20 by default), `--mep-mode dmf`, and `--gsm-param energy`; or run a fine scan between the two ends of that step and take a TS-like frame. Other starting structures are the HEI of the unsplit MEP (`hei.xyz` of `path-opt`), the HEI of each segment (`hei_seg_NN.xyz` of `path-search`; under `_work/path_opt/` in `all`, or `_work/path_search/` with `--refine-path`), or a scan frame near the top. If the forming and breaking bonds are both long at the HEI, or the HEI lies after the bond broke, the TS optimization tends to lose the reaction mode. Once a TS has n_imag = 1 and its IRC ends are the intended R and P, the node count of the MEP that seeded it does not change the result; in a batch of snapshots, allow one such path retry per snapshot, then move to the next snapshot.
- **An extra mode that also moves the reacting bonds**: two stages may overlap in one candidate; run them as separate stages.
- **Split the reaction differently**: run each plausible stage order and compare the energy diagrams; to choose between stepwise and concerted, see [Staged vs concerted scan](#6-staged-vs-concerted-scan). When R and P can be prepared, give both to an MEP search instead of a scan.
- **ML region**: when a residue, water, or cofactor that takes part lies outside the ML region, enlarge it with a larger `-r` or a hand-built `--model-pdb` ([Enlarge when the model is too small](../mlmm-model-setup/SKILL.md#enlarge-when-the-model-is-too-small)), and recheck the charge.

The full guide is [`docs/mechanism-tips.md`](../../docs/mechanism-tips.md).

## 8. Controlled comparisons

- Form each barrier within one system (TS − R or TS − P), then compare the barriers: ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT. Do not subtract total energies of systems with different compositions.
- Give R and P of each system in Endpoint mode, so that R is the chemical reactant; G_TS − G_R is `post_segments[].gibbs_mlip.barrier_kcal`.
- Keep the backend and model, precision, force field, convergence criteria, restraints, thermochemistry settings, and temperature the same.
- Check each TS on its own: one imaginary mode along the reaction, its displacement, and both IRC ends before naming R and P.
- The same ML region and layers for both systems, the added or deleted atoms, and the charge and multiplicity of each system are in [Same atoms across states and variants](../mlmm-model-setup/SKILL.md#same-atoms-across-states-and-variants). Two radius-based selections can differ at the boundary, so compare the two `ml_region.pdb` files.
- **Mechanism candidates**: use one parm7 and one ML region (`--parm7`, `--model-pdb`) for every hypothesis, so that all barriers share the same R. Build each candidate explicitly (its own `-s` stages or endpoints, concerted and stepwise); an MEP or TS search that ends on the saddle of one mechanism, or a model with a reacting fragment removed, is not evidence against another. Rank the candidates by barriers whose TS has n_imag = 1 and the intended IRC ends, after `--tsopt --thermo` on every candidate that could be the lowest, not by MEP-level barriers (`segments[].barrier_kcal`, or `rate_limiting_step` with `method` `MEP`): the verified barrier can differ by many kcal/mol and reverse the order.
- **MD snapshots**: one structure gives one static barrier. Run the same protocol and flags on several snapshots that sample the reactive arrangement ([building their models](../mlmm-model-setup/SKILL.md#build-the-ml-region)), and compare barriers, not total energies. Before counting, deduplicate TSs by energy and geometry: several scan orders or settings from one snapshot often reach the same TS, and two inputs built from the same endpoints are one candidate. Report the TSs, the snapshots with a validated TS, and the fully connected paths as separate counts, with the lowest and median barrier rather than a single value.
- **Repeat runs**: runs are not bitwise repeatable by default. On the same input and GPU model, optimizer paths can separate within a few tens of cycles and end at different structures, or one TS run can converge while the other stops. Do not attribute a difference between two single runs to the setting you changed, and do not discard a candidate after one failed run; repeat the run, or compare with `--deterministic` on the same hardware and software.

## See also

- [cli/tsopt.md](../mlmm-cli/tsopt.md), [cli/scan.md](../mlmm-cli/scan.md), [cli/path.md](../mlmm-cli/path.md), [cli/define-layer.md](../mlmm-cli/define-layer.md): running and judging each command.
- [cli/all-ts-only.md](../mlmm-cli/all-ts-only.md): TS-only mode.
- [outputs.md](outputs.md): R/TS/P paths and bond changes.
- [mlmm-hpc](../mlmm-hpc/SKILL.md): job templates and CPU or GPU resources.
- [`docs/backends.md`](../../docs/backends.md#determinism-and-reproducibility): fp64 and `--deterministic`.
