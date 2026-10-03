# `mlmm irc`

Integrates the intrinsic reaction coordinate (IRC) from a TS candidate in both
directions with EulerPC in mass-weighted Cartesians, and writes the path and
its two endpoint candidates. Run
`mlmm irc -i <ts> --parm7 <real.parm7> -q <charge> --out-json`.
The result is usable when both endpoint candidates, optimized with `opt`,
reach the intended R and P.

## When to use

- After `tsopt` gives n_imag = 1, to check which minima the TS candidate
  connects.
- The endpoints are raw IRC frames, not optimized minima. Run `mlmm opt` on
  them separately.

## Minimal run

```bash
mlmm irc -i result_tsopt/final_geometry.xyz --parm7 real.parm7 \
    --ref-pdb enzyme_layered.pdb -q 0 -m 1 -b uma --out-json -o result_irc
```

A smaller step and a longer trace for a shallow surface:

```bash
mlmm irc -i ts.xyz --parm7 real.parm7 --ref-pdb enzyme_layered.pdb -q -1 -m 1 \
    --max-cycles 250 --step-size 0.05 \
    -b uma -o result_irc_long
```

## Judge success

IRC has no independent scientific success verdict. A standalone `irc` does
not know which end is the reactant or the product. Judge it in three steps:

1. The start was a TS: the console line
   `Transition vector is mode 0 with wavenumber … cm⁻¹.` shows a negative
   wavenumber.
2. Each requested direction records its frame count (`n_frames_forward`,
   `n_frames_backward`) and `*_integration_stop_reason`.
   `*_integration_converged` describes whether the RMS-gradient stationarity
   criterion fired, so `--never-stop` leaves it false. This field and
   `*_downhill_departure_valid` are diagnostics, not endpoint-optimization
   gates. Finite retained endpoints can proceed to optimization after a
   predictor-budget or max-cycle stop. Missing or non-finite coordinates and
   execution errors must still be reported.
3. Optimize both endpoint candidates and compare them with the intended R and
   P. Even if the IRC does not converge, the result is usable when the
   optimized endpoints reach the intended R and P. The direction forward or
   backward does not decide which one is R.

```bash
mlmm opt -i result_irc/forward_first.xyz --ref-pdb enzyme_layered.pdb \
    --parm7 real.parm7 -q 0 -m 1 -o result_opt_forward
mlmm opt -i result_irc/backward_last.xyz --ref-pdb enzyme_layered.pdb \
    --parm7 real.parm7 -q 0 -m 1 -o result_opt_backward
```

Files in `result_irc/`:

```
finished_irc_trj.xyz      # whole path: first frame -> TS -> last frame
finished_first.xyz        # first frame of the whole path
finished_last.xyz         # last frame of the whole path
forward_irc_trj.xyz       # forward branch, from the TS
backward_irc_trj.xyz      # backward branch, from the TS
forward_first.xyz         # end of the forward branch (endpoint candidate)
backward_last.xyz         # end of the backward branch (endpoint candidate)
result.json               # with --out-json
```

The `.pdb` copies are written for PDB/mmCIF input or with `--ref-pdb`. With a
non-empty YAML `irc.prefix`, EulerPC inserts one underscore before each
filename (`prefix: trial` → `trial_finished_irc_trj.xyz`); read the normalized
names from `files` in `result.json`.

Read `energy_first_hartree` / `energy_last_hartree` and assign R and P after
inspecting or matching the endpoint structures. `never_stop` records whether
the opt-in mode was enabled; `never_stop_energy_bypasses` is the observed
bypass count.

```python
import json
d = json.load(open("result_irc/result.json"))
print(d["n_frames_forward"], d["n_frames_backward"])
print(d["energy_first_hartree"], d["energy_ts_hartree"], d["energy_last_hartree"])
print(d["execution_status"], d["scientific_status"])
print(d["forward_requested"], d["backward_requested"])
print(d["forward_integration_converged"], d["backward_integration_converged"])
print(d["forward_integration_stop_reason"], d["backward_integration_stop_reason"])
print(d["never_stop"], d["never_stop_energy_bypasses"])
```

`--read-hess` checks only the size, symmetry, and finiteness of the `.npy`
file, so pass a Hessian computed for the same geometry, charge,
multiplicity, layers, and calculator. Frozen atoms are treated as in
[freq.md](freq.md#phva).

## Bond changes

`bond_changes` records the directed difference from `finished_first.xyz` to
`finished_last.xyz` according to a 1.20× covalent-radius cutoff. This is the
same algorithm used by [bond-summary](utilities.md#bond-summary) and the
`path-search` segmentation. The direction is not a chemical R→P assignment;
for the R/TS/P conventions of `all`, see
[outputs.md](../mlmm-overview/outputs.md#oriented-rtsp-paths).

```python
import json
bc = json.load(open("result_irc/result.json")).get("bond_changes")
if bc is None:
    raise RuntimeError("IRC endpoint comparison was not available")
for b in bc["formed"]: print("FORMED ", b)
for b in bc["broken"]: print("BROKEN ", b)
```

## Pitfalls and recovery

- IRC starts from a TS with a single imaginary mode. If `tsopt` left several,
  IRC may follow the wrong one, so re-optimize the TS first.
- A branch that stops almost at once prints
  `[irc] IRC stopped after only a few frames in …`. Reduce `--step-size`
  first, for example to 0.05. Use `--never-stop` (or `all --irc-never-stop`)
  when tracing to the cycle cap is intended. It ignores the gradient and
  energy endpoint criteria; numerical or integration failures still stop the
  branch. It is off by default. Always inspect both branches and the bond
  connectivity.
- `--max-cycles 125` is enough for most systems. A branch that hits the cap
  still leaves a finite endpoint that goes on to endpoint `opt`; raise
  `--max-cycles` only when the branch must be followed further.
- The bond-change detector is geometry-based (covalent-radius cutoff), not
  physics-based. Metal–ligand bonds may flicker on the borderline.

## Next step

- Optimize the endpoints with [opt.md](opt.md); then [freq.md](freq.md) and
  [dft.md](dft.md).
- The IRC start comes from [tsopt.md](tsopt.md).
