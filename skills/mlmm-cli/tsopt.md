# `mlmm tsopt`

Optimizes a TS candidate on the ML/MM model, then computes the Hessian at the
final geometry and counts its imaginary modes (n_imag). Run
`mlmm tsopt -i <candidate> --parm7 <real.parm7> -q <charge> --out-json`.
A successful TS optimization gives one imaginary mode along the reaction
coordinate: `optimization_status` is `converged` and `saddle_validation` is
`first_order`.

## When to use

- Refine a HEI candidate from `path-opt` / `path-search` ([path.md](path.md))
  or from a [scan](scan.md) into a TS candidate.
- Check a TS candidate built elsewhere.

The default optimizer is RS-P-RFO (`--opt-mode hess`); RS-I-RFO, TRIM, and
the Hessian-guided Dimer are alternatives.

## Minimal run

```bash
mlmm tsopt -i hei.xyz --parm7 real.parm7 --ref-pdb enzyme_layered.pdb \
    -q 0 -m 1 -b uma --out-json -o result_tsopt
```

An XYZ candidate needs `--ref-pdb`; the other ML/MM flags are in
[SKILL.md](SKILL.md#shared-mlmm-conventions).

Dimer:

```bash
mlmm tsopt -i hei.xyz --parm7 real.parm7 --ref-pdb enzyme_layered.pdb -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

RS-I-RFO on MACE, with the ML-region charge from residue charges:

```bash
mlmm tsopt -i hei.xyz --parm7 real.parm7 --ref-pdb enzyme_layered.pdb \
    -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Judge success

How the run ended decides whether n_imag exists:

- `converged`: the console prints `[microiter] Converged!` or
  `[tsopt] Numerical optimization converged.`, then
  `[tsopt] Wrote N final imaginary mode(s).`. The final Hessian is computed.
- `stalled`: the energy-plateau stop of `--stop-plateau` (off by default)
  fired, and the console prints `Stalled (energy plateau; not converged)`. The
  Hessian is still computed, so n_imag is reported for an unconverged geometry.
- `not_converged`: `--max-cycles` was reached
  (`[microiter] Reached max macro iterations (M).` or
  `[tsopt] Reached max cycles (N/M).`). No Hessian is computed, so
  `hessian_status` is `skipped`, `n_imaginary_modes` is null, and the saddle
  order is unknown.

`saddle_validation` is `first_order` for n_imag = 1, `higher_order` for 2 or
more, `no_imaginary` for 0, and `unavailable` without a Hessian. A mode counts
as imaginary when ν < −5.00 cm⁻¹. Standalone `tsopt` judges convergence only:
`converged` gives `scientific_status` `success` and exit 0 whatever n_imag is,
and `stalled` or `not_converged` gives `failed` and exit 1. A failed final
Hessian (`[tsopt] ERROR: Terminal PHVA failed.`) sets `hessian_status` to
`failed` and exits 1. Read n_imag yourself:

```python
import json
d = json.load(open("result_tsopt/result.json"))
status = d["optimization_status"]        # converged / stalled / not_converged
n = d["n_imaginary_modes"]
if d["hessian_status"] != "completed":
    print(status, "no n_imag:", d["hessian_status"], d["hessian_error"])
elif n == 1:
    print(status, "single imaginary mode at", d["imaginary_frequencies_cm"][0], "cm-1")
elif n == 0:
    print(status, "no imaginary mode: collapsed to a minimum")
else:
    print(status, "multiple imaginary modes; inspect vib/imag_*")
print(d["energy_hartree"], d["files"]["final_geometry_xyz"])
```

A converged run with n_imag = 1 is still a candidate until [IRC](irc.md)
shows that it connects the expected R and P.

## Choosing --opt-mode

- `hess` / `rsprfo` (default): RS-P-RFO, a partitioned restricted-step
  treatment. Memory and runtime depend on the active DOFs and the backend.
- `rsirfo`: RS-I-RFO, an image-function alternative on the same
  microiteration driver. `trim` selects TRIM.
- `grad` / `dimer`: Hessian-guided Dimer. Its initial and periodic orientation
  Hessians make it more robust than a random initial direction for large
  systems; convergence still depends on the seed and the system.

If Dimer stalls, inspect the followed mode and the step diagnostics, then
compare RS-P-RFO or RS-I-RFO on the same seed rather than using a universal
cycle threshold.

## Pitfalls and recovery

- `tsopt` always forces `reject_uphill=False`, regardless of optimizer mode or
  YAML. Uphill trial steps can be part of saddle-point mode following. The
  `--reject-uphill/--no-reject-uphill` toggle belongs only to `opt` and to the
  endpoint optimization after IRC in `all`.
- n_imag = 0 is a failed TS optimization, even when the optimizer stopped
  normally. `--flatten` can remove surplus imaginary modes but cannot create a
  missing reaction direction. Improve the MEP or the starting guess;
  `all --refine-path` is opt-in because recursive refinement can split a poor
  path into several costly segments.
- n_imag of 2 or more: watch each `vib/imag_*_trj.xyz` to decide whether the
  extra modes are spurious or a real higher-order saddle point, and re-optimize
  with `--flatten`; see
  [ts-strategy.md](../mlmm-overview/ts-strategy.md#3-wrong-n_imag-after-ts-optimization).
- `--ref-mode` is an advanced input; `all` supplies the MEP tangent this way
  by default; with `all --no-tsopt-from-mep-tan`, the root comes from
  the Hessian modes of the starting structure. Standalone runs omit it. Supply
  it by hand only when the non-zero 3N vector uses exactly the same atom order
  as the TS input.
- `--max-cycles` is a safety cap, not evidence of correctness. On repeated
  non-convergence, inspect the TS seed, the followed mode, the optimizer
  diagnostics, and the backend and model behavior; see
  [ts-strategy.md](../mlmm-overview/ts-strategy.md#6-when-the-ts-does-not-come-out).
- The backend and model change the curvature surface. Validate every candidate
  by exactly one imaginary mode, its displacement, and the intended IRC
  connectivity.

## Outputs

```
result_tsopt/
├── final_geometry.{xyz,pdb}    # final geometry; check result.json status
├── vib/imag_*_trj.xyz, .pdb    # animation of each imaginary mode
├── optimization_all_trj.xyz    # with --dump
└── result.json                 # with --out-json
```

## Next step

- n_imag = 1: run [irc.md](irc.md); [freq.md](freq.md) for thermochemistry.
- A new candidate: [path.md](path.md) or [scan.md](scan.md).
- Backends for the TS step:
  [UMA](../mlmm-install/backends.md#uma),
  [MACE](../mlmm-install/backends.md#mace-separate-environment).
