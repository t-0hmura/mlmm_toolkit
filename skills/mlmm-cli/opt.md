# `mlmm opt`

Relaxes one structure on the ML/MM model to its nearest local minimum with
L-BFGS (default) or RFO. Run
`mlmm opt -i <structure> --parm7 <real.parm7> -q <charge> --out-json`.
Success is the console line `[opt] Converged!` and `optimization_status`
`converged`; then confirm n_imag = 0 with `freq`.

## When to use

- Relax a starting structure before `path-opt` / `path-search`.
- Refine the IRC endpoints into R and P.
- Not a TS optimizer; for a TS, use [tsopt.md](tsopt.md).

## Minimal run

```bash
mlmm opt -i my.pdb --parm7 real.parm7 -l 'SAM:1' -b uma --out-json -o result_opt
```

RFO for stiffer convergence:

```bash
mlmm opt -i my.xyz --parm7 real.parm7 --ref-pdb topology.pdb -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

Pre-relax endpoints before `path-opt`:

```bash
mlmm opt -i 1.R.pdb --parm7 real.parm7 -l '...' -o /tmp/relax_R
mlmm opt -i 3.P.pdb --parm7 real.parm7 -l '...' -o /tmp/relax_P
mlmm path-opt -i /tmp/relax_R/final_geometry.xyz /tmp/relax_P/final_geometry.xyz \
    --parm7 real.parm7 --ref-pdb 1.R.pdb \
    -l 'SAM:1,GPP:-3' -o result_path_opt
```

## Judge success

- `converged`: `[opt] Converged!`, or `[microiter] Converged!` for RFO with
  microiteration. `scientific_status` is `success` and the exit code 0.
- `not_converged`: `[opt] Reached max cycles (N/M).` `failed`, exit 1.
- `stalled`: only with `--stop-plateau` (off by default),
  `[opt] Stalled (energy plateau; not converged)`. `failed`, exit 1.

With `--out-json`, `result.json` holds `optimization_status`, `n_opt_cycles`,
`energy_hartree`, `final_max_force`, and `files.final_geometry_xyz`. The
outputs are `final_geometry.{xyz,pdb}`, plus `optimization_trj.xyz` and
`optimization_all_trj.xyz` with `--dump`.

With `--thresh baker`, convergence requires ALL of `max(|force|) <= 3e-4`,
`rms(force) <= 2e-4`, `max(|step|) <= 3e-4`, `rms(step) <= 2e-4` and
`|delta E| < 1e-6`. This is a deliberately tightened variant of the published
criterion.

## Choosing --opt-mode

- `grad` / `lbfgs` (default): L-BFGS. Fast and robust for most
  well-conditioned minima.
- `hess` / `rfo`: RFO with Hessian updates. Stiffer convergence; useful when
  L-BFGS oscillates.

## Pitfalls and recovery

- `--reject-uphill/--no-reject-uphill` (off by default) opts in to rejecting
  an energy-raising RFO trial above `1e-4` Hartree, restoring the lower-energy
  geometry and shrinking the trust radius. At the smallest trust radius, `opt`
  runs one final convergence check on the retained geometry. It is ignored in
  L-BFGS mode.
- L-BFGS occasionally walks past a saddle on shallow surfaces. If `freq` shows
  imaginary frequencies, re-run with `--opt-mode rfo`, or add `--flatten`,
  which removes extra imaginary modes after the optimization. With
  `--flatten`, frozen atoms are treated as in [freq.md](freq.md#phva).
- `--mm-only` skips the MLIP component and minimizes on the MM force field
  only, as a cheap pre-relaxation before the ML/MM optimization. It supports
  only `--opt-mode grad` and turns microiteration off.
- Less common settings such as step limits and the trust radius go in a
  `--config` YAML, under its `opt`, `lbfgs`, and `rfo` sections.

## Next step

- [freq.md](freq.md): verify the optimized minimum (n_imag = 0).
- [path.md](path.md): MEP between relaxed endpoints.
- [tsopt.md](tsopt.md): the TS counterpart.
