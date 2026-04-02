# Quickstart: `mlmm tsopt`

## Goal

Optimize a TS candidate and verify that it is a first-order saddle point.

## Prerequisites

- TS candidate geometry: `ts_guess.pdb`
- MM topology: `real.parm7`
- ML region definition: `ml_region.pdb`

## 1. TS optimization

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

`tsopt` performs a final Hessian evaluation and imaginary-frequency check automatically at the end of optimization. Check the terminal output for lines like:

```
[Imaginary modes] n=1  ([-593.1])
```

## What to check

- `result_tsopt/final_geometry.pdb` — optimized TS structure
- `result_tsopt/vib/` — animation files for the imaginary-frequency normal mode (`imag_*.xyz`, `.pdb`)
- Terminal output: **n=1** with a sufficiently large imaginary frequency (|ν| ≥ 100 cm⁻¹) indicates a good TS candidate

## 2. (Optional) Separate frequency analysis

A standalone `freq` run is useful when you want full vibrational frequency output or thermochemistry corrections (`--thermo` in the `all` command). If you only need the imaginary-frequency check, the `tsopt` output above is sufficient.

```bash
mlmm freq -i ./result_tsopt/final_geometry.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## Tips

- Use `--hessian-calc-mode Analytical` when VRAM is sufficient.
- To use a different MLIP backend, add `-b orb` (or `mace`, `aimnet2`). Default is `uma`.
- Add `--embedcharge` to enable xTB point-charge embedding for MM-to-ML environmental corrections.
- Check full options with `mlmm tsopt --help-advanced` and `mlmm freq --help-advanced`.

## Next step

- Trace the path with [irc](irc.md), or run the end-to-end route with [all](all.md).
