# Quickstart: `mlmm tsopt` -> `mlmm freq`

## Goal

Optimize a TS candidate and validate it by frequency analysis.

## Prerequisites

- TS candidate geometry: `ts_guess.pdb`
- MM topology: `real.parm7`
- ML region definition: `ml_region.pdb`

## 1. TS optimization

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

## 2. Frequency check on optimized TS

```bash
mlmm freq -i ./result_tsopt/final_geometry.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## What to check

- `result_tsopt/final_geometry.pdb`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz` and `result_freq/mode_*.pdb`

For a valid first-order saddle, frequencies should contain exactly one imaginary mode (negative cm^-1).

## Tips

- Use `--hessian-calc-mode Analytical` when VRAM is sufficient.
- Check full options with `mlmm tsopt --help-advanced` and `mlmm freq --help-advanced`.

## Next step

- Trace the path with [irc](irc.md), or run the end-to-end route with [all](all.md).
