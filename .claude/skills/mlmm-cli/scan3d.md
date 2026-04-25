# `mlmm scan3d`

## Purpose

Three-bond distance scan. Output is a 3D grid of relaxed geometries.
Rare in practice; usually you can decompose to multiple 2D scans or
sequential 1D scans (see `all-scan-list.md`). Provided for completeness
when an inherently 3D coupled mechanism is mapped (e.g. proton-coupled
electron transfer with two protons + one electron-donor distance).

## Synopsis

```bash
mlmm scan3d -i input.pdb \
    -c <bond1> --target1 ... \
    -c <bond2> --target2 ... \
    -c <bond3> --target3 ... \
    [-b uma] [-o ./result_scan3d/]
```

## Caveats

- Output volume scales as `n_steps1 × n_steps2 × n_steps3` — keep grid
  sizes modest (5×5×5 = 125 points already).
- Wall-time grows quickly; consider GPU multi-worker (`-b uma` with
  `workers > 1`) for large grids.

## See also

- `scan.md`, `scan2d.md` — 1D / 2D analogs.
- `all-scan-list.md` — staged scans without the grid blow-up.
