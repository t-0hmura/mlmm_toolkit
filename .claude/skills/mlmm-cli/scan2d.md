# `mlmm scan2d`

## Purpose

Two-bond distance scan with harmonic restraints. Generates a 2D grid
of relaxed geometries along two simultaneously-varied bonds. Useful
for mapping a 2D potential surface (e.g., concerted vs. stepwise
distinction in a SN2 attack + leaving-group departure).

## Synopsis

```bash
mlmm scan2d -i input.pdb \
    -c 'A1 R1 N1' 'A2 R2 N2' --target1 1.6 \
    -c 'A3 R3 N3' 'A4 R4 N4' --target2 0.9 \
    [-b uma] [-o ./result_scan2d/]
```

## Key flags

| flag | type | description |
|---|---|---|
| `-c` (×2) | atom-spec pair | Each `-c` defines one bond to scan |
| `--target1` / `--target2` | float | Target distances (Å) for each bond |
| `--n-steps1` / `--n-steps2` | int | Grid resolution per axis |
| `-q` / `-l` / `-m` / `-b` / `-o` | — | Common conventions |

## Examples

```bash
mlmm scan2d -i 1.R.pdb \
    -c 'CS1 SAM 320' 'C7 GPP 321' --target1 1.6 \
    -c 'GPP 321 H11' 'GLU 186 OE2' --target2 0.9 \
    --n-steps1 10 --n-steps2 10 \
    -l 'SAM:1,GPP:-3' -b uma -o result_scan2d
```

## Output

Grid of `scan2d_NNxMM.xyz` plus a 2D energy map.
`result.json` stores `{ "grid": [[..., ...], ...] }` for downstream
plotting.

## See also

- `scan.md`, `scan3d.md` — 1D / 3D analogs.
- Defaults: `import mlmm.defaults as d; print(d.BIAS_KW)`
