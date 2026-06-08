# `mlmm trj2fig`

## Purpose

Plot an energy profile from an XYZ trajectory. Reads the first
floating-point number found in each frame's comment line (interpreted as
the Hartree energy) and exports the figure or CSV. Output format is
inferred from the filename suffix (`.png` / `.html` / `.svg` / `.pdf` /
`.csv`). Useful for quickly visualizing IRC, MEP, or scan output.

## Synopsis

```bash
mlmm trj2fig -i trajectory.xyz [-o out.png] [-o out.html] [--unit kcal|hartree]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | XYZ trajectory with energy in comment line |
| `-o, --out` | path (multi) | `energy.png` | Output file(s); repeat `-o` or list extras positionally. Format from suffix. |
| `--unit` | choice | `kcal` | Energy unit (`kcal` or `hartree`). |
| `-r, --reference` | str | `init` | Reference: `init` (initial frame), `None` (absolute E), or integer index. |
| `-q, --charge` | int | — | Total charge; recompute energies when supplied. |
| `-m, --multiplicity` | int | — | Spin multiplicity (2S+1); recompute energies when supplied. |
| `--reverse-x` / `--no-reverse-x` | flag | off | Reverse the x-axis (last frame on the left). |

## Examples

### Static PNG

```bash
mlmm trj2fig -i finished_irc_trj.xyz -o irc_profile.png
```

### Interactive HTML

```bash
mlmm trj2fig -i scan_trj.xyz -o mep.html
```

### Multiple outputs in one call

```bash
mlmm trj2fig -i scan_trj.xyz -o profile.png -o profile.csv
```

## Caveats

- The XYZ comment line must contain the energy as a decimal or scientific
  float (e.g. `-1234.56` or `-1.23e3`); the first such number is used. A bare
  integer is rejected to avoid mistaking a frame index for an energy.
- For a labeled energy diagram (R / TS / IM / P), use `energy-diagram.md`
  instead.

## See also

- `energy-diagram.md` — composed energy diagrams from explicit values.
- `irc.md`, `path-search.md`, `scan.md` — produce trajectories that
  feed `trj2fig`.
