# `mlmm energy-diagram`

## Purpose

Plot an energy diagram from numeric inputs. Use this when you have
energy values from one or more `mlmm-toolkit` runs (e.g. R from one
run, TS / IM from another) and want a single composite figure.

## Synopsis

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" \
    [-o energy_diagram.png]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | numeric sequence | required | Energy values. Accepts a Python-list literal (`-i "[0, 12.5, 4.3]"`), repeated `-i` calls (`-i 0 -i 12.5 -i 4.3`), or a bare space-separated list after one flag (`-i 0 12.5 4.3`) — all three forms are equivalent |
| `-o, --output` | path | `energy_diagram.png` | Output image path |
| `--help-advanced` | flag | — | Reveal x-axis labels (`--label-x`) and y-axis label (`--label-y`) |

State labels and y-axis units are exposed via `--help-advanced`. Without
labels, points are plotted in input order.

## Examples

### Five points (two-step mechanism)

```bash
mlmm energy-diagram -i 0.0 -i 21.5 -i -0.7 -i 2.2 -i -18.2 -o diagram.png
```

### Bracketed list literal

```bash
mlmm energy-diagram -i "[0.0, 21.5, -0.7, 2.2, -18.2]" -o diagram.png
```

## Caveats

- Energies are taken **as-is**; the command does no unit conversion.
  Make sure all values share a unit (kcal/mol is typical) before
  calling.
- For a profile along a continuous trajectory (XYZ frames with energies
  in the comment line), use `trj2fig.md`.
- Per-state x-axis labels (`--label-x`) and the y-axis label (`--label-y`)
  live behind `--help-advanced`; refer to that for production figures.

## See also

- `trj2fig.md` — plot from a trajectory file.
- `../mlmm-workflows-output/SKILL.md` — extracting per-segment
  energies from `summary.json` to feed this command.