# `mlmm trj2fig`

## Purpose

Plot an energy profile from an XYZ trajectory. Reads ASE-style `energy=...`
metadata from the comment line of each frame and produces a static PNG
or HTML plot. Useful for quickly visualizing IRC, MEP, or scan output.

## Synopsis

```bash
mlmm trj2fig -i trajectory.xyz [-o out.png] [--html]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | XYZ trajectory with energy in comment line |
| `-o, --output` | path | `out.png` | Output figure path |
| `--html` | flag | off | Write interactive Plotly HTML instead of static PNG |
| `--xlabel` / `--ylabel` / `--title` | str | sensible defaults | Plot labels |

## Examples

### Static PNG

```bash
mlmm trj2fig -i finished_irc_trj.xyz -o irc_profile.png
```

### Interactive HTML

```bash
mlmm trj2fig -i mep.xyz --html -o mep.html
```

## Caveats

- The XYZ comment line must encode the energy. ASE format `... energy=-1234.56`
  is recognized. Pure XYZ without energy fails silently with a flat plot.
- For a labeled energy diagram (R / TS / IM / P), use `energy-diagram.md`
  instead.

## See also

- `energy-diagram.md` — composed energy diagrams from explicit values.
- `irc.md`, `path-search.md`, `scan.md` — produce trajectories that
  feed `trj2fig`.
