# `mlmm trj2fig`

## Purpose

Plot an energy profile from an XYZ trajectory. Reads a strict energy
record from each frame's comment line (interpreted as Hartree) and exports
the figure or CSV. Accepted records are a lone decimal/scientific value,
an explicit `E=<value>` or `Energy=<value>` field, or the pysisyphus
`<value> , ...` form. Output format is
inferred from the filename suffix (`.png` / `.jpg` / `.jpeg` / `.html` /
`.svg` / `.pdf` / `.csv`). Useful for quickly visualizing IRC, MEP, or scan
output.

## Synopsis

```bash
mlmm trj2fig -i trajectory.xyz [-o out.png] [-o out.html] [--unit kcal|hartree] [--out-json]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | XYZ trajectory with energy in comment line |
| `-o, --out` | path (multi) | `energy.png` | Output file(s); repeat `-o`, or give additional output paths as positional arguments. Format from suffix. |
| `--unit` | choice | `kcal` | Energy unit (`kcal` or `hartree`). |
| `-r, --reference` | str | `init` | Reference: `init` (initial frame), `None` (absolute E), or integer index. |
| `-q, --charge` | int | — | Total charge; recompute energies when supplied. |
| `-m, --multiplicity` | int | — | Spin multiplicity (2S+1); recompute energies when supplied. |
| `-b, --backend` | choice | `uma` | Recomputation backend: `uma`, `orb`, `mace`, or `aimnet2`. |
| `--backend-model` | str | backend default | Model variant used for recomputation. |
| `--precision` | choice | backend default | Case-insensitive `fp32` or `fp64`. |
| `--out-json` / `--no-out-json` | flag | off | Write `result.json` beside the first output. |
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

### Recompute and record provenance

```bash
mlmm trj2fig -i trajectory.xyz -q 0 -m 1 -b uma \
    --backend-model uma-s-1p2 --precision fp32 -o profile.png --out-json
```

## Caveats

- Use a lone decimal/scientific energy (for example `-1234.56` or
  `-1.23e3`), an explicit `E=<value>` / `Energy=<value>` field, or the
  pysisyphus `<value> , ...` form. Other surrounding text is rejected as
  ambiguous, and a bare integer is rejected as a possible frame index.
- Either `-q` or `-m` recomputes every frame with the selected MLIP; omitted values resolve to charge 0 and multiplicity 1. This is a direct MLIP rescore, not an ONIOM energy.
- Comment-mode JSON records null backend/model/precision/charge/multiplicity; recomputation records the resolved provenance.
- For a labeled energy diagram (R / TS / IM / P), use `energy-diagram.md`
  instead.

## See also

- `energy-diagram.md` — composed energy diagrams from explicit values.
- `irc.md`, `path-search.md`, `scan.md` — produce trajectories that
  feed `trj2fig`.
