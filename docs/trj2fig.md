# `trj2fig`

`mlmm trj2fig` reads the Hartree energies encoded in each frame's comment line of an XYZ trajectory, converts them to kcal/mol or Hartree, optionally references all values to a chosen frame, and exports the resulting series as static/interactive figures and CSV tables. Supplying `-q/--charge` or `-m/--multiplicity` instead recomputes every frame with the selected MLIP backend. This recomputation is a direct MLIP frame rescore, not an ML/MM ONIOM energy. The figure uses bold ticks, consistent fonts, markers, and a smoothed spline curve (no title).

## Examples

Default PNG, relative energy with respect to the first frame:

```bash
# Default PNG, relative energy with respect to the first frame
mlmm trj2fig -i traj.xyz
```

CSV + SVG with reference frame #5, reported in Hartree:

```bash
# CSV + SVG with reference frame #5, reported in Hartree
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree
```

Multiple outputs in one run with x-axis reversed:

```bash
# Multiple outputs in one run with x-axis reversed
mlmm trj2fig -i traj.xyz -o energy.png energy.html energy.pdf --reverse-x
```

Recompute every frame with an explicit backend configuration and write JSON provenance:

```bash
mlmm trj2fig -i traj.xyz -q 0 -m 1 -b uma --backend-model uma-s-1p2 \
    --precision fp32 -o energy.png energy.csv --out-json
```

## Workflow

1. Parse the XYZ trajectory. With neither `-q/--charge` nor `-m/--multiplicity`, extract Hartree energies from each frame's comment line. Supplying either option recomputes every frame with the selected backend; omitted recomputation values resolve to charge 0 and multiplicity 1.
2. Normalize the reference specification:
    - `init` -- frame `0` (or the last frame when `--reverse-x` is active).
    - `None`/`none`/`null` -- absolute energies (no referencing).
    - Integer literal -- the corresponding 0-based frame index.
3. Convert energies to either kcal/mol (default) or Hartree and, when a
    reference is active, subtract the reference value to produce delta-E.
4. Build the Plotly figure (bold ticks, spline interpolation, markers, no
    title) and export it to every requested extension.
5. Optionally emit a CSV table of the per-frame energies (see Outputs for the column layout).

## Outputs

```
<output>.[png|jpg|jpeg|html|svg|pdf] # Plotly export for every requested extension (defaults to energy.png)
<output>.csv # Optional energy table when CSV is requested
result.json # Machine-readable result and energy provenance with --out-json
summary.json # Identical machine-readable mirror with --out-json
```
- When no `-o` or positional outputs are provided, a single `energy.png` is written
  to the current directory.
- CSV exports include `frame`, `energy_hartree`, and either a delta-E column
  (`delta_kcal`/`delta_hartree`) or an absolute column (`energy_kcal`/`energy_hartree`
  when no reference is applied).
- PNG uses Plotly's PNG export with `scale=2` for higher resolution.
- In comment mode, JSON records `energy_source: trajectory_comment` and null `mlip_backend`, `mlip_model`, `mlip_precision`, `charge`, and `multiplicity`. Recomputed output records `energy_source: mlip_recomputed` and the resolved values.
- With `--out-json`, `result.json` and `summary.json` contain identical payloads.
- In that payload, consume the ordered `output_files` list. The legacy
  basename-keyed `files` map is retained for compatibility and cannot represent
  two outputs with the same basename in different directories.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | XYZ trajectory whose per-frame comment line stores energies. | Required |
| `-o, --out PATH` | Repeatable output filenames; supports `.png`, `.jpg`/`.jpeg`, `.html`, `.svg`, `.pdf`, `.csv`. | `energy.png` |
| _extra arguments_ | Positional filenames listed after options; merged with the `-o` list. | _None_ |
| `--unit {kcal,hartree}` | Target unit for the plotted/exported values. | `kcal` |
| `-r, --reference TEXT` | Reference specification (`init`, `None`, or 0-based integer). | `init` |
| `-q, --charge INT` | Total charge used for MLIP recomputation. Triggers recomputation when supplied. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) used for MLIP recomputation. Triggers recomputation when supplied. | _None_ |
| `-b, --backend {uma,orb,mace,aimnet2}` | MLIP backend used only for recomputation. | `uma` |
| `--backend-model TEXT` | Model variant for the selected backend. | Backend default |
| `--precision {fp32,fp64}` | Backend-neutral recomputation precision; values are case-insensitive. | Backend default |
| `--out-json/--no-out-json` | Write `result.json` beside the first output. | `False` |
| `--reverse-x/--no-reverse-x` | Reverse the x-axis so the last frame appears on the left (and `init` becomes the last frame). | `False` |

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [path-search](path-search.md) — Recursive MEP search (produces XYZ trajectories suitable for trj2fig)
- [irc](irc.md) — IRC from TS (produces trajectories for energy profiling)
- [all](all.md) — End-to-end workflow
