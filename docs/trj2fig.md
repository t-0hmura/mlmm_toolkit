# `trj2fig`

## Overview

> **Summary:** Extract energies from XYZ trajectory comment lines (or recompute with UMA), compute relative or absolute energy profiles, and export Plotly figures and CSV tables.

### Quick reference
- **Input:** An XYZ trajectory whose second line stores Hartree energies, or recompute energies with `-q/--charge` and/or `-m/--multiplicity`.
- **Reference modes:** first frame (`init`), no reference (`None`), or an explicit 0-based frame index.
- **Output formats:** PNG (default), JPEG, HTML, SVG, PDF, CSV.
- **Units:** kcal/mol (default) or Hartree.
- **X-axis flip:** `--reverse-x` reverses the axis so the last frame appears on the left.

`mlmm trj2fig` reads the Hartree energies encoded in each frame's comment line of an XYZ trajectory, converts them to kcal/mol or Hartree, optionally references all values to a chosen frame, and exports the resulting series as static/interactive figures and CSV tables. The figure uses bold ticks, consistent fonts, markers, and a smoothed spline curve (no title).

## Usage
```bash
mlmm trj2fig -i TRAJECTORY.xyz [-o OUTPUTS...] [-r REFERENCE] [--unit {kcal|hartree}] \
 [-q CHARGE] [-m MULTIPLICITY] [--reverse-x]
```

### Examples
```bash
# Default PNG, relative energy with respect to the first frame
mlmm trj2fig -i traj.xyz

# CSV + SVG with reference frame #5, reported in Hartree
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree

# Multiple outputs in one run with x-axis reversed
mlmm trj2fig -i traj.xyz -o energy.png energy.html energy.pdf --reverse-x
```

## Workflow
1. Parse the XYZ trajectory. By default, Hartree energies are extracted from each frame's comment line.
 If `-q/--charge` or `-m/--multiplicity` is provided, energies are recomputed with UMA (`uma-s-1p1`) instead.
2. Normalize the reference specification:
 - `init` -- frame `0` (or the last frame when `--reverse-x` is active).
 - `None`/`none`/`null` -- absolute energies (no referencing).
 - Integer literal -- the corresponding 0-based frame index.
3. Convert energies to either kcal/mol (default) or Hartree and, when a
 reference is active, subtract the reference value to produce delta-E.
4. Build the Plotly figure (strong ticks, spline interpolation, markers, no
 title) and export it to every requested extension.
5. Optionally emit a CSV table with columns `frame`, `energy_hartree`, and the
 appropriate delta-E or absolute-E column in the requested unit.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | XYZ trajectory whose second line stores energies. | Required |
| `-o, --out PATH` | Repeatable output filenames; supports `.png`, `.jpg`/`.jpeg`, `.html`, `.svg`, `.pdf`, `.csv`. | `energy.png` |
| _extra arguments_ | Positional filenames listed after options; merged with the `-o` list. | _None_ |
| `--unit {kcal,hartree}` | Target unit for the plotted/exported values. | `kcal` |
| `-r, --reference TEXT` | Reference specification (`init`, `None`, or 0-based integer). | `init` |
| `-q, --charge INT` | Total charge used for UMA recomputation. Triggers recomputation when supplied. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) used for UMA recomputation. Triggers recomputation when supplied. | _None_ |
| `--reverse-x` | Reverse the x-axis so the last frame appears on the left (and `init` becomes the last frame). | `False` |

## Outputs
```
<output>.[png|jpg|jpeg|html|svg|pdf] # Plotly export for every requested extension (defaults to energy.png)
<output>.csv # Optional energy table when CSV is requested
```
- When no `-o` or positional outputs are provided, a single `energy.png` is written
 to the current directory.
- CSV exports include `frame`, `energy_hartree`, and either a delta-E column
 (`delta_kcal`/`delta_hartree`) or absolute column (`energy_kcal`/`energy_hartree`
 when no reference is applied).
- PNG uses Plotly's PNG export with `scale=2` for higher resolution.

## Notes
- Energies are taken from the first decimal number in each comment; malformed
 comments raise an error.
- When recomputation is triggered, omitted charge/multiplicity values default to `0` and `1`.
- Unsupported file extensions in `-o` cause an error.
- `--reverse-x` flips both the axis direction and the behavior of `-r init` so the plotted pathway is read in reverse direction.
- The `--output-peak` option has been removed.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [path_search](path_search.md) -- Recursive MEP search (produces XYZ trajectories suitable for trj2fig)
- [irc](irc.md) -- IRC from TS (produces trajectories for energy profiling)
- [all](all.md) -- End-to-end workflow
