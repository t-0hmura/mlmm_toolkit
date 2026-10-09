# `trj2fig` (energy profile of a trajectory)

`trj2fig` **plots the energy along an XYZ trajectory**, such as one written by `opt`, `scan`, `path-opt`, `path-search`, or `irc`. It reads the energy of each frame from its comment line and exports figures and a CSV table.

---

## What it is for

* **Energy profiles of paths**: minimum energy path (MEP), scan, and intrinsic reaction coordinate (IRC) trajectories as ΔE plots.
* **Checking an optimization**: how the energy fell over the optimization cycles.
* **Data for your own plots**: per-frame energies as CSV.

---

## Examples

### 1. Default PNG

Plot ΔE relative to the first frame and write `energy.png`.

```bash
mlmm trj2fig -i traj.xyz --out-json
```

The console prints `[trj2fig] Saved figure -> energy.png`, and `result.json` has `n_frames`, the number of frames read.

### 2. CSV and SVG relative to frame 5, in hartree

Write a table and a figure in hartree, relative to frame index 5 (counted from 0, so the sixth frame).

```bash
mlmm trj2fig -i traj.xyz -o energy.csv energy.svg -r 5 --unit hartree
```

### 3. Several formats, x-axis reversed

Put the last frame on the left; ΔE is then relative to the last frame.

```bash
mlmm trj2fig -i traj.xyz --reverse-x -o energy.png energy.html energy.pdf
```

### 4. Recompute energies with an MLIP

Recompute every frame as a neutral singlet with a chosen [UMA](backends.md) model and precision instead of reading the comments, and record the backend in `result.json`.

```bash
mlmm trj2fig -i traj.xyz -q 0 -m 1 -b uma --backend-model uma-s-1p2 \
    --precision fp32 -o energy.png energy.csv --out-json
```

---

## How it works

1. **Reading the energies**:
Each frame's comment line gives its energy. Trajectories written by mlmm, such as `optimization_trj.xyz` (`opt --dump`), `scan_trj.xyz`, `mep_trj.xyz`, and `finished_irc_trj.xyz`, are read as they are. In other files, write `E=<value>` with an optional unit (`Ha`, `Eh`, `hartree`, `eV`, `kcal/mol`); without a unit it is hartree, except `energy=` in extended XYZ (a comment with `Properties=` or `Lattice=`), which is eV. With `-q` or `-m`, the MLIP (machine-learning interatomic potential) backend recomputes every frame instead.
2. **Choosing the reference**:
`-r init` (the default) is the frame at the left end: the first frame, or the last with `--reverse-x`; an integer is a 0-based frame index; `none` plots absolute energies.
3. **Converting the unit**:
The energies are converted to kcal/mol (default) or hartree, and the reference is subtracted to give ΔE. The y-axis reads `ΔE (kcal/mol)`, or `E (…)` for absolute energies.
4. **Exporting**:
Each output gets its format from its extension: `.png`, `.jpg`, `.jpeg`, `.svg`, `.pdf`, and `.html` are figures, and `.csv` is a table. PNG is written at twice the pixel width and height of the figure.

---

## Output files

```text
energy.png      # Figure (default when no output is given)
energy.csv      # Energy table (when a .csv output is given)
result.json     # Summary (with --out-json)
summary.json    # Copy of result.json; read result.json (with --out-json)
```

* **CSV columns**: `frame`, `energy_hartree`, and the plotted value in the `--unit` unit. The third column is named `delta_kcal` or `delta_hartree` with a reference, and `energy_kcal` or `energy_hartree` with `-r none`.
* **`result.json`** is written in the directory of the first output. It has `n_frames`, `min_energy_hartree`, `max_energy_hartree`, `energy_source` (`trajectory_comment` or `mlip_recomputed`), and `output_files`, the written files in order. `mlip_backend`, `mlip_model`, `mlip_model_label`, `mlip_task`, and `mlip_precision` are filled only when the energies were recomputed, and are null otherwise. See [JSON Output Reference](json-output.md#trj2fig).

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | XYZ trajectory |
| `-o, --output` | path | `energy.png` | Output files (`.png`, `.jpg`, `.jpeg`, `.html`, `.svg`, `.pdf`, `.csv`); repeat `-o`, or list more file names after it (`-o energy.csv energy.svg`) |
| `--unit` | `kcal` / `hartree` | `kcal` | Unit of the plotted and exported values |
| `-r, --reference` | text | `init` | Reference: `init`, `none`, or a 0-based frame index |
| `-q, --charge` | integer | `None` | Total charge of the frame; recomputes the energies with the MLIP when given |
| `-m, --multiplicity` | integer | `None` | Spin multiplicity (2S+1); recomputes the energies with the MLIP when given |
| `--reverse-x/--no-reverse-x` | flag | `False` | Put the last frame on the left |
| `-b, --backend` | text | `uma` | MLIP for recomputation (`uma`, `orb`, `mace`, `aimnet2`; see [MLIP Backends](backends.md)) |
| `--backend-model` | text | `None` | Model of the selected backend (e.g. `uma-s-1p2`); without it, the backend's default model |
| `--precision` | `fp32` / `fp64` | per backend | Precision of the recomputation (UMA `fp32`; ORB and MACE `fp64`); AIMNet2 accepts only `fp32` |
| `--out-json/--no-out-json` | flag | `False` | Write `result.json` and `summary.json` |

See the [generated CLI reference](reference/commands/trj2fig.md) for every option.

---

## Notes

* **Recomputation**: with only one of `-q` and `-m`, the other is taken as charge 0 or multiplicity 1.
* **Recomputed energies are MLIP energies**: the MLIP computes every atom of the frame, with no ML/MM split, so the result is not the ML/MM ONIOM energy. To plot the ONIOM energy of an mlmm trajectory, read its comment lines (no `-q` or `-m`).
* **Comments without a readable energy** stop the run with an error that names the frame.
* **Unsupported extensions** stop the run with an error.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [path-search](path-search.md) — MEP trajectories to plot
* [irc](irc.md) — IRC trajectories to plot
* [energy-diagram](energy-diagram.md) — a state energy diagram from numbers you give
* [all](all.md) — the full workflow
* [Troubleshooting](troubleshooting.md) — what to do when a run fails; for a failed figure export, see {ref}`Installation / environment <installation-environment-problems>`
