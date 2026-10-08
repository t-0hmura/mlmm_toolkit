# `energy-diagram` (state energy diagram)

`energy-diagram` **draws a state energy diagram** from numbers you give it. It reads no structure files and runs no calculation. It suits energies you already have, for example the {ref}`energy_diagrams <summary-json-path-search-all>` of the `summary.json` written by `all` or `path-search`, or a table in a paper. It saves the diagram as an image file.

---

## What it is for

* **Figures from known energies**: energies of the reactant (R), transition state (TS), intermediate (IM), and product (P) taken from a workflow or from DFT.
* **Figures for papers and slides**: vector output as SVG or PDF.
* **Quick checks**: a diagram of a few values without writing plotting code.

---

## Examples

### 1. Values as one list

Pass all values as one quoted list.

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" -o energy.png --out-json
```

The console prints `[energy-diagram] Saved -> …/energy.png`, and `result.json` next to the image has `n_points: 3`.

### 2. One `-i` per value

Repeat `-i` once for each value.

```bash
mlmm energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

### 3. State and axis labels

Name the states on the x-axis and set the y-axis label.

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" \
  --label-x "['R','TS','P']" --label-y "ΔE (kcal/mol)" -o energy.png
```

---

## How it works

1. **Reading the values**:
Values come from `-i`, repeated once per value or given as one list-like string (`"[0, 12.5, 4.3]"` or `"0, 12.5, 4.3"`). Any numbers are drawn as given, so absolute energies work as well as relative ones.
2. **Labels**:
`--label-x` gives one label per state, repeated or as one list-like string. Without it the states are named `S1`, `S2`, ….
3. **Drawing**:
Each state is drawn as a short horizontal bar at its energy, neighboring bars are joined by dotted lines, and a light gray dotted line marks the energy of the first state.
4. **Saving**:
The extension of `-o` chooses the format. A path without an extension gets `.png`, and missing parent directories are created.

---

## Output files

```text
energy_diagram.png   # The diagram (default name; set with -o)
result.json          # execution_status, scientific_status, n_points, and files (with --out-json)
summary.json         # Copy of result.json; read result.json (with --out-json)
```

`result.json` and `summary.json` are written in the directory of the image. They record the number of points and the image path, not the values or the labels.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | text | (required) | Energies: repeat `-i` once per value, or give one list-like string |
| `-o, --output` | path | `energy_diagram.png` | Output image (`.png`, `.jpg`, `.jpeg`, `.svg`, `.pdf`) |
| `--label-x` | text | `S1, S2, …` | State labels on the x-axis: repeat once per state, or give one list-like string |
| `--label-y` | text | `ΔE (kcal/mol)` | Y-axis label |
| `--out-json/--no-out-json` | flag | `False` | Write `result.json` and `summary.json` next to the image |

See the [generated CLI reference](reference/commands/energy_diagram.md) for every option.

---

## Notes

* **At least two values**: fewer stop the run with `Provide at least two numeric values with -i/--input.`
* **Label count**: the number of `--label-x` labels must equal the number of values.
* **Order**: the input order is the order on the x-axis.
* **Units**: the values are drawn unchanged, so state their unit in `--label-y`.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [trj2fig](trj2fig.md) — energy profile from the frames of a trajectory
* [all](all.md) — the full workflow, which draws its own energy diagrams
* [JSON Output Reference](json-output.md#energy-diagram) — the fields of `result.json`
* [Troubleshooting](troubleshooting.md) — what to do when a run fails; for a failed image export, see {ref}`Installation / environment <installation-environment-problems>`
