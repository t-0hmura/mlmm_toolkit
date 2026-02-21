# `energy-diagram`

## Overview

> **Summary:** Draw a state energy diagram directly from numeric values (no structure file, no ML/MM calculation).

### Quick reference
- **Input:** Numeric energies passed via `-i/--input` (multiple tokens or list-like string).
- **Output:** One image file (`.png`, `.jpg`, `.jpeg`, `.svg`, or `.pdf`).
- **Default output:** `energy_diagram.png`.
- **State labels:** Optional `--label-x`; defaults to `S1`, `S2`, ...
- **Use when:** You already have energies and only need a formatted diagram.

`mlmm energy-diagram` only visualizes numbers you provide. It does not read PDB/XYZ structures and does not run thermochemistry (`--thermo`) or DFT (`--dft`) steps.

## Usage
```bash
mlmm energy-diagram -i VALUES... [-o OUTPUT] [--label-x ...] [--label-y ...]
```

## Examples
```bash
# Multiple numeric arguments
mlmm energy-diagram -i 0 12.5 4.3 -o energy.png

# List string
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png

# X/Y labels
mlmm energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## Workflow
1. Collect values from `-i/--input` (supports repeated flags, multiple values after one flag, and list-like strings).
2. Parse all input values as floats and fail early if fewer than two values are provided.
3. Parse optional `--label-x` values. If omitted, labels are auto-generated as `S1`, `S2`, ...
4. Validate label count (`--label-x`) against value count, then render the diagram.
5. Save the image to `-o/--output` and print the saved path.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT ...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT ...` | X-axis state labels. Count must match input value count. | `S1, S2, ...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |

## Outputs
```
OUTPUT.(png|jpg|jpeg|svg|pdf)
```
- If `-o/--output` is omitted, `energy_diagram.png` is written to the current directory.
- When output has no extension, `.png` is appended automatically.
- Parent directories are created automatically when needed.

## Notes
- Input order is used directly as plotting order.
- At least two numeric values are required.
- This command does not read structure files and does not compute energies.

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [trj2fig](trj2fig.md) -- Plot profile from trajectory energies
- [all](all.md) -- End-to-end workflow with built-in energy diagram output
