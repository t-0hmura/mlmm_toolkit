# `energy-diagram`

## Overview

> **Summary:** Draw a state energy diagram directly from numeric values (no structure file, no ML/MM calculation).

`mlmm energy-diagram` only visualizes numbers you provide. It does not read PDB/XYZ structures and does not run thermochemistry (`--thermo`) or DFT (`--dft`) steps.

## Minimal example

```bash
mlmm energy-diagram -i 0 12.5 4.3 -o energy.png
```

## Output checklist

- `OUTPUT.(png|jpg|jpeg|svg|pdf)` -- the rendered energy diagram image

## Common examples

1. List string input.

```bash
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png
```

2. X/Y labels.

```bash
mlmm energy-diagram -i 0 12.5 4.3 --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## Workflow
1. Collect values from `-i/--input` (supports repeated flags, multiple values after one flag, and list-like strings).
2. Parse all input values as floats and fail early if fewer than two values are provided.
3. Parse optional `--label-x` values. If omitted, labels are auto-generated as `S1`, `S2`,...
4. Validate label count (`--label-x`) against value count, then render the diagram.
5. Save the image to `-o/--output` and print the saved path.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT...` | X-axis state labels. Count must match input value count. | `S1, S2,...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [trj2fig](trj2fig.md) -- Plot profile from trajectory energies
- [all](all.md) -- End-to-end workflow with built-in energy diagram output
