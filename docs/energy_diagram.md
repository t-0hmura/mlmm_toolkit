# `energy-diagram`

## Overview
`energy-diagram` draws an energy diagram from numeric values only.

- No structure input
- No energy calculation
- No `--thermo` / `--dft`

`-i/--input` accepts either:
- multiple numeric arguments
- a list-like numeric string

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

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT ...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT ...` | X-axis state labels. Count must match input value count. | `S1, S2, ...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |

## Notes
- Input order is used directly as plotting order.
- At least two numeric values are required.

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide
- [trj2fig](trj2fig.md) -- Plot profile from trajectory energies
- [all](all.md) -- End-to-end workflow with built-in energy diagram output
