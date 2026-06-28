# `energy-diagram`

Draw a state energy diagram directly from numeric values (no structure file, no ML/MM calculation). It does not read PDB/XYZ structures and does not run thermochemistry (`--thermo`) or DFT (`--dft`) steps.

## Examples

Relative energies as a list-like string:

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" -o energy.png
```

Absolute energies (any numeric values):

```bash
mlmm energy-diagram -i "[-205.1, -190.4, -198.7]" -o energy.png
```

Values via repeated `-i` flags:

```bash
mlmm energy-diagram -i 0 -i 12.5 -i 4.3 -o energy.png
```

Custom x-axis state labels and y-axis label:

```bash
mlmm energy-diagram -i "[0, 12.5, 4.3]" --label-x R TS P --label-y "ΔE (kcal/mol)" -o energy.png
```

## Workflow
1. Collect values from `-i/--input` (supports repeated flags, multiple values after one flag, and list-like strings).
2. Parse all input values as floats and fail early if fewer than two values are provided.
3. Parse optional `--label-x` values. If omitted, labels are auto-generated as `S1`, `S2`,...
4. Validate label count (`--label-x`) against value count, then render the diagram.
5. Save the image to `-o/--output` and print the saved path.

## Outputs

- `OUTPUT.(png|jpg|jpeg|svg|pdf)` -- the rendered energy diagram image

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input TEXT...` | Numeric values (multiple args or list-like string). | Required |
| `-o, --output PATH` | Output image path (`.png/.jpg/.jpeg/.svg/.pdf`). | `energy_diagram.png` |
| `--label-x TEXT...` | X-axis state labels. Count must match input value count. | `S1, S2,...` |
| `--label-y TEXT` | Y-axis label. | `ΔE (kcal/mol)` |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` next to the output image. | `--no-out-json` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [trj2fig](trj2fig.md) — Plot profile from trajectory energies
- [all](all.md) — End-to-end workflow with built-in energy diagram output
