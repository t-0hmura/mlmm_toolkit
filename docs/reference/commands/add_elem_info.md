# `mlmm add-elem-info`

```text
mlmm-toolkit ver. 0.2.5.dev18

Usage: mlmm add-elem-info [OPTIONS]

  Add/repair element columns (77–78) in a PDB using Biopython.

Options:
  --help-advanced               Show all options (including advanced settings)
                                and exit.
  -i, --input FILE              Input PDB filepath  [required]
  -o, --out FILE                Output PDB filepath (omit to overwrite input)
  --overwrite / --no-overwrite  Re-infer and overwrite element fields even if
                                present (by default, existing values are
                                preserved).  [default: no-overwrite]
  -h, --help                    Show this message and exit.
```
