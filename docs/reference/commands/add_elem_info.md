# `mlmm add-elem-info`

```text
Usage: mlmm add-elem-info [OPTIONS]

  Add/repair element columns (77–78) in a PDB using Biopython.

Options:
  -v, --verbose LEVEL           Console verbosity 0-3 (default 2). 0=silent;
                                1=milestones only; 2=+detailed step logging and
                                deliverable paths; 3=everything (full config
                                blocks, per-file paths, DEBUG logging).
                                [0<=x<=3]
  --help-advanced               Show all options (including advanced settings)
                                and exit.
  -i, --input FILE              Input PDB filepath  [required]
  -o, --out FILE                Output PDB filepath (default:
                                <input>_add_elem.pdb; overrides --inplace)
  --inplace / --no-inplace      Replace the input file when -o/--out is omitted.
                                [default: no-inplace]
  --overwrite / --no-overwrite  Re-infer and overwrite element fields even if
                                present (by default, existing values are
                                preserved).  [default: no-overwrite]
  -h, --help                    Show this message and exit.
```
