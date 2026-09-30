# `mlmm add-elem-info`

```text
Usage: mlmm add-elem-info [OPTIONS]

  Add/repair element columns (77–78) in a PDB.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input PDB filepath  [required]
  -o, --output, --out FILE        Output PDB filepath (default: replace ".pdb"
                                  with "_add_elem.pdb"; when provided,
                                  --overwrite is ignored unless this is the
                                  input file, which requires it).
  --overwrite / --no-overwrite    Overwrite the input file in-place when
                                  -o/--out is omitted.  [default: no-overwrite]
  --overwrite-elem / --no-overwrite-elem
                                  Also re-infer element columns that already
                                  hold a valid element (off: only empty or
                                  invalid columns are filled).  [default: no-
                                  overwrite-elem]
  -h, --help                      Show this message and exit.
```
