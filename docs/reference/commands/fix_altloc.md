# mlmm fix-altloc

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli fix-altloc [OPTIONS]

  Blank PDB altLoc column (col 17) without shifting, and keep one altLoc per
  atom by default rule: highest occupancy, then earliest appearance.

Options:
  --help-advanced               Show all options (including advanced settings)
                                and exit.
  -i, --input PATH              Input PDB file or directory.  [required]
  -o, --out PATH                Output file (if input is a file) or output
                                directory (if input is a directory).
  --recursive / --no-recursive  When input is a directory, process *.pdb
                                recursively (including subdirectories).
                                [default: no-recursive]
  --inplace / --no-inplace      Overwrite input file(s) in place (creates .bak
                                next to each file).  [default: no-inplace]
  --overwrite / --no-overwrite  Allow overwriting existing output files.
                                [default: no-overwrite]
  --force / --no-force          Process files even if no altLoc is detected
                                (default: skip files without altLoc).  [default:
                                no-force]
  -h, --help                    Show this message and exit.
```
