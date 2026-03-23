# mlmm define-layer

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli define-layer [OPTIONS]

  Define 3-layer ML/MM system based on distance from ML region.

Options:
  --help-advanced       Show all options (including advanced settings) and exit.
  -i, --input FILE      Input PDB file containing the full system.  [required]
  --model-pdb FILE      PDB file defining atoms in the ML region.
  --model-indices TEXT  Comma-separated atom indices for ML region (e.g.,
                        '1,2,3,4' or '1-10,15,20-25'). Takes precedence over
                        --model-pdb.
  -o, --output FILE     Output PDB file with B-factors set to layer values.
                        Defaults to '<input>_layered.pdb'.
  -h, --help            Show this message and exit.
```
