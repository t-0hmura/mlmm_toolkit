# `mlmm define-layer`

```text

Usage: mlmm define-layer [OPTIONS]

  Define 3-layer ML/MM system based on distance from ML region.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input PDB file containing the full system.
                                  [required]
  --model-pdb FILE                PDB file defining atoms in the ML region.
  --model-indices TEXT            Comma-separated atom indices for ML region
                                  (e.g., '1,2,3,4' or '1-10,15,20-25'). Takes
                                  precedence over --model-pdb.
  --radius-freeze FLOAT           Distance cutoff (Å) from ML region for movable
                                  MM atoms. Atoms beyond this distance are
                                  frozen.  [default: 8.0]
  -o, --output FILE               Output PDB file with B-factors set to layer
                                  values. Defaults to '<input>_layered.pdb'.
  --one-based / --zero-based      Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: one-based]
  -h, --help                      Show this message and exit.
```
