# `mlmm oniom-import`

```text
Usage: mlmm oniom-import [OPTIONS]

  Import ONIOM input (Gaussian g16 or ORCA) and reconstruct XYZ + B-factor
  layered PDB.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input ONIOM file (.gjf/.com for g16, .inp for
                                  ORCA).  [required]
  --mode [g16|orca]               Input mode. If omitted, inferred from input
                                  suffix.
  -o, --out-prefix PATH           Output prefix. Defaults to input stem in the
                                  current working directory.
  --ref-pdb FILE                  Reference PDB to preserve atom naming/residue
                                  metadata (atom count must match).
  --allow-unverified-ref-order / --no-allow-unverified-ref-order
                                  Allow legacy positional --ref-pdb mapping when
                                  repeated elements make atom identity
                                  unverifiable. Use only after independently
                                  checking order.  [default: no-allow-
                                  unverified-ref-order]
  -h, --help                      Show this message and exit.
```
