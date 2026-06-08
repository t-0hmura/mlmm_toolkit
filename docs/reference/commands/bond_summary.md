# `mlmm bond-summary`

```text
Usage: mlmm bond-summary [OPTIONS] [EXTRA_INPUTS]...

  Detect bond changes between consecutive structures.

Options:
  -v, --verbose LEVEL         Console verbosity 0-3 (default 2). 0=silent;
                              1=milestones only; 2=+optimizer cycle tables, per-
                              stage timing, VRAM, deliverable paths;
                              3=everything (full config blocks, per-file paths,
                              DEBUG logging).  [0<=x<=3]
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  -i, --input TEXT            Input structure files (XYZ/PDB/GJF). Repeat -i for
                              each file.
  --device TEXT               Compute device for distance calculations.
                              [default: cpu]
  --bond-factor FLOAT         Scaling factor for covalent radii sum.  [default:
                              1.2]
  --one-based / --zero-based  Use 1-based atom indices in output.  [default:
                              one-based]
  --json / --no-json          Write the JSON report to stdout; no file is
                              created.  [default: no-json]
  -h, --help                  Show this message and exit.
```
