# `mlmm energy-diagram`

```text
Usage: mlmm energy-diagram [OPTIONS]

  Plot an energy diagram from numeric inputs only. Supply values via repeated -i
  or a single list-like string.

Options:
  -v, --verbose LEVEL         Console verbosity 0-3 (default 2). 0=silent;
                              1=milestones only; 2=+optimizer cycle tables, per-
                              stage timing, VRAM, deliverable paths;
                              3=everything (full config blocks, per-file paths,
                              DEBUG logging).  [0<=x<=3]
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  -i, --input TEXT            Numeric sequence. Give one value per -i (-i 0 -i
                              12.5 -i 4.3), or a single list-like string (-i
                              "[0, 12.5, 4.3]").  [required]
  -o, --output FILE           Output image path.  [default: energy_diagram.png]
  --label-x TEXT              State labels on x-axis. Give one label per
                              --label-x (--label-x R --label-x TS --label-x P),
                              or a single list-like string (--label-x
                              "['R','TS','P']").
  --label-y TEXT              Y-axis label.  [default: ΔE (kcal/mol)]
  --out-json / --no-out-json  Write machine-readable result.json next to the
                              output image.  [default: no-out-json]
  -h, --help                  Show this message and exit.
```
