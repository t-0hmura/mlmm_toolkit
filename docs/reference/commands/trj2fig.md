# `mlmm trj2fig`

```text
Usage: mlmm trj2fig [OPTIONS] [EXTRA_OUTS]...

  Plot ΔE or E from an XYZ trajectory and export figure/CSV.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                XYZ trajectory file  [required]
  -o, --out FILE                  Output file(s). You can repeat -o, and/or list
                                  extra filenames after options
                                  (.png/.jpg/.jpeg/.html/.svg/.pdf/.csv). If
                                  nothing is given, defaults to energy.png.
  --unit [kcal|hartree]           Energy unit.  [default: kcal]
  -r, --reference TEXT            Reference: "init" (initial frame; last frame
                                  if --reverse-x), "None" (absolute E), or an
                                  integer index.  [default: init]
  -q, --charge INTEGER            Total charge. Recompute energies when
                                  supplied.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1). Recompute energies
                                  when supplied.  [default: (1)]
  --reverse-x / --no-reverse-x    Reverse the x-axis (last frame on the left).
                                  [default: no-reverse-x]
  -b, --backend [uma|orb|mace|aimnet2]
                                  MLIP backend used when energies are
                                  recomputed.  [default: uma]
  --backend-model TEXT            Model variant for the selected backend;
                                  defaults to its built-in model.  [default:
                                  (the selected backend's own model)]
  --precision [fp32|fp64]         Backend-neutral precision used when energies
                                  are recomputed.  [default: (per backend: uma
                                  fp32; orb, mace fp64)]
  --out-json / --no-out-json      Write machine-readable result.json next to the
                                  first output.  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
