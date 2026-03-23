# mlmm trj2fig

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli trj2fig [OPTIONS] [EXTRA_OUTS]...

  Plot ΔE or E from an XYZ trajectory and export figure/CSV.

Options:
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  -i, --input FILE            XYZ trajectory file  [required]
  -o, --out FILE              Output file(s). You can repeat -o, and/or list
                              extra filenames after options
                              (.png/.html/.svg/.pdf/.csv). If nothing is given,
                              defaults to energy.png.
  --unit [kcal|hartree]       Energy unit.
  -r, --reference TEXT        Reference: "init" (initial frame; last frame if
                              --reverse-x), "None" (absolute E), or an integer
                              index.
  -q, --charge INTEGER        Total charge. Recompute energies when supplied.
  -m, --multiplicity INTEGER  Spin multiplicity (2S+1). Recompute energies when
                              supplied.
  -h, --help                  Show this message and exit.
```
