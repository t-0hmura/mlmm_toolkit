# mlmm freq

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli freq [OPTIONS]

  ML/MM vibrational frequency analysis (PHVA-compatible).

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Enzyme complex PDB used by both geom_loader
                                  and the ML/MM calculator.  [required]
  --parm FILE                     Amber parm7 topology for the full enzyme
                                  complex.  [required]
  --model-pdb FILE                PDB defining atoms belonging to the ML region.
                                  Optional when --detect-layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.
  --temperature FLOAT             Temperature (K) for thermochemistry summary.
                                  [default: 298.15]
  --pressure FLOAT                Pressure (atm) for thermochemistry summary.
                                  [default: 1.0]
  -o, --out-dir TEXT              Output directory.  [default: ./result_freq/]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  -h, --help                      Show this message and exit.
```
