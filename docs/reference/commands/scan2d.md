# mlmm scan2d

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli scan2d [OPTIONS]

  2D distance scan with harmonic restraints using the ML/MM calculator.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input enzyme complex PDB (required).
                                  [required]
  --parm FILE                     Amber parm7 topology for the enzyme
                                  (required).  [required]
  --model-pdb FILE                PDB defining the ML region. Optional when
                                  --detect-layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            ML-region total charge. Required unless
                                  --ligand-charge is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path.
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan2d/]
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
