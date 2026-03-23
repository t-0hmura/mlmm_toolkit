# mlmm opt

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli opt [OPTIONS]

  ML/MM geometry optimization with LBFGS (light) or RFO (heavy).

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (PDB, XYZ). XYZ provides
                                  higher coordinate precision. If XYZ, use
                                  --ref-pdb to specify PDB topology for atom
                                  ordering and output conversion.  [required]
  --parm FILE                     Amber parm7 topology covering the whole enzyme
                                  complex.  [required]
  --model-pdb FILE                PDB defining atoms that belong to the ML
                                  (high-level) region. Optional when --detect-
                                  layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.
  -o, --out-dir TEXT              Output directory.  [default: ./result_opt/]
  --opt-mode [grad|hess|light|heavy|lbfgs|rfo]
                                  Optimization mode: grad (lbfgs) or hess (rfo).
                                  Aliases light/heavy and lbfgs/rfo are
                                  accepted.  [default: grad]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
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
