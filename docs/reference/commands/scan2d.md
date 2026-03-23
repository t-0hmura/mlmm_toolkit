# `mlmm scan2d`

```text
mlmm-toolkit ver. 0.2.5.dev18

Usage: mlmm scan2d [OPTIONS]

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
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            ML-region total charge. Required unless
                                  --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., "1,3,5").
  --hess-cutoff FLOAT             Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer.
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path.
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based (default) or 0-based.  [default: one-
                                  based]
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  --scan-lists.  [default: no-print-parsed]
  --max-step-size FLOAT           Maximum spacing between successive distance
                                  targets [Å].  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2].  [default:
                                  300.0]
  --relax-max-cycles INTEGER      Maximum LBFGS cycles per biased relaxation
                                  (also used for preopt).  [default: 10000]
  --dump / --no-dump              Write inner d2 scan TRJs per d1 slice.
                                  [default: no-dump]
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan2d/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset.  [default: baker]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --ref-pdb FILE                  Reference PDB topology to use when --input is
                                  XYZ (keeps XYZ coordinates).
  --preopt / --no-preopt          Run an unbiased pre-optimization.  [default:
                                  no-preopt]
  --baseline [min|first]          Reference for relative energy (kcal/mol):
                                  'min' or 'first' (i=0,j=0).  [default: min]
  --zmin FLOAT                    Lower bound of the contour color scale
                                  (kcal/mol).
  --zmax FLOAT                    Upper bound of the contour color scale
                                  (kcal/mol).
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from ML region for MM
                                  point charges in xTB embedding. Default: 12.0
                                  Å. Only used when --embedcharge is enabled.
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  --mm-backend [hessian_ff|openmm]
                                  MM backend: hessian_ff (analytical Hessian,
                                  default) or openmm (finite-difference Hessian,
                                  slower).
  --cmap / --no-cmap              Enable CMAP (backbone cross-map) terms in
                                  model parm7. Default: disabled (Gaussian
                                  ONIOM-compatible).
  -h, --help                      Show this message and exit.
```
