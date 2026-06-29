# `mlmm scan2d`

```text
Usage: mlmm scan2d [OPTIONS]

  2D distance scan with harmonic restraints using the ML/MM calculator.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
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
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  --scan-lists.  [default: no-print-parsed]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
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
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based (default) or 0-based.  [default: one-
                                  based]
  --max-step-size FLOAT           Maximum spacing between successive distance
                                  targets [Å].  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2]. Defaults to
                                  YAML bias.k (BIAS_KW['k']=300 in defaults.py)
                                  when omitted; explicit CLI value overrides
                                  YAML.
  --relax-max-cycles INTEGER      Maximum LBFGS cycles per biased relaxation
                                  (also used for preopt).  [default: 10000]
  --dump / --no-dump              Write inner d2 scan TRJs per d1 slice.
                                  [default: no-dump]
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan2d/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset.  [default: baker]
  --ref-pdb FILE                  Reference PDB topology to use when --input is
                                  XYZ (keeps XYZ coordinates).
  --preopt / --no-preopt          Run an unbiased pre-optimization.  [default:
                                  no-preopt]
  --baseline [min|first]          Reference for relative energy (kcal/mol):
                                  'min' or 'first' (i=0,j=0).  [default: min]
  --zmin FLOAT                    Lower bound of the color scale (kcal/mol).
  --zmax FLOAT                    Upper bound of the color scale (kcal/mol).
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [x>=1]
  --precision [fp32|fp64]         MLIP backend precision: fp32 (default) or
                                  fp64. Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p1 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / MACE-OFF23_small for mace). Default: the
                                  backend's built-in model.
  --deterministic / --no-deterministic
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). For an intentional open-shell or
                                  covalently-cut ML region.
  -h, --help                      Show this message and exit.
```
