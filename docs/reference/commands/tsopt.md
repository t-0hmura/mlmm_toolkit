# `mlmm tsopt`

```text

Usage: mlmm tsopt [OPTIONS]

  TS optimization: grad (Dimer) or hess (RS-I-RFO) for the ML/MM calculator.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Starting geometry (PDB or XYZ). XYZ provides
                                  higher coordinate precision. If XYZ, use
                                  --ref-pdb to specify PDB topology for atom
                                  ordering and output conversion.  [required]
  --ref-pdb FILE                  Reference PDB topology when input is XYZ. XYZ
                                  coordinates are used (higher precision) while
                                  PDB provides atom ordering and residue
                                  information for output conversion.
  --parm FILE                     Amber parm7 topology for the whole enzyme (MM
                                  region).  [required]
  --model-pdb FILE                PDB containing the ML-region atoms. Optional
                                  when --detect-layer is enabled.
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
  -q, --charge INTEGER            Total charge of the ML region. Required unless
                                  --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
  --freeze-atoms TEXT             Comma-separated 1-based indices to freeze
                                  (e.g., '1,3,5').
  --radius-hessian, --hess-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms. Default 0.0 means
                                  ML-only partial Hessian.  [default: 0.0]
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.
  --hessian-calc-mode [analytical|finitedifference]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); overrides
                                  calc.hessian_calc_mode from YAML. Default:
                                  'FiniteDifference'. Use 'Analytical' when VRAM
                                  is sufficient.
  --max-cycles INTEGER            Maximum total optimization cycles.  [default:
                                  10000]
  --dump / --no-dump              Write concatenated trajectory
                                  'optimization_all_trj.xyz'.  [default: no-
                                  dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_tsopt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset.
  --opt-mode [grad|hess|light|heavy|dimer|rsirfo]
                                  grad (dimer) or hess (rsirfo). Aliases
                                  light/heavy and dimer/rsirfo are accepted.
                                  [default: hess]
  --microiter / --no-microiter    Enable microiteration: alternate ML 1-step
                                  (RS-I-RFO) and MM relaxation (LBFGS with MM-
                                  only forces). Only effective in --opt-mode
                                  hess. Ignored in grad mode.  [default:
                                  microiter]
  --partial-hessian-flatten / --full-hessian-flatten
                                  Use partial (active-block) Hessian for
                                  imaginary mode detection in flatten loop.
                                  [default: partial-hessian-flatten]
  --flatten / --no-flatten        Enable/disable extra imaginary-mode flattening
                                  loop. --flatten uses the default
                                  flatten_max_iter (50); --no-flatten forces it
                                  to 0. When not provided, the value is
                                  determined by the YAML config or defaults.
  --ml-only-hessian-dimer / --no-ml-only-hessian-dimer
                                  Use ML-region-only Hessian (no MM Hessian
                                  contribution) for dimer orientation in grad
                                  mode. Faster but less accurate for mode
                                  direction.  [default: no-ml-only-hessian-
                                  dimer]
  --active-dof-mode [all|ml-only|partial|unfrozen]
                                  Active DOF selection for final frequency
                                  analysis: all (all atoms), ml-only (ML only),
                                  partial (ML + MovableMM, default), unfrozen
                                  (all except frozen layer).  [default: partial]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running TS optimization.  [default:
                                  no-dry-run]
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
  --skip-final-freq / --no-skip-final-freq
                                  Skip the post-convergence frequency analysis
                                  and imaginary-mode flattening. Useful for
                                  large unfrozen systems where the final Hessian
                                  diagonalization is expensive.  [default: no-
                                  skip-final-freq]
  -h, --help                      Show this message and exit.
```
