# `mlmm irc`

```text

Usage: mlmm irc [OPTIONS]

  Run an IRC calculation with EulerPC. Only the documented CLI options are
  accepted; all other settings come from YAML.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  etc.).  [required]
  --parm FILE                     Amber parm7 topology for the whole enzyme (MM
                                  region). If omitted, must be provided in YAML
                                  as calc.real_parm7.
  --model-pdb FILE                PDB defining atoms belonging to the ML region.
                                  Optional when --detect-layer is enabled.
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
  -q, --charge INTEGER            Total charge; overrides calc.charge from YAML.
                                  Required unless --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1); overrides calc.spin
                                  from YAML.
  --max-cycles INTEGER            Maximum number of IRC steps; overrides
                                  irc.max_cycles from YAML.
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates). Default: 0.10 Bohr. Overrides
                                  irc.step_length from YAML.
  --root INTEGER                  Imaginary mode index used for the initial
                                  displacement; overrides irc.root from YAML.
  --forward / --no-forward        Run the forward IRC; overrides irc.forward
                                  from YAML.
  --backward / --no-backward      Run the backward IRC; overrides irc.backward
                                  from YAML.
  -o, --out-dir TEXT              Output directory; overrides irc.out_dir from
                                  YAML.  [default: ./result_irc/]
  --hessian-calc-mode [analytical|finitedifference]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); overrides
                                  calc.hessian_calc_mode from YAML. Default:
                                  'FiniteDifference'. Use 'Analytical' when VRAM
                                  is sufficient.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running IRC.  [default: no-dry-run]
  --ref-pdb FILE                  Reference PDB topology to use when --input is
                                  XYZ (keeps XYZ coordinates).
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
  --hess-device [auto|cuda|cpu]   Device for initial Hessian storage and IRC
                                  operations (auto/cuda/cpu). Use 'cpu' for
                                  large unfrozen systems to avoid VRAM limits.
                                  [default: auto]
  --read-hess FILE                Read initial Hessian from a .npz file
                                  (produced by 'mlmm freq --dump-hess'). Takes
                                  priority over the hessian_cache and fresh
                                  computation.
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
