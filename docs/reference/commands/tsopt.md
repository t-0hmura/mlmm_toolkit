# `mlmm tsopt`

```text
Usage: mlmm tsopt [OPTIONS]

  TS optimization: grad (Dimer) or hess (RS-I-RFO) for the ML/MM calculator.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Starting geometry (PDB/mmCIF or XYZ). XYZ
                                  provides higher coordinate precision. If XYZ,
                                  use --ref-pdb to specify PDB topology for atom
                                  ordering and output conversion.  [required]
  --ref-mode FILE                 Advanced path-mode hint for Hessian TS root
                                  selection (.npy or whitespace Cartesian 3N
                                  text). 'mlmm all' supplies this from its MEP;
                                  ordinary standalone tsopt runs normally omit
                                  it.
  --ref-pdb FILE                  Reference PDB topology when input is XYZ. XYZ
                                  coordinates are used (higher precision) while
                                  PDB provides atom ordering and residue
                                  information for output conversion.
  --parm FILE                     Amber parm7 topology for the whole enzyme (MM
                                  region).  [required]
  --model-pdb FILE                ML-only, link-H-free PDB subset; atom
                                  identity/order must match the full PDB/parm7.
                                  When provided, it defines ML membership;
                                  --detect-layer still reads valid
                                  movable/frozen MM B-factors.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  -q, --charge INTEGER            Total charge of the ML region. Required unless
                                  --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  [default: (1)]
  --freeze-atoms TEXT             Comma-separated 1-based indices to freeze
                                  (e.g., '1,3,5').
  --radius-hessian, --hess-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms. Unset includes
                                  every required movable MM atom; 0.0 requests
                                  an ML-only Hessian and should be paired with
                                  --active-dof-mode ml-only for final frequency
                                  validation.  [default: (all movable MM atoms)]
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.  [default: (use freeze_atoms)]
  --hessian-calc-mode [analytical|finitedifference]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); overrides
                                  calc.hessian_calc_mode from YAML. Default:
                                  'FiniteDifference'. Runtime and memory depend
                                  on the backend and system; compare both modes
                                  on a representative pilot.  [default:
                                  (FiniteDifference)]
  --max-cycles INTEGER            Maximum total optimization cycles.  [default:
                                  10000]
  --dump / --no-dump              Write concatenated trajectory
                                  'optimization_all_trj.xyz'.  [default: no-
                                  dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_tsopt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset.  [default: (baker)]
  --opt-mode [grad|hess|light|heavy|dimer|rsirfo|trim|rsprfo]
                                  grad/dimer/light → Hessian Guided Dimer;
                                  hess/rsirfo/heavy → RS-I-RFO; trim → TRIM
                                  (Helgaker); rsprfo → RS-P-RFO (Banerjee). All
                                  three Hessian TS optimizers
                                  (rsirfo/rsprfo/trim) are microiter-capable.
                                  [default: hess]
  --microiter / --no-microiter    Enable microiteration: alternate a 1-step
                                  macro TS move (RS-I-RFO / RS-P-RFO / TRIM) and
                                  MM relaxation (L-BFGS with MM-only forces).
                                  Effective in any Hessian --opt-mode
                                  (hess/rsirfo/rsprfo/trim); ignored in
                                  grad/dimer mode.  [default: microiter]
  --partial-hessian-flatten / --full-hessian-flatten
                                  Use partial (active-block) Hessian for
                                  imaginary mode detection in flatten loop.
                                  [default: partial-hessian-flatten]
  --flatten / --no-flatten        Enable/disable extra imaginary-mode flattening
                                  loop. --flatten uses the default
                                  flatten_max_iter (50); --no-flatten forces it
                                  to 0. When not provided, the loop is disabled
                                  unless YAML/config enables it.  [default: (no-
                                  flatten)]
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
                                  (default: uma).  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Enable the experimental, computationally
                                  expensive xTB point-charge delta correction
                                  for MLIP/MM.  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from the ML region for MM
                                  point charges used by the xTB delta
                                  correction.  [default: (12.0)]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
                                  [default: (scaled)]
  --mm-backend [hessian_ff|openmm]
                                  MM backend (default: hessian_ff). MM Hessians
                                  use finite differences by default; set
                                  calc.mm_fd: false for the hessian_ff
                                  analytical path.  [default: (hessian_ff)]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers. Default: enabled when present in
                                  parm7.  [default: (cmap)]
  --skip-final-freq / --no-skip-final-freq
                                  Skip the post-convergence frequency analysis
                                  and imaginary-mode flattening. Useful for
                                  large unfrozen systems where the final Hessian
                                  diagonalization is expensive.  [default: no-
                                  skip-final-freq]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) when explicit ML membership is
                                  absent. With explicit membership, retain valid
                                  movable/frozen MM B-factor layers.  [default:
                                  detect-layer]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --workers INTEGER               MLIP predictor workers (UMA). >1 uses a
                                  parallel predictor (fairchem-core[extras]);
                                  combining it with an analytical Hessian is an
                                  error. Default 1.  [default: (1)]
  --workers-per-node INTEGER      Workers per node when the parallel MLIP
                                  predictor is used (--workers > 1).  [default:
                                  (1)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace). Default: the backend's
                                  built-in model.  [default: (the selected
                                  backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator used as the ML-region backend
                                  (overrides --backend). Couples GFN-xTB / DFTB+
                                  / any ASE engine. See --calc-file-func-name.
  --calc-file-func-name TEXT      Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
                                  [default: (get_calculator)]
  --deterministic / --no-deterministic
                                  Request deterministic algorithms for
                                  controlled operations; verify exact
                                  reproducibility on the complete target stack.
                                  [default: no-deterministic]
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the default;
                                  command-specific choices are listed here.
                                  [default: (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). An open-shell ML region needs a
                                  matching multiplicity; use this only for an
                                  intentional nonstandard input such as a
                                  covalently-cut region.
  --stop-plateau / --no-stop-plateau
                                  Stop when the energy stops changing while the
                                  convergence criteria are still unmet, and
                                  report the run as stalled. It never signals
                                  convergence; --max-cycles remains the real
                                  bound. The MM micro iterations are never
                                  stopped this way.  [default: no-stop-plateau]
  --stop-plateau-thresh FLOAT     Energy range (hartree) below which --stop-
                                  plateau treats the window as flat.  [default:
                                  (1e-4)]
  --stop-plateau-window INTEGER   Number of consecutive cycles --stop-plateau
                                  inspects.  [default: (50)]
  -h, --help                      Show this message and exit.
```
