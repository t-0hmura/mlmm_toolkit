# `mlmm opt`

```text
Usage: mlmm opt [OPTIONS]

  ML/MM geometry optimization with L-BFGS (light) or RFO (heavy).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (PDB, XYZ). XYZ provides
                                  higher coordinate precision. If XYZ, use
                                  --ref-pdb to specify PDB topology for atom
                                  ordering and output conversion.  [required]
  --ref-pdb FILE                  Reference PDB topology when input is XYZ. XYZ
                                  coordinates are used (higher precision) while
                                  PDB provides atom ordering and residue
                                  information for output conversion.
  --parm FILE                     Amber parm7 topology covering the whole enzyme
                                  complex.  [required]
  --model-pdb FILE                ML-only, link-H-free PDB subset; atom
                                  identity/order must match the full PDB/parm7.
                                  When provided, it defines ML membership;
                                  --detect-layer still reads valid
                                  movable/frozen MM B-factors.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --radius-partial-hessian, --hess-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer. `--hess-cutoff`
                                  is a compatibility alias.
  --radius-freeze, --movable-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --radius-freeze disables --detect-
                                  layer and uses distance-based layer
                                  assignment. `--movable-cutoff` is a
                                  compatibility alias.
  --dist-freeze TEXT              Distance restraints: inline Python literal
                                  (e.g. '[(1,5,1.4)]') or a YAML/JSON spec file
                                  path. Format: (i,j,target_Å) triples. Target
                                  may be omitted to freeze at the current
                                  distance: (i,j).
  --one-based / --zero-based      Interpret --dist-freeze indices as 1-based
                                  (default) or 0-based.  [default: one-based]
  --bias-k FLOAT                  Harmonic restraint strength k [eV/Å^2] for
                                  --dist-freeze. Defaults to BIAS_KW['k']=300
                                  (in defaults.py) when omitted.
  --max-cycles INTEGER            Maximum number of optimization cycles.
                                  [default: 10000]
  --dump / --no-dump              Write optimization trajectories
                                  ('optimization_trj.xyz' and
                                  'optimization_all_trj.xyz').  [default: no-
                                  dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_opt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset.
  --opt-mode [grad|hess|light|heavy|lbfgs|rfo]
                                  Optimization mode: grad (lbfgs) or hess (rfo).
                                  Aliases light/heavy and lbfgs/rfo are
                                  accepted.  [default: grad]
  --microiter / --no-microiter    Enable microiteration: alternate ML 1-step
                                  (RFO) and MM relaxation (L-BFGS with MM-only
                                  forces). Only effective in --opt-mode hess
                                  (RFO). Ignored in grad mode.  [default:
                                  microiter]
  --flatten / --no-flatten        Enable/disable imaginary-mode flatten loop
                                  after optimization.  [default: no-flatten]
  --reject-uphill / --no-reject-uphill
                                  Opt in to rejecting uphill RFO trials in hess
                                  mode (tolerance: 1e-4 Hartree) and final-check
                                  the retained geometry at the emergency floor.
                                  Ignored in grad/lbfgs mode.  [default: no-
                                  reject-uphill]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running optimization.  [default: no-
                                  dry-run]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Unavailable in v0.3.3; retained so older
                                  commands fail with an actionable diagnostic.
                                  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Unavailable in v0.3.3 together with the
                                  retired electronic-embedding path.
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  --mm-backend [hessian_ff|openmm]
                                  MM backend (default: hessian_ff). MM Hessians
                                  use finite differences by default; set
                                  calc.mm_fd: false for the hessian_ff
                                  analytical path.
  --mm-only / --no-mm-only        Skip the MLIP component entirely and minimize
                                  using only the MM force field on the full
                                  system. Layers (movable/frozen) are still
                                  honored via B-factor encoding or --radius-
                                  freeze. Only --opt-mode grad (L-BFGS) is
                                  supported in this mode; microiteration is
                                  automatically disabled.  [default: no-mm-only]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers. Default: enabled when present in
                                  parm7.
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --detect-layer                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) when explicit ML membership is
                                  absent. With explicit membership, retain valid
                                  movable/frozen MM B-factor layers.  [default:
                                  True]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.  [x>=1]
  --coord-type [cart|redund|dlc|tric]
                                  Optimization coordinate system
                                  (cart|redund|dlc|tric). cart is the default;
                                  command-specific choices are listed here.
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [x>=1]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --workers INTEGER               MLIP predictor workers (UMA). >1 uses a
                                  parallel predictor (fairchem-core[extras]);
                                  combining it with an analytical Hessian is an
                                  error. Default 1.
  --workers-per-node INTEGER      Workers per node when the parallel MLIP
                                  predictor is used (--workers > 1).
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace). Default: the backend's
                                  built-in model.
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator used as the ML-region backend
                                  (overrides --backend). Couples GFN-xTB / DFTB+
                                  / any ASE engine. See --calc-factory.
  --calc-factory TEXT             Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
  --deterministic / --no-deterministic
                                  Request deterministic algorithms for
                                  controlled operations; verify exact
                                  reproducibility on the complete target stack.
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). For an intentional open-shell or
                                  covalently-cut ML region.
  -h, --help                      Show this message and exit.
```
