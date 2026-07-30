# `mlmm freq`

```text
Usage: mlmm freq [OPTIONS]

  ML/MM vibrational frequency analysis (PHVA-compatible).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Enzyme complex PDB used by both geom_loader
                                  and the ML/MM calculator.  [required]
  --parm FILE                     Amber parm7 topology for the full enzyme
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
  --tr-projection [constrained|legacy-active]
                                  Rigid-mode treatment for PHVA. 'constrained'
                                  removes only full-system rigid motions
                                  compatible with frozen anchors (default);
                                  'legacy-active' is deprecated comparison-only
                                  behavior and must not be used for pass/HOSP
                                  transition-state certification.
  --hess-cutoff FLOAT             Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer.
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.
  --hessian-calc-mode [analytical|finitedifference]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); overrides
                                  calc.hessian_calc_mode from YAML. Default:
                                  'FiniteDifference'. Runtime and memory depend
                                  on the backend and system; compare both modes
                                  on a representative pilot.
  --max-write INTEGER             Maximum number of modes to export.  [default:
                                  10]
  --amplitude-ang FLOAT           Mode animation amplitude (Å).  [default: 0.8]
  --n-frames INTEGER              Frames per vibrational mode animation.
                                  [default: 20]
  --sort [value|abs]              Sort modes by signed value or absolute value.
                                  [default: value]
  --temperature FLOAT             Temperature (K) for thermochemistry summary.
                                  [default: 298.15]
  --pressure FLOAT                Pressure (atm) for thermochemistry summary.
                                  [default: 1.0]
  --symmetry-number INTEGER RANGE
                                  External rotational symmetry number used in
                                  the thermochemistry partition function.
                                  [default: 1; x>=1]
  --dump / --no-dump              Write 'thermoanalysis.yaml' alongside the
                                  console summary.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default: ./result_freq/]
  --active-dof-mode [all|ml-only|partial|unfrozen]
                                  Active DOF selection for frequency analysis:
                                  all (all atoms), ml-only (ML only), partial
                                  (ML + MovableMM, default), unfrozen (all non-
                                  frozen atoms).  [default: partial]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running frequency analysis.  [default:
                                  no-dry-run]
  --ref-pdb FILE                  Reference PDB topology to use when --input is
                                  XYZ (keeps XYZ coordinates).
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  --hess-device [auto|cuda|cpu]   Device for Hessian assembly and
                                  diagonalization (auto/cuda/cpu). Use 'cpu' to
                                  avoid VRAM issues with large unfrozen systems.
                                  ML model inference always uses ml_device
                                  (typically GPU).  [default: auto]
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
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers. Default: enabled when present in
                                  parm7.
  --dump-hess FILE                Save the computed Hessian and geometry/active-
                                  basis identity to a compressed .npz file for a
                                  matching 'mlmm irc --read-hess' run. The file
                                  also identifies model charge and multiplicity.
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --detect-layer / --no-detect-layer
                                  Without --model-pdb/--model-indices, detect
                                  ML/MM layers from input PDB B-factors (ML=0,
                                  MovableMM=10, FrozenMM=20). With explicit ML
                                  membership, retain valid movable/frozen MM
                                  B-factor layers. If disabled, explicit
                                  membership is required.  [default: detect-
                                  layer]
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
