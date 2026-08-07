# `mlmm irc`

```text
Usage: mlmm irc [OPTIONS]

  Run an IRC calculation with EulerPC. Only the documented CLI options are
  accepted; all other settings come from YAML.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .cif, .mmcif,
                                  .xyz, _trj.xyz, etc.).  [required]
  --parm FILE                     Amber parm7 topology for the whole enzyme (MM
                                  region). If omitted, must be provided in YAML
                                  as calc.real_parm7.
  --model-pdb FILE                ML-only, link-H-free PDB subset; atom
                                  identity/order must match the full PDB/parm7.
                                  When provided, it defines ML membership;
                                  --detect-layer still reads valid
                                  movable/frozen MM B-factors.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  -q, --charge INTEGER            Net charge of the ML region/model system;
                                  overrides calc.model_charge from YAML.
                                  Required unless --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
                                  omitted (requires PDB/mmCIF input or --ref-
                                  pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1); overrides
                                  calc.model_mult from YAML.  [default: (1)]
  --max-cycles INTEGER            Maximum number of IRC steps; overrides
                                  irc.max_cycles from YAML.  [default: (125)]
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates). Default: 0.10 Bohr. Overrides
                                  irc.step_length from YAML.  [default: (0.10)]
  --root INTEGER                  Imaginary mode index used for the initial
                                  displacement; overrides irc.root from YAML.
                                  [default: (0)]
  --forward / --no-forward        Run the forward IRC; overrides irc.forward
                                  from YAML.  [default: (forward)]
  --backward / --no-backward      Run the backward IRC; overrides irc.backward
                                  from YAML.  [default: (backward)]
  --never-stop / --no-never-stop  Ignore RMS-gradient, hard-gradient, energy-
                                  rise, and energy-change stops and trace until
                                  max-cycles. Numerical/integration failures and
                                  external interruption still stop the run;
                                  default off.  [default: (no-never-stop)]
  -o, --out-dir TEXT              Output directory; overrides irc.out_dir from
                                  YAML.  [default: ./result_irc/]
  --hessian-calc-mode [analytical|finitedifference]
                                  How the ML backend builds the Hessian
                                  (Analytical or FiniteDifference); overrides
                                  calc.hessian_calc_mode from YAML. Default:
                                  'FiniteDifference'. Runtime and memory depend
                                  on the backend and system; compare both modes
                                  on a representative pilot.  [default:
                                  (FiniteDifference)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running IRC.  [default: no-dry-run]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when
                                  --input is XYZ (keeps XYZ coordinates).
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Unavailable in v0.3.3; retained so older
                                  commands fail with an actionable diagnostic.
                                  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Unavailable in v0.3.3 together with the
                                  retired electronic-embedding path.  [default:
                                  (12.0)]
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
  --hess-device [auto|cuda|cpu]   Device for initial Hessian storage and IRC
                                  operations (auto/cuda/cpu). Use 'cpu' for
                                  large unfrozen systems to avoid VRAM limits.
                                  [default: auto]
  --read-hess FILE                Read an identified initial Hessian from 'mlmm
                                  freq --dump-hess'. Geometry, atom order,
                                  active-DOF basis, charge, and multiplicity
                                  must match; the file takes priority over
                                  hessian_cache and fresh computation.
  --allow-unverified-hess-state / --no-allow-unverified-hess-state
                                  Allow a schema-1 Hessian file whose charge and
                                  multiplicity cannot be verified. Use only
                                  after independently checking the electronic
                                  state.  [default: no-allow-unverified-hess-
                                  state]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
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
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
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
                                  [default: no-deterministic]
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). An open-shell ML region needs a
                                  matching multiplicity; use this only for an
                                  intentional nonstandard input such as a
                                  covalently-cut region.
  --irc-pos-def / --no-irc-pos-def
                                  Require pos-def Hessian at IRC convergence
                                  (blocks shoulder false-convergence).
  -h, --help                      Show this message and exit.
```
