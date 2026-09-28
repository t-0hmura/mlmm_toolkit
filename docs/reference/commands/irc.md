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
  --parm7, --parm FILE            Amber parm7 topology for the whole enzyme (MM
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
  --max-cycles INTEGER RANGE      Maximum number of IRC steps.  [default: (125);
                                  x>=1]
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates); overrides irc.step_length from
                                  YAML.  [default: (0.10)]
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
                                  calc.hessian_calc_mode from YAML. Runtime and
                                  memory depend on the backend and system;
                                  compare both modes on a representative pilot.
                                  [default: (FiniteDifference)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print the loaded YAML file and its top-level
                                  keys, then continue.  [default: no-show-
                                  config]
  --dry-run / --no-dry-run        Validate options and inputs without running
                                  IRC.  [default: no-dry-run]
  --ref-pdb FILE                  Reference PDB/mmCIF topology to use when
                                  --input is XYZ (keeps XYZ coordinates).
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2|dft]
                                  High-level backend for the ONIOM model region.
                                  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Enable electrostatic embedding. MLIP backends
                                  use the computationally expensive xTB point-
                                  charge delta correction; dft uses native PySCF
                                  MM point charges.  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from the ML region for MM
                                  point charges used by embedding.  [default:
                                  (12.0)]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor) or
                                  fixed (legacy 1.09/1.01 Å).  [default:
                                  (scaled)]
  --mm-backend [hessian_ff|openmm]
                                  MM backend. MM Hessians use finite differences
                                  by default; set calc.mm_fd: false for the
                                  hessian_ff analytical path.  [default:
                                  (hessian_ff)]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers when present in parm7.  [default:
                                  (cmap)]
  --hess-device [auto|cuda|cpu]   Device for initial Hessian storage and IRC
                                  operations (auto/cuda/cpu). Use 'cpu' for
                                  large unfrozen systems to avoid VRAM limits.
                                  [default: auto]
  --read-hess FILE                Start from the Hessian in this .npy file (e.g.
                                  from freq or tsopt --dump-hess): the Cartesian
                                  Hessian of the input geometry in
                                  Hartree/bohr^2, for all atoms or only the
                                  atoms in the Hessian calculation. The file
                                  takes priority over a Hessian from an earlier
                                  stage and fresh computation.  [default:
                                  (None)]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) when explicit ML membership is
                                  absent. With explicit membership, retain valid
                                  movable/frozen MM B-factor layers.  [default:
                                  detect-layer]
  --precision [fp32|fp64]         MLIP backend precision: fp32 or fp64. Unset
                                  defaults per backend (uma: fp32; orb, mace:
                                  fp64). Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.  [default: (per backend: uma fp32;
                                  orb, mace fp64)]
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers (UMA). >1 uses a
                                  parallel predictor (fairchem-core[extras]);
                                  combining it with an analytical Hessian is an
                                  error. Default 1.  [default: (1)]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when the parallel MLIP
                                  predictor is used (--workers > 1).  [default:
                                  (1)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / off:small for mace).  [default: (the
                                  selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator used as the ML-region backend
                                  (overrides --backend). Couples GFN-xTB / DFTB+
                                  / any ASE engine. See --calc-file-func-name.
  --calc-factory, --calc-file-func-name TEXT
                                  Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
                                  [default: (get_calculator)]
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
                                  [default: (no-irc-pos-def (rms-only
                                  criterion))]
  --func-basis TEXT               High-level method as FUNCTIONAL/BASIS;
                                  HF/BASIS is accepted.  [default:
                                  (wb97m-v/def2-svp)]
  --dft-engine, --engine [gpu|cpu]
                                  PySCF execution engine used by --backend dft.
                                  [default: (gpu)]
  --save-scf-checkpoint / --no-save-scf-checkpoint
                                  Persist a structure-bound PySCF checkpoint.
                                  [default: (disabled)]
  --scf-checkpoint FILE           Load/save the optional structure-bound PySCF
                                  checkpoint at PATH.
  --dft-low-memory, --lowmem / --no-dft-low-memory, --no-lowmem
                                  Use GPU4PySCF rks_lowmem for closed-shell GPU
                                  DFT; open-shell GPU and CPU use standard
                                  direct JK. --no-lowmem enables density
                                  fitting.  [default: (lowmem)]
  --dft-nprocs INTEGER RANGE      PySCF/OpenMP CPU threads; GPU count is
                                  unaffected.  [default: (auto); x>=1]
  --dft-memory, --dft-mem TEXT    PySCF host RAM limit (for example 64GB or
                                  120000MB).  [default: (auto)]
  -h, --help                      Show this message and exit.
```
