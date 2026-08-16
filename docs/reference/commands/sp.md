# `mlmm sp`

```text
Usage: mlmm sp [OPTIONS]

  Compute a single-point ML/MM ONIOM energy + forces (and optionally Hessian).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Layered PDB/mmCIF, or XYZ with --ref-pdb,
                                  defining the ML/MM/Frozen system.  [required]
  --ref-pdb FILE                  Full-system PDB/mmCIF topology required when
                                  --input is XYZ.
  --parm, --real-parm7 FILE       Amber parm7 of the full enzyme (canonical flag
                                  is --parm; --real-parm7 retained as alias).
                                  [required]
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
  --hess-cutoff FLOAT             Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms; combinable with
                                  --detect-layer.  [default: (all movable MM
                                  atoms)]
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  [default: (use freeze_atoms)]
  -q, --charge INTEGER            ML region total charge.
  -l, --ligand-charge TEXT        Per-ligand charge mapping, e.g.
                                  'SAM:1,GPP:-3'.
  -m, --multiplicity INTEGER      ML region spin multiplicity (2S+1).  [default:
                                  (1)]
  -o, --out-dir TEXT              Output directory.  [default: ./result_sp/]
  --hess / --no-hess              Also compute the active-coordinate ONIOM
                                  Hessian and save to hessian.npy.  [default:
                                  no-hess]
  --hessian-calc-mode [analytical|finitedifference]
                                  Hessian backend when --hess is set. Analytical
                                  is supported by UMA, ORB, MACE, and AIMNet2;
                                  custom calculators use FiniteDifference.
                                  Analytical cannot be combined with --workers >
                                  1.  [default: (FiniteDifference)]
  --config FILE                   YAML config file with sections (calc:, geom:,
                                  …).
  --show-config / --no-show-config
                                  Print effective merged config and exit.
                                  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the plan without
                                  running.  [default: no-dry-run]
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
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
                                  Link-atom positioning: scaled (g-factor) or
                                  fixed (1.09/1.01 Å).  [default: (scaled)]
  --mm-backend [hessian_ff|openmm]
                                  MM backend (default: hessian_ff).  [default:
                                  (hessian_ff)]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers. Default: enabled when present in
                                  parm7.  [default: (cmap)]
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
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). An open-shell ML region needs a
                                  matching multiplicity; use this only for an
                                  intentional nonstandard input such as a
                                  covalently-cut region.
  -h, --help                      Show this message and exit.
```
