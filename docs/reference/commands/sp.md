# `mlmm sp`

```text
Usage: mlmm sp [OPTIONS]

  Compute a single-point ML/MM ONIOM energy + forces (and optionally Hessian).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Layered PDB (or XYZ) defining the ML/MM/Frozen
                                  system.  [required]
  --parm, --real-parm7 FILE       Amber parm7 of the full enzyme (canonical flag
                                  is --parm; --real-parm7 retained as alias).
                                  [required]
  --model-pdb FILE                PDB defining atoms that belong to the ML
                                  (high-level) region. Optional when --detect-
                                  layer is enabled.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --radius-partial-hessian, --hess-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms; combinable with
                                  --detect-layer.
  --radius-freeze, --movable-cutoff FLOAT
                                  Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
  -q, --charge INTEGER            ML region total charge.
  -l, --ligand-charge TEXT        Per-ligand charge mapping, e.g.
                                  'SAM:1,GPP:-3'.
  -m, --multiplicity INTEGER      ML region spin multiplicity (2S+1).
  -o, --out-dir TEXT              Output directory.  [default: ./result_sp/]
  --hess / --no-hess              Also compute the full ONIOM Hessian and save
                                  to hessian.npy.  [default: no-hess]
  --hessian-calc-mode [analytical|finitedifference]
                                  Hessian backend when --hess is set. Analytical
                                  only works for UMA; other backends fall back
                                  to FiniteDifference.
  --convert-files / --no-convert-files
                                  Auto-convert output XYZ-like files into
                                  matching PDB beside them.  [default: convert-
                                  files]
  --config FILE                   YAML config file with sections (calc:, geom:,
                                  …).
  --show-config / --no-show-config
                                  Print effective merged config and exit.
  --dry-run / --no-dry-run        Validate options and print the plan without
                                  running.
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM->ML environmental effects.  [default:
                                  no-embedcharge]
  --embedcharge-cutoff FLOAT      Embed-charge cutoff radius (Å) around the ML
                                  region.
  --link-atom-method [scaled|fixed]
                                  Link-atom positioning: scaled (g-factor) or
                                  fixed (1.09/1.01 Å).
  --mm-backend [hessian_ff|openmm]
                                  MM backend (default: hessian_ff).
  --cmap / --no-cmap              Enable CMAP (backbone cross-map) terms in
                                  model parm7.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  --precision [fp32|fp64]         MLIP backend precision: fp32 (default) or
                                  fp64. Routed to backend-specific kwargs (UMA
                                  precision / ORB precision / MACE
                                  default_dtype). aimnet2: fp32 no-op; fp64
                                  rejected.
  --deterministic / --no-deterministic
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [x>=1]
  -h, --help                      Show this message and exit.
```
