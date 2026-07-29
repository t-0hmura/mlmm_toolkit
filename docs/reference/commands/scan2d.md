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
  -i, --input FILE                Input enzyme complex PDB.  [required]
  --parm FILE                     Amber parm7 topology for the enzyme.
                                  [required]
  --model-pdb FILE                ML-only, link-H-free PDB subset; atom
                                  identity/order must match the full PDB/parm7.
                                  When provided, it defines ML membership;
                                  --detect-layer still reads valid
                                  movable/frozen MM B-factors.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  -q, --charge INTEGER            ML-region total charge. Required unless
                                  --ligand-charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
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
                                  YAML/JSON spec file path.  [required]
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  --scan-lists.  [default: no-print-parsed]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running the scan.  [default: no-dry-
                                  run]
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
  --relax-max-cycles INTEGER      Maximum L-BFGS cycles per biased relaxation
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
                                  / MACE-OFF23_small for mace). Default: the
                                  backend's built-in model.
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator used as the ML-region backend
                                  (overrides --backend). Couples GFN-xTB / DFTB+
                                  / any ASE engine. See --calc-factory.
  --calc-factory TEXT             Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance). CLI overrides config
                                  YAML; otherwise defaults to get_calculator.
  --deterministic / --no-deterministic
                                  Enable strict deterministic GPU algorithms and
                                  the index_reduce_ shim. Slower; raises if
                                  unsupported. Default off.
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). For an intentional open-shell or
                                  covalently-cut ML region.
  -h, --help                      Show this message and exit.
```
