# `mlmm scan`

```text
Usage: mlmm scan [OPTIONS]

  Bond-length driven scan with staged harmonic restraints and relaxation
  (ML/MM).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Full-system PDB/mmCIF, or XYZ with --ref-pdb,
                                  used by the ML/MM calculator.  [required]
  --parm FILE                     Amber parm7 topology covering the entire
                                  enzyme complex.  [required]
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
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer.  [default: (all
                                  movable MM atoms)]
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.  [default: (use freeze_atoms)]
  -s, --scan-lists TEXT           Scan targets: inline Python literal (e.g.
                                  '[(1,5,1.4)]') or a YAML/JSON spec file path.
                                  Multiple inline literals define sequential
                                  stages.
  --one-based / --zero-based      Interpret (i,j) indices in --scan-lists as
                                  1-based (default) or 0-based.  [default: one-
                                  based]
  --print-parsed / --no-print-parsed
                                  Print parsed scan targets after resolving
                                  -s/--scan-lists.  [default: no-print-parsed]
  --max-step-size FLOAT           Maximum change in any scanned bond length per
                                  step [Å].  [default: 0.2]
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2]. Defaults to
                                  YAML bias.k (BIAS_KW['k']=300 in defaults.py)
                                  when omitted; explicit CLI value overrides
                                  YAML.  [default: (300.0)]
  --opt-mode [grad|hess|lbfgs|rfo|light|heavy]
                                  Compatibility option for mlmm all forwarding.
                                  Scan relaxations always use L-BFGS; values
                                  other than grad/lbfgs/light emit a warning.
  --max-cycles INTEGER            Maximum L-BFGS cycles per biased step and per
                                  (pre|end)opt stage.  [default: 10000]
  --relax-max-cycles INTEGER      Compatibility alias of --max-cycles (overrides
                                  it when provided).
  --dump / --no-dump              Write per-step optimizer trajectory files.
                                  scan_trj.xyz is always written per-stage and
                                  as a combined file in out-dir; scan.pdb
                                  companions are written when --convert-files is
                                  enabled.  [default: no-dump]
  -o, --out-dir TEXT              Base output directory.  [default:
                                  ./result_scan/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for relaxations.  [default:
                                  (gau)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --ref-pdb FILE                  Reference PDB topology to use when --input is
                                  XYZ (keeps XYZ coordinates).
  --preopt / --no-preopt          Pre-optimize initial structure without bias
                                  before the scan.  [default: no-preopt]
  --endopt / --no-endopt          After each stage, run an additional unbiased
                                  optimization of the stage result.  [default:
                                  no-endopt]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running the scan.  [default: no-dry-
                                  run]
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
                                  Defaults to 1 when omitted.  [default: (1);
                                  x>=1]
  --coord-type [cart|redund|dlc|tric]
                                  Compatibility input for composite workflows.
                                  ML/MM restrained scan relaxation always uses
                                  Cartesian coordinates; non-cart values are
                                  accepted with a notice and resolved to cart.
  --print-every INTEGER RANGE     Print optimizer status every N cycles (debug
                                  knob).  [default: (100); x>=1]
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
                                  / any ASE engine. See --calc-file-func-name.
  --calc-file-func-name, --calc-factory TEXT
                                  Name of the callable in --calc-file that
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
  -h, --help                      Show this message and exit.
```
