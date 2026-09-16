# `mlmm path-search`

```text
Usage: mlmm path-search [OPTIONS]

  Multistep MEP search via recursive GSM/DMF segmentation.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more PDB/mmCIF structures, or XYZ files
                                  with corresponding --ref-pdb entries, in
                                  reaction order. Either repeat '-i' (e.g., '-i
                                  A -i B -i C') or use a single '-i' followed by
                                  multiple space-separated paths (e.g., '-i A B
                                  C').  [required]
  --parm FILE                     Amber parm7 topology covering the full enzyme
                                  complex.  [required]
  --model-pdb FILE                ML-only, link-H-free PDB subset; atom
                                  identity/order must match the full PDB/parm7.
                                  When provided, it defines ML membership;
                                  --detect-layer still reads valid
                                  movable/frozen MM B-factors.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1).  [default: (1)]
  --mep-mode [gsm|dmf]            MEP method: gsm (Growing String) or dmf
                                  (Direct Max Flux).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  --refine-mode [peak|minima]     Refinement seed around the highest-energy
                                  image: 'peak' uses HEI±1, 'minima' uses
                                  nearest local minima. Defaults to peak for gsm
                                  and minima for dmf.  [default: (peak for gsm,
                                  minima for dmf)]
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.  [default: (use freeze_atoms)]
  --max-nodes INTEGER             Number of movable internal images per GSM or
                                  DMF segment (total images = max_nodes + 2
                                  endpoints); recursive segments may override it
                                  with YAML search.max_nodes_segment.  [default:
                                  20]
  --max-depth INTEGER RANGE       Number of recursive subdivision levels allowed
                                  while splitting a multistep path. 0 performs
                                  no subdivision, returning each input pair as
                                  one MEP segment (none when its HEI sits at an
                                  endpoint). Reaching the limit is not an error.
                                  Any segment retained at a positive cap is
                                  tagged seg_NNN_maxdepth and is not guaranteed
                                  to be a single elementary step. When not
                                  given, YAML search.max_depth applies.
                                  [default: (10); x>=0]
  --gsm-param [equi|energy]       GSM node parameterization after string growth.
                                  The energy scheme concentrates nodes in high-
                                  energy regions and may be tried when an
                                  equidistant path skips the reaction-coordinate
                                  region near the HEI.  [default: (equi)]
  --max-cycles-gsm INTEGER RANGE  Maximum GSM string-optimizer cycles for the
                                  MEP stage.  [default: (300); x>=1]
  --max-cycles-dmf INTEGER RANGE  Maximum IPOPT iterations for the DMF MEP
                                  stage. This is a solver iteration count, not a
                                  string-optimizer cycle count.  [default:
                                  (300); x>=1]
  --climb / --no-climb            Enable transition-state search after path
                                  growth.  [default: climb]
  --dump / --no-dump              Dump GSM/single-optimization trajectories
                                  during the run.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_search/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single L-BFGS runs
                                  only. The MEP itself keeps --thresh-gsm /
                                  --thresh-dmf.  [default: (gau)]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer (gau_loose|gau|gau_tight|gau_vtight|
                                  baker|never).  [default: (gau_loose)]
  --thresh-dmf TEXT               IPOPT dual-infeasibility tolerance for the DMF
                                  path optimizer: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float. This is not
                                  a Gaussian preset.  [default: (tight)]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path search.  [default: no-
                                  dry-run]
  --preopt / --no-preopt          If True, run initial single-structure
                                  optimizations of inputs.  [default: preopt]
  --align / --no-align            After optional preoptimization, align adjacent
                                  inputs in sequence and match frozen-atom
                                  positions while relaxing the remaining atoms.
                                  [default: align]
  --ref-pdb FILE                  Full-size template PDBs in the same reaction
                                  order as --input. Required when using XYZ
                                  inputs to provide topology and B-factor
                                  information.
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region.
                                  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Enable the experimental, computationally
                                  expensive xTB point-charge delta correction
                                  for MLIP/MM.  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from the ML region for MM
                                  point charges used by the xTB delta
                                  correction.  [default: (12.0)]
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
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) when explicit ML membership is
                                  absent. With explicit membership, retain valid
                                  movable/frozen MM B-factor layers.  [default:
                                  detect-layer]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based or
                                  0-based.  [default: model-indices-one-based]
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
                                  / off:small for mace).  [default: (the
                                  selected backend's own model)]
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
