# `mlmm scan`

```text
Usage: mlmm scan [OPTIONS]

  Internal-coordinate scan with harmonic restraints and relaxation (ML/MM).

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
  --parm7, --parm FILE            Amber parm7 topology covering the entire
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
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.  [default: (use freeze_atoms)]
  -s, --scan-lists TEXT           Distance targets (i,j,target), or scan ranges:
                                  distance (i,j,low,high), angle
                                  (i,j,k,low,high), or dihedral
                                  (i,j,k,l,low,high). Multiple inline literals
                                  define sequential stages.
  --one-based / --zero-based      Interpret atom indices in --scan-lists as
                                  1-based or 0-based.  [default: one-based]
  --max-step-size FLOAT           Maximum scanned distance change per step [Å].
                                  [default: 0.2]
  --max-angle-step-size FLOAT RANGE
                                  Maximum scanned angle change per step
                                  [degree].  [default: 5.0; x>0.0]
  --max-dihedral-step-size FLOAT RANGE
                                  Maximum scanned dihedral change per step
                                  [degree].  [default: 10.0; x>0.0]
  --restraint-k, --bias-k FLOAT   Harmonic well strength k [eV/Å^2 for
                                  distances; eV/rad^2 for angles]. YAML bias.k
                                  applies when this option is omitted; explicit
                                  CLI wins.  [default: (300.0)]
  --relax-max-cycles INTEGER RANGE
                                  Maximum optimizer cycles per biased step and
                                  per (pre|end)opt stage.  [default: (100000);
                                  x>=1]
  --opt-mode [grad|hess]          Relaxation mode: grad (=LBFGS) or hess (=RFO).
                                  [default: grad]
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
  --dry-run / --no-dry-run        Resolve and validate options (input,
                                  charge/spin, --scan-lists parse) and print the
                                  planned scan, then exit without running any
                                  optimization.  [default: no-dry-run]
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
                                  fixed (1.09/1.01 Å).  [default: (scaled)]
  --mm-backend [hessian_ff|openmm]
                                  MM backend. MM Hessians use finite differences
                                  by default; set calc.mm_fd: false for the
                                  hessian_ff analytical path.  [default:
                                  (hessian_ff)]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers when present in parm7.  [default:
                                  (cmap)]
  --out-json / --no-out-json      Write result.json to the output directory.
                                  [default: no-out-json]
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) when explicit ML membership is
                                  absent. With explicit membership, retain valid
                                  movable/frozen MM B-factor layers.  [default:
                                  detect-layer]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -l, --ligand-charge TEXT        Total charge for unknown ligand residues or a
                                  per-resname mapping (e.g., GPP:-3,SAM:1), used
                                  to derive the ML-region charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER RANGE
                                  Spin multiplicity (2S+1) for the ML region.
                                  [default: (1); x>=1]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
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
  --scf-stepwise-grid / --no-scf-stepwise-grid
                                  Converge the first SCF on a coarse grid, then
                                  on the final grid (later SCFs reuse the
                                  previous density as usual).  [default:
                                  (enabled)]
  --dft-nprocs INTEGER RANGE      PySCF/OpenMP CPU threads; GPU count is
                                  unaffected.  [default: (auto); x>=1]
  --dft-memory, --dft-mem TEXT    PySCF host RAM limit (for example 64GB or
                                  120000MB).  [default: (auto)]
  -h, --help                      Show this message and exit.
```
