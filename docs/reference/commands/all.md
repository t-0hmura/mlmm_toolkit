# `mlmm all`

```text
Usage: mlmm all [OPTIONS]

  Run pocket extraction → (optional scan-defined single-structure route) → MEP
  search in one command. If exactly one input is provided: (a) with --scan-
  lists, stage results feed into path-opt (or path_search with --refine-path);
  (b) with --tsopt and no --scan-lists, run TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more full PDB/mmCIF structures, or XYZ
                                  files sharing one matching --ref-pdb, in
                                  reaction order (reactant [intermediates ...]
                                  product); one full structure is allowed with
                                  --scan-lists or --tsopt. A single '-i' may be
                                  followed by multiple space-separated files
                                  (for example, '-i A.pdb B.pdb C.pdb').
                                  [required]
  -c, --center TEXT               Substrate specification for the extractor: a
                                  PDB path, a residue-ID list like '123,124' or
                                  'A:123,B:456' (insertion codes OK: '123A' /
                                  'A:123A'), or a residue-name list like
                                  'GPP,MMT'. When omitted, extraction is skipped
                                  and full structures are used directly.
  -o, --out-dir DIRECTORY         Top-level output directory for the pipeline.
                                  [default: result_all]
  -r, --radius FLOAT RANGE        Inclusion cutoff (Å) around substrate atoms.
                                  Zero is accepted and evaluated internally as
                                  0.001 Å (effectively off for ordinary radius-
                                  based neighbors).  [default: 2.6; x>=0.0]
  --radius-het2het FLOAT RANGE    Independent hetero–hetero cutoff (Å) for
                                  non‑C/H pairs.  [default: 0.0; x>=0.0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL)
                                  in the pocket.  [default: include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Remove backbone atoms on non‑substrate amino
                                  acids (with PRO/HYP safeguards).  [default:
                                  no-exclude-backbone]
  --add-linkh / --no-add-linkh    Add extractor-only link H to scratch pocket
                                  PDBs. The ML/MM model selection remains link-
                                  free; runtime link H are generated from parm7
                                  boundary bonds.  [default: no-add-linkh]
  --selected-resn TEXT            Force-include residues using IDs ('123',
                                  'A:123A'), names ('SAM'), or chain-qualified
                                  names ('A:SAM', 'A:SAM:123'); comma/space
                                  separated.  [default: ""]
  --modified-residue TEXT         Comma-separated modified-residue names with
                                  integer charges for backbone truncation and
                                  charge assignment. A known catalog residue may
                                  omit its charge. Example: 'HD1:0,SEP'.
                                  [default: ""]
  -l, --ligand-charge TEXT        Either a total charge (number) to distribute
                                  across unknown residues or a mapping like
                                  'GPP:-3,MMT:-1'.
  -q, --charge INTEGER            Override the net charge of the ML region/model
                                  atoms. Highest priority over the charge
                                  derived by the workflow.
  --parm FILE                     Pre-built AMBER parm7 topology file. When
                                  provided, mm_parm generation is skipped.
  --model-pdb FILE                ML-only atom-selection PDB. It must be an
                                  unchanged, link-H-free subset of the full
                                  PDB/parm7 in the same atom order. It takes
                                  precedence over ML membership from -c/--center
                                  or input B-factors.
  --auto-mm-ff-set [ff19sb|ff14sb]
                                  Force-field set forwarded to mm_parm (ff19SB
                                  uses OPC3; ff14SB uses TIP3P).  [default:
                                  ff19SB]
  --auto-mm-add-ter / --auto-mm-no-add-ter
                                  Control mm_parm TER insertion around
                                  ligand/water/ion blocks and disconnected
                                  peptide blocks.  [default: auto-mm-add-ter]
  --auto-mm-keep-temp             Keep the mm_parm temporary working directory
                                  (for debugging).
  --auto-mm-ligand-mult TEXT      Spin multiplicity mapping forwarded to mm_parm
                                  (e.g., 'GPP:2,SAM:1'). If omitted, mm_parm
                                  defaults to 1 for all ligands.  [default: (1)]
  --auto-mm-disulfide / --auto-mm-no-disulfide
                                  Forwarded to mm_parm: detect disulfides from
                                  SG-SG geometry across CYS/CYX and bond them
                                  (renaming a bonded CYS to CYX). With --auto-
                                  mm-no-disulfide only residues already named
                                  CYX are bonded and CYS is left untouched.
                                  [default: auto-mm-disulfide]
  -m, --multiplicity INTEGER      Multiplicity (2S+1).  [default: 1]
  --freeze-atoms TEXT             Comma-separated 1-based full-system atom
                                  indices to freeze throughout scan, MEP, TSOPT,
                                  endpoint optimization, IRC, and frequency
                                  stages (for example, '1,3,5'). Merged with
                                  YAML geom.freeze_atoms and the automatically
                                  detected Frozen-MM layer.
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  --max-nodes INTEGER             Max internal nodes per GSM/DMF segment
                                  (max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-depth INTEGER RANGE       Maximum recursive subdivision levels; requires
                                  --refine-path. 0 disables subdivision.
                                  Intervals retained at a positive cap use
                                  seg_NNN_maxdepth and may contain multiple
                                  steps.  [default: (10); x>=0]
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
                                  (3000); x>=1]
  --climb / --no-climb            Enable transition-state climbing after growth
                                  for the *first* segment in each pair.
                                  [default: climb]
  --opt-mode [grad|hess]          Fallback optimizer mode for TSOPT and post-IRC
                                  endpoint optimization: grad (=L-BFGS/Dimer) or
                                  hess (=RFO/RS-P-RFO). --opt-mode-post takes
                                  precedence.  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode for TSOPT and post-IRC endpoint
                                  optimizations. Takes precedence over --opt-
                                  mode for these stages.  [default: hess]
  --dump / --no-dump              Dump MEP / single-structure trajectories
                                  during the run, forwarding the same flag to
                                  scan/tsopt/freq.  [default: no-dump]
  --refine-path / --no-refine-path
                                  When disabled, run single-pass path-opt with
                                  the selected MEP optimizer between each
                                  adjacent pair and concatenate the segments (no
                                  path-search); when enabled, run recursive
                                  path-search on the full ordered series: it
                                  proposes multistep reaction paths and also
                                  refines a single-step MEP, which can improve a
                                  poor HEI or TS estimate.  [default: no-refine-
                                  path]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single-structure
                                  optimizations and scan relaxations (gau_loose|
                                  gau|gau_tight|gau_vtight|baker|never). The MEP
                                  stage keeps its own --thresh-gsm / --thresh-
                                  dmf.  [default: (gau)]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer of the MEP stage (gau_loose|gau|gau_
                                  tight|gau_vtight|baker|never).  [default:
                                  (gau_loose)]
  --thresh-dmf TEXT               IPOPT dual-infeasibility tolerance for the DMF
                                  MEP stage: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float. This is not
                                  a Gaussian preset.  [default: (tight)]
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for TS and post-IRC
                                  endpoint optimizations (gau_loose|gau|gau_tigh
                                  t|gau_vtight|baker|never).  [default: baker]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Run input preparation and preflight checks in
                                  a temporary directory, print the execution
                                  plan, and skip calculation stages.  [default:
                                  no-dry-run]
  --preopt / --no-preopt          Run initial single-structure optimizations of
                                  the pocket inputs.  [default: preopt]
  --hessian-calc-mode [analytical|finitedifference]
                                  Common MLIP Hessian mode forwarded to tsopt
                                  and freq. Runtime and memory depend on the
                                  backend and system; compare both modes on a
                                  representative pilot.  [default:
                                  (FiniteDifference)]
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) in downstream tools.  [default:
                                  detect-layer]
  --tsopt / --no-tsopt            TS optimization + EulerPC IRC per reactive
                                  segment (or TSOPT-only mode for single-
                                  structure), and build energy diagrams.
                                  [default: no-tsopt]
  --tsopt-from-mep-tan / --no-tsopt-from-mep-tan
                                  Guide Hessian-based TS root identity from MEP
                                  tangent candidate(s) at the highest-energy
                                  image. The CPU/file cache is not created or
                                  used when disabled. Dimer does not consume
                                  this Hessian reference mode.  [default: tsopt-
                                  from-mep-tan]
  --thermo / --no-thermo          Run freq on (R,TS,P) per reactive segment (or
                                  TSOPT-only mode) and build a Gibbs free-energy
                                  diagram (ML/MM).  [default: no-thermo]
  --dft / --no-dft                Run DFT single-point on (R,TS,P) and build a
                                  DFT energy diagram. With --thermo, also
                                  generate a DFT//MLIP/MM Gibbs diagram.
                                  [default: no-dft]
  --tsopt-max-cycles INTEGER RANGE
                                  Override tsopt --max-cycles.  [default:
                                  (100000); x>=1]
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt (grad: dimer loop, hess: post-
                                  RS-P-RFO); --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --reject-uphill / --no-reject-uphill
                                  Opt in to rejecting uphill RFO trials during
                                  post-IRC endpoint re-optimization only
                                  (tolerance: 1e-4 Hartree). Does not affect TS
                                  optimization or path search.  [default: no-
                                  reject-uphill]
  --stop-plateau / --no-stop-plateau
                                  Stop when the energy stops changing while the
                                  convergence criteria are still unmet, and
                                  report the run as stalled. It never signals
                                  convergence; --max-cycles remains the real
                                  bound. The MM micro iterations are never
                                  stopped this way.  [default: no-stop-plateau]
  --stop-plateau-thresh FLOAT     Energy range (hartree) below which --stop-
                                  plateau treats the window as flat.  [default:
                                  (1e-4)]
  --stop-plateau-window INTEGER   Number of consecutive cycles --stop-plateau
                                  inspects.  [default: (50)]
  --irc-step-size FLOAT           Override IRC --step-size (Bohr). If an IRC
                                  stops after only a few frames, retry with a
                                  smaller value such as 0.05.  [default: (0.10)]
  --irc-max-cycles INTEGER RANGE  Cycle cap for each post-TS IRC.  [default:
                                  (125); x>=1]
  --irc-never-stop / --no-irc-never-stop
                                  Forward IRC never-stop mode to every post-TS
                                  IRC. It ignores gradient and energy endpoint
                                  criteria and traces to the cycle cap;
                                  numerical/integration failures still stop.
                                  Default follows irc.never_stop (off).
                                  [default: (no-irc-never-stop)]
  --skip-final-freq / --no-skip-final-freq
                                  Skip terminal PHVA/frequency analysis in
                                  tsopt. The TS structure is retained with
                                  unverified saddle order, and all stops before
                                  IRC because no imaginary reaction direction
                                  can be validated.  [default: no-skip-final-
                                  freq]
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
                                  [default: (<segment>/ts)]
  --freq-out-dir DIRECTORY        Override freq output base directory (relative
                                  paths resolved against the default).
                                  [default: (<tsopt dir>/freq)]
  --freq-max-write INTEGER        Override freq --max-write value.  [default:
                                  (10)]
  --freq-amplitude-ang FLOAT      Override freq --amplitude-ang (Å).  [default:
                                  (0.8)]
  --freq-n-frames INTEGER         Override freq --n-frames value.  [default:
                                  (20)]
  --freq-sort [value|abs]         Override freq mode sorting.  [default:
                                  (value)]
  --freq-temperature FLOAT        Override freq thermochemistry temperature (K).
                                  [default: (298.15)]
  --freq-pressure FLOAT           Override freq thermochemistry pressure (atm).
                                  [default: (1.0)]
  --dft-out-dir DIRECTORY         Override dft output base directory (relative
                                  paths resolved against the default).
                                  [default: (<tsopt dir>/dft)]
  --dft-func-basis TEXT           Override dft --func-basis value.  [default:
                                  (wb97m-v/def2-tzvpd)]
  --dft-max-cycle INTEGER RANGE   Override dft --max-cycle value.  [default:
                                  (100); x>=1]
  --dft-conv-tol FLOAT            Override dft --conv-tol value.  [default:
                                  (1e-09)]
  --dft-grid-level INTEGER        Override dft --grid-level value.  [default:
                                  (3)]
  --dft-engine [gpu|cpu]          Override dft --engine value.  [default: (gpu)]
  -s, --scan-lists TEXT           Scan targets: inline Python literals
                                  containing (i,j,target) triples. Multiple
                                  literals define sequential stages, e.g.
                                  "[(12,45,1.35)]"
                                  "[(10,55,2.20),(23,34,1.80)]". Use standalone
                                  mlmm scan for YAML/JSON or bidirectional
                                  4-tuples. Indices refer to the original full
                                  PDB (1-based) or PDB atom selectors like
                                  "TYR,285,CA".
  --scan-out-dir DIRECTORY        Override the scan output directory (default:
                                  <out-dir>/_work/scan). Relative paths are
                                  resolved against the default parent.
                                  [default: (<out-dir>/_work/scan)]
  --scan-one-based / --scan-zero-based
                                  Override scan indexing interpretation (one-
                                  based or zero-based).  [default: (True (one-
                                  based))]
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å).  [default:
                                  (0.2)]
  --scan-bias-k FLOAT             Override scan harmonic bias strength k
                                  (eV/Å^2).  [default: (300.0)]
  --scan-relax-max-cycles INTEGER RANGE
                                  Override scan relaxation max cycles per step.
                                  [default: (100000); x>=1]
  --scan-preopt / --no-scan-preopt
                                  Override scan --preopt flag.  [default:
                                  (inherits --preopt)]
  --scan-endopt / --no-scan-endopt
                                  Override scan --endopt flag.  [default:
                                  (inherits --endopt)]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs to PDB format using
                                  reference topology; forwarded to all
                                  subcommands.  [default: convert-files]
  --ref-pdb FILE                  Reference PDB for topology/B-factor layer
                                  information when -i provides XYZ inputs. Used
                                  for define-layer, mm_parm, ml_region, and
                                  forwarded to downstream tools (tsopt, irc,
                                  freq, path_search) as --ref-pdb.
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region.
                                  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Enable experimental point-charge treatment.
                                  MLIP/MM stages use the computationally
                                  expensive xTB delta correction; DFT/MM stages
                                  embed MM charges in the PySCF Hamiltonian.
                                  [default: no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from the ML region for
                                  embedded MM point charges.  [default: (12.0)]
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
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  cart is the default; command-specific choices
                                  are listed here.  [default: (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
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
