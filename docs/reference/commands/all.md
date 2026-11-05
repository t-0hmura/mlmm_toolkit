# `mlmm all`

```text
Usage: mlmm all [OPTIONS]

  Run ML-region selection → optional staged scan → MEP search in one command. If
  exactly one input is provided: (a) with --scan-lists, run staged scan and use
  stage results as inputs for the MEP search; (b) with --tsopt and no --scan-
  lists, run TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, output paths;
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
  -c, --center TEXT               Extraction centers (normally substrate +
                                  catalytic residues): a PDB path, a residue-ID
                                  list like '123,124' or 'A:123,B:456'
                                  (insertion codes OK: '123A' / 'A:123A'), or a
                                  residue-name list like 'GPP,MMT'. Each match
                                  starts radius expansion. When omitted,
                                  extraction is skipped and full structures are
                                  used directly.
  -o, --out-dir DIRECTORY         Top-level output directory.  [default:
                                  result_all]
  -r, --radius FLOAT RANGE        Inclusion cutoff (Å) around center atoms. 0
                                  disables radius-based expansion.  [default:
                                  2.6; x>=0.0]
  --radius-het2het FLOAT RANGE    Independent hetero-hetero cutoff (Å) for
                                  non-C/H pairs.  [default: 0.0; x>=0.0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL)
                                  in the pocket.  [default: include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Exclude amino-acid backbone atoms from the ML
                                  region, except between peptide-bonded
                                  extraction centers.  [default: no-exclude-
                                  backbone]
  --add-linkh / --no-add-linkh    Add extractor-only link H to scratch pocket
                                  PDBs. The ML/MM model selection remains link-
                                  free; runtime link H are generated from parm7
                                  boundary bonds.  [default: no-add-linkh]
  --selected-resn TEXT            Force-include residues without radius
                                  expansion using IDs ('123', 'A:123A'), names
                                  ('SAM'), or chain-qualified names ('A:SAM',
                                  'A:SAM:123'); comma/space separated.
                                  [default: ""]
  --modified-residue TEXT         Comma-separated modified-residue names with
                                  integer charges for backbone truncation and
                                  charge assignment. A known catalog residue may
                                  omit its charge. Example: 'HD1:0,SEP'.
                                  [default: ""]
  -l, --ligand-charge TEXT        Either a total charge (number) to distribute
                                  across unknown residues or a mapping like
                                  'GPP:-3,MMT:-1'.
  -q, --charge INTEGER            Net charge of the ML region. Overrides the
                                  charge derived by the workflow.
  --parm7, --parm FILE            Pre-built AMBER parm7 topology file. When
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
  --auto-mm-add-ter / --no-auto-mm-add-ter
                                  Control mm_parm TER insertion around
                                  ligand/water/ion blocks and disconnected
                                  peptide blocks.  [default: auto-mm-add-ter]
  --auto-mm-keep-temp             Keep the mm_parm temporary working directory
                                  (for debugging).
  --auto-mm-ligand-mult TEXT      Spin multiplicity mapping forwarded to mm_parm
                                  (e.g., 'GPP:2,SAM:1'). If omitted, mm_parm
                                  defaults to 1 for all ligands.  [default: (1)]
  --auto-mm-disulfide / --no-auto-mm-disulfide
                                  Forwarded to mm_parm: detect disulfides from
                                  SG-SG geometry across CYS/CYX and bond them
                                  (renaming a bonded CYS to CYX). With --no-
                                  auto-mm-disulfide only residues already named
                                  CYX are bonded and CYS is left untouched.
                                  [default: auto-mm-disulfide]
  -m, --multiplicity INTEGER RANGE
                                  Multiplicity (2S+1).  [default: 1; x>=1]
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
                                  --refine-path. 0 disables subdivision. With a
                                  positive limit, an interval still undivided at
                                  the limit is named seg_NNN_maxdepth and may
                                  contain several steps.  [default: (10); x>=0]
  --gsm-param [equi|energy]       GSM node parameterization after string growth.
                                  The energy scheme concentrates nodes in high-
                                  energy regions and may be tried when an
                                  equidistant path skips the reaction-coordinate
                                  region near the HEI.  [default: (equi)]
  --max-cycles-gsm INTEGER RANGE  Maximum GSM string-optimizer cycles for the
                                  MEP stage.  [default: (300); x>=1]
  --dmf-max-iterations, --max-cycles-dmf INTEGER RANGE
                                  Maximum IPOPT iterations for the DMF MEP
                                  stage.  [default: (3000); x>=1]
  --climb / --no-climb            Enable climbing image for standard GSM
                                  segments (bridge segments always disable
                                  climbing).  [default: climb]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/tsopt and
                                  used for single optimizations: grad
                                  (=LBFGS/Dimer) or hess (=RFO for scan/opt; RS-
                                  P-RFO for tsopt).  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode for TSOPT and post-IRC endpoint
                                  optimizations. Takes precedence over --opt-
                                  mode for these stages.  [default: hess]
  --dump / --no-dump              Dump MEP / single-structure trajectories
                                  during the run, forwarding the same flag to
                                  scan/tsopt/freq.  [default: no-dump]
  --refine-path / --no-refine-path
                                  By default, run a single-pass path-opt with
                                  the selected GSM/DMF method between each
                                  adjacent pair and concatenate the segments.
                                  Use --refine-path to run recursive path-search
                                  on the full ordered series: it proposes
                                  multistep reaction paths and also refines a
                                  single-step MEP, which can improve a poor HEI
                                  or TS estimate.  [default: no-refine-path]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for single-structure
                                  optimizations and scan relaxations. The MEP
                                  stage keeps its own --thresh-gsm / --dmf-tol.
                                  [default: (gau)]
  --thresh-gsm [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the GSM string
                                  optimizer of the MEP stage.  [default:
                                  (gau_loose)]
  --dmf-tol, --thresh-dmf TEXT    IPOPT dual-infeasibility tolerance for the DMF
                                  MEP stage: tight (0.04) | middle (0.10) |
                                  loose (0.20) or a positive float.  [default:
                                  (tight)]
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for TS and post-IRC
                                  endpoint optimizations.  [default: baker]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --resume / --no-resume          Resume from the first unfinished segment using
                                  the verified MEP in --out-dir.  [default: no-
                                  resume]
  --resume-segment INTEGER RANGE  Reuse the verified MEP in --out-dir and rerun
                                  post-processing from segment N. Repeat the
                                  original path/extraction/calculator options;
                                  post-processing options other than
                                  --precision-post may be changed.  [default:
                                  (disabled); x>=1]
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Run input preparation and preflight checks in
                                  a temporary directory, print the execution
                                  plan, and skip calculation stages.  [default:
                                  no-dry-run]
  --preopt / --no-preopt          Run initial single-structure optimizations of
                                  input structures.  [default: preopt]
  --hessian-calc-mode [analytical|finitedifference]
                                  How the backend computes the Hessian in tsopt,
                                  irc and freq.  [default: (FiniteDifference)]
  --detect-layer / --no-detect-layer
                                  Automatically detect ML/MM layers from input
                                  PDB B-factors (ML=0, MovableMM=10,
                                  FrozenMM=20) in downstream tools.  [default:
                                  detect-layer]
  --tsopt / --no-tsopt            Run TS optimization and IRC per reactive
                                  segment (or TSOPT-only mode for a single
                                  structure) and build energy diagrams.
                                  [default: no-tsopt]
  --tsopt-from-mep-tan / --no-tsopt-from-mep-tan
                                  Use the MEP tangent at the highest-energy
                                  image to choose which Hessian mode RS-P-RFO
                                  (hess) follows as the reaction mode. Not used
                                  with Dimer (grad).  [default: tsopt-from-mep-
                                  tan]
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
                                  loop in tsopt; --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --reject-uphill / --no-reject-uphill
                                  Reject uphill RFO steps during post-IRC
                                  endpoint re-optimization only (tolerance: 1e-4
                                  Hartree). Does not affect TS optimization or
                                  path search.  [default: no-reject-uphill]
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
                                  Skip the final PHVA/frequency analysis in
                                  tsopt. The TS structure is retained with
                                  unverified saddle order, and IRC and the later
                                  stages of that segment are skipped; other
                                  segments continue.  [default: no-skip-final-
                                  freq]
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
                                  [default: (<segment>/ts)]
  --freq-out-dir DIRECTORY        Override freq output base directory (relative
                                  paths resolved against the default).
                                  [default: (<segment dir>/freq)]
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
                                  [default: (<segment dir>/dft)]
  --func-basis, --dft-func-basis TEXT
                                  Override dft --func-basis value.  [default:
                                  (wb97m-v/def2-svp)]
  --scf-max-cycles, --dft-max-cycle INTEGER RANGE
                                  Override the DFT SCF iteration limit.
                                  [default: (100); x>=1]
  --scf-tol, --dft-conv-tol FLOAT
                                  Override the DFT SCF convergence tolerance.
                                  [default: (1e-09)]
  --dft-grid-level INTEGER        Override dft --grid-level value.  [default:
                                  (3)]
  --dft-engine, --engine [gpu|cpu]
                                  Override the DFT execution engine.  [default:
                                  (gpu)]
  -s, --scan-lists TEXT           Scan targets: distance (i,j,target), angle
                                  (i,j,k,target), or dihedral (i,j,k,l,target).
                                  Multiple inline literals define sequential
                                  stages. Distances use Å; angles and dihedrals
                                  use degrees. Indices refer to the original
                                  full PDB (1-based) or PDB atom selectors.
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
  --scan-restraint-k, --scan-bias-k FLOAT
                                  Override scan harmonic bias strength k
                                  (eV/Å^2).  [default: (300.0)]
  --scan-relax-max-cycles INTEGER RANGE
                                  Override scan relaxation max cycles per step.
                                  [default: (100000); x>=1]
  --scan-preopt / --no-scan-preopt
                                  Override scan --preopt flag.  [default:
                                  (inherits --preopt)]
  --scan-endopt / --no-scan-endopt
                                  Override scan --endopt flag.  [default:
                                  (False)]
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs to PDB format using
                                  reference topology; forwarded to all
                                  subcommands. --no-convert-files cannot be
                                  combined with --tsopt, --thermo, or --dft.
                                  [default: convert-files]
  --ref-pdb FILE                  Reference PDB for topology/B-factor layer
                                  information when -i provides XYZ inputs.
  -b, --backend [uma|orb|mace|aimnet2|dft]
                                  High-level backend for the ONIOM model region.
                                  [default: (uma)]
  --embedcharge / --no-embedcharge
                                  Enable point-charge treatment. MLIP/MM stages
                                  use the computationally expensive xTB delta
                                  correction; DFT/MM stages embed MM charges in
                                  the PySCF Hamiltonian.  [default: no-
                                  embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from the ML region for the
                                  embedded MM point charges, measured on the -i
                                  PDB (or the --ref-pdb structure when given) so
                                  that all states share one point-charge set. A
                                  cutoff lowers the cost, but the default (none:
                                  every MM atom except those one or two bonds
                                  from the ML region) is the correct
                                  electrostatic embedding and is recommended.
                                  [default: (none)]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor) or
                                  fixed (1.09/1.01 Å).  [default: (scaled)]
  --mm-backend [hessian_ff|openmm]
                                  MM backend: 'hessian_ff' (the MM engine
                                  bundled with mlmm-toolkit; default) or
                                  'openmm'. MM Hessians use finite differences
                                  by default; set calc.mm_fd: false for the
                                  hessian_ff analytical path.  [default:
                                  (hessian_ff)]
  --cmap / --no-cmap              Preserve CMAP terms in both real and model MM
                                  layers when present in parm7.  [default:
                                  (cmap)]
  --coord-type [cart|dlc]         Optimization coordinate system.  [default:
                                  (cart)]
  --print-every INTEGER RANGE     Print optimizer status every N cycles.
                                  [default: (100); x>=1]
  --precision [tf32|fp32|fp64]    MLIP backend precision: tf32, fp32 (IEEE) or
                                  fp64. tf32 uses TF32 matmul for speed (uma,
                                  orb only); avoid it for Hessians. aimnet2:
                                  fp32 only. --precision-post overrides it for
                                  TSOPT, IRC, endpoint optimization and freq.
                                  [default: (per backend: uma fp32; orb, mace
                                  fp64)]
  --precision-post [tf32|fp32|fp64]
                                  MLIP precision for TSOPT, IRC, the post-IRC
                                  endpoint optimizations and freq. Unset uses
                                  --precision.  [default: (inherits
                                  --precision)]
  --compile / --no-compile        torch.compile the UMA or ORB model: slower
                                  first call, slightly faster afterwards; worth
                                  it only for long runs.  [default: (off)]
  --uma-workers, --workers INTEGER
                                  MLIP predictor workers (UMA). >1 uses a
                                  parallel predictor (needs fairchem-core[ray]);
                                  combining it with an analytical Hessian is an
                                  error. Default 1.  [default: (1)]
  --uma-workers-per-node, --workers-per-node INTEGER
                                  Workers per node when the parallel MLIP
                                  predictor is used (--workers > 1).  [default:
                                  (1)]
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p2p1 / uma-s-1p2p1-ft-d4rev / uma-m-1p1
                                  for uma, orbmol-v2 / orb_v3_conservative_omol
                                  for orb, MACE-OMOL-0 / MACE-POLAR-1-M for
                                  mace, aimnet2 / aimnet2-nse for aimnet2).
                                  [default: (the selected backend's own model)]
  --calc-file FILE                Python file exposing get_calculator(...) -> an
                                  ASE Calculator for the ML region (overrides
                                  --backend). Works with GFN-xTB, DFTB+ or any
                                  ASE engine. See --calc-file-func-name.
  --calc-factory, --calc-file-func-name TEXT
                                  Name of the callable in --calc-file that
                                  returns an ASE Calculator (or a module-level
                                  Calculator instance).  [default:
                                  (get_calculator)]
  --deterministic / --no-deterministic
                                  Make repeated runs on the same software
                                  versions and hardware give identical results
                                  (PyTorch deterministic algorithms; AM1-BCC
                                  ligand charges are not covered, so reuse one
                                  parm7). Slower; raises an error on unsupported
                                  ops; custom calculators (--calc-file) and
                                  AIMNet2 are not supported: the run stops.
                                  [default: no-deterministic]
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). An open-shell ML region needs a
                                  matching multiplicity; use this only for an
                                  intentional nonstandard input such as a
                                  covalently-cut region.
  --save-scf-checkpoint / --no-save-scf-checkpoint
                                  Save a structure-bound PySCF checkpoint.
                                  [default: (disabled)]
  --scf-checkpoint FILE           Load/save the optional structure-bound PySCF
                                  checkpoint.
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
