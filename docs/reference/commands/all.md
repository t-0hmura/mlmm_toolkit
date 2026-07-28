# `mlmm all`

```text
Usage: mlmm all [OPTIONS]

  Run pocket extraction → (optional single-structure staged scan) → MEP search
  in one shot. If exactly one input is provided: (a) with --scan-lists, stage
  results feed into path-opt (or path_search with --refine-path); (b) with
  --tsopt and no --scan-lists, run TSOPT-only mode.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more full PDB/mmCIF structures in
                                  reaction order (reactant [intermediates ...]
                                  product), or one full structure with --scan-
                                  lists or --tsopt. A single '-i' may be
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
  -r, --radius FLOAT              Inclusion cutoff (Å) around substrate atoms.
                                  [default: 2.6]
  --radius-het2het FLOAT          Independent hetero–hetero cutoff (Å) for
                                  non‑C/H pairs.  [default: 0.0]
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
  --selected-resn TEXT            Force-include residues (comma/space separated;
                                  chain/insertion codes allowed).  [default: ""]
  --modified-residue TEXT         Comma-separated residue names (with optional
                                  charge) to treat as amino acids for backbone
                                  truncation and charge assignment. Examples:
                                  'HD1,HD2,HD3' or 'HD1:0,SEP:-2'.  [default:
                                  ""]
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
                                  precedence over an ML selection produced by
                                  -c/--center.
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
                                  defaults to 1 for all ligands.
  --auto-mm-disulfide / --auto-mm-no-disulfide
                                  Forwarded to mm_parm: detect disulfides from
                                  SG-SG geometry across CYS/CYM/CYX and bond
                                  them (renaming a bonded CYS to CYX). With
                                  --auto-mm-no-disulfide only residues already
                                  named CYX are bonded and CYS is left
                                  untouched.  [default: auto-mm-disulfide]
  -m, --multiplicity INTEGER      Multiplicity (2S+1).  [default: 1]
  --tr-projection [constrained|legacy-active]
                                  Rigid translation/rotation treatment forwarded
                                  to TSopt, IRC, freq, and flatten PHVA. The
                                  default respects frozen anchors; 'legacy-
                                  active' is deprecated and must not be used for
                                  pass/HOSP transition-state certification.
                                  [default: constrained]
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory error, retry with cpu.
                                  [default: gpu]
  --max-nodes INTEGER             Max internal nodes per GSM/DMF segment
                                  (max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-cycles INTEGER            Maximum MEP optimization cycles.  [default:
                                  300]
  --climb / --no-climb            Enable transition-state climbing after growth
                                  for the *first* segment in each pair.
                                  [default: climb]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/path-search
                                  and used for single optimizations: grad
                                  (=L-BFGS/Dimer) or hess (=RFO/RSIRFO).
                                  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode for TSOPT and post-IRC endpoint
                                  optimizations. Takes precedence over --opt-
                                  mode for these stages.  [default: hess]
  --dump / --no-dump              Dump MEP / single-structure trajectories
                                  during the run, forwarding the same flag to
                                  scan/tsopt/freq.  [default: no-dump]
  --refine-path / --no-refine-path
                                  If False (default), run single-pass path-opt
                                  with the selected MEP optimizer between each
                                  adjacent pair and concatenate the segments (no
                                  path_search); if True, run recursive
                                  path_search on the full ordered series for
                                  automatic multistep discovery.  [default: no-
                                  refine-path]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau_loose'
                                  for path-opt, 'gau' for scan.
  --thresh-post [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for post-IRC endpoint
                                  optimizations (gau_loose|gau|gau_tight|gau_vti
                                  ght|baker|never).  [default: baker]
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
                                  Common MLIP Hessian calculation mode forwarded
                                  to tsopt and freq. Default:
                                  'FiniteDifference'. Use 'Analytical' when VRAM
                                  is sufficient.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20) in
                                  downstream tools. If disabled, downstream
                                  tools require --model-pdb or --model-indices.
                                  [default: detect-layer]
  --tsopt / --no-tsopt            TS optimization + EulerPC IRC per reactive
                                  segment (or TSOPT-only mode for single-
                                  structure), and build energy diagrams.
                                  [default: no-tsopt]
  --thermo / --no-thermo          Run freq on (R,TS,P) per reactive segment (or
                                  TSOPT-only mode) and build Gibbs free-energy
                                  diagram (MLIP).  [default: no-thermo]
  --dft / --no-dft                Run DFT single-point on (R,TS,P) and build a
                                  DFT energy diagram. With --thermo, also
                                  generate a DFT//MLIP/MM Gibbs diagram.
                                  [default: no-dft]
  --tsopt-max-cycles INTEGER      Override tsopt --max-cycles value.
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt (grad: dimer loop, hess: post-
                                  RSIRFO); --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
  --reject-uphill / --no-reject-uphill
                                  Reject uphill RFO trials during post-IRC
                                  endpoint re-optimization only and final-check
                                  the retained endpoint at the emergency floor.
                                  Does not affect TS optimization or path
                                  search.  [default: reject-uphill]
  --irc-step-size FLOAT           Override IRC --step-size (Bohr). If an IRC
                                  stops after only a few frames, retry with a
                                  smaller value such as 0.05.
  --irc-never-stop / --no-irc-never-stop
                                  Forward IRC never-stop mode to every post-TS
                                  IRC. It ignores energy-rise/plateau stops but
                                  retains physical/integrator stops; default
                                  follows irc.never_stop (off).
  --skip-final-freq / --no-skip-final-freq
                                  Skip post-convergence frequency analysis in
                                  tsopt. Useful for large unfrozen systems.
                                  [default: no-skip-final-freq]
  --tsopt-out-dir DIRECTORY       Override tsopt output subdirectory (relative
                                  paths are resolved against the default).
  --freq-out-dir DIRECTORY        Override freq output base directory (relative
                                  paths resolved against the default).
  --freq-max-write INTEGER        Override freq --max-write value.
  --freq-amplitude-ang FLOAT      Override freq --amplitude-ang (Å).
  --freq-n-frames INTEGER         Override freq --n-frames value.
  --freq-sort [value|abs]         Override freq mode sorting.
  --freq-temperature FLOAT        Override freq thermochemistry temperature (K).
  --freq-pressure FLOAT           Override freq thermochemistry pressure (atm).
  --freq-symmetry-number INTEGER RANGE
                                  Use one rotational symmetry number for every
                                  R/TS/P frequency job. When omitted, each child
                                  follows its YAML/default setting.  [x>=1]
  --dft-out-dir DIRECTORY         Override dft output base directory (relative
                                  paths resolved against the default).
  --dft-func-basis TEXT           Override dft --func-basis value.
  --dft-max-cycle INTEGER         Override dft --max-cycle value.
  --dft-conv-tol FLOAT            Override dft --conv-tol value.
  --dft-grid-level INTEGER        Override dft --grid-level value.
  --dft-engine [gpu|cpu]          Override dft --engine value.
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. Multiple inline
                                  literals define sequential stages, e.g.
                                  "[(12,45,1.35)]"
                                  "[(10,55,2.20),(23,34,1.80)]". Indices refer
                                  to the original full PDB (1-based) or PDB atom
                                  selectors like "TYR,285,CA"; they are auto-
                                  mapped to the pocket after extraction.
  --scan-out-dir DIRECTORY        Override the scan output directory (default:
                                  <out-dir>/scan/). Relative paths are resolved
                                  against the default parent.
  --scan-one-based / --scan-zero-based
                                  Override scan indexing interpretation (one-
                                  based or zero-based).
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å).
  --scan-bias-k FLOAT             Override scan harmonic bias strength k
                                  (eV/Å^2).
  --scan-relax-max-cycles INTEGER
                                  Override scan relaxation max cycles per step.
  --scan-preopt / --no-scan-preopt
                                  Override scan --preopt flag.
  --scan-endopt / --no-scan-endopt
                                  Override scan --endopt flag.
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
  --coord-type [cart|dlc]         Optimization coordinate system (cart|dlc).
                                  cart is the robust default used in published
                                  numbers; dlc speeds up torsion-rich opts.
                                  mlmm-specific caveats: DLC + link atom and DLC
                                  + 3-layer frozen MM are numerically
                                  unverified.
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
                                  Strict bit-reproducible GPU runs
                                  (deterministic algorithms + index_reduce_
                                  shim). Slower; raises if unsupported. Default
                                  off.
  --allow-charge-mult-mismatch    Skip the ML-region charge/multiplicity
                                  electron-parity check (logs that it was
                                  skipped). For an intentional open-shell or
                                  covalently-cut ML region.
  -h, --help                      Show this message and exit.
```
