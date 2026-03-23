# `mlmm all`

```text
mlmm-toolkit ver. 0.2.5.dev18

Usage: mlmm all [OPTIONS]

  Run pocket extraction → (optional single-structure staged scan) → MEP search
  in one shot. If exactly one input is provided: (a) with --scan-lists, stage
  results feed into path_search; (b) with --tsopt True and no --scan-lists, run
  TSOPT-only mode.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more **full** PDBs in reaction order
                                  (reactant [intermediates ...] product), or a
                                  single **full** PDB (with --scan-lists or with
                                  --tsopt True). You may pass a single '-i'
                                  followed by multiple space-separated files
                                  (e.g., '-i A.pdb B.pdb C.pdb').  [required]
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
  --include-h2o BOOLEAN           Include waters (HOH/WAT/TIP3/SOL) in the
                                  pocket.  [default: True]
  --exclude-backbone BOOLEAN      Remove backbone atoms on non‑substrate amino
                                  acids (with PRO/HYP safeguards).  [default:
                                  False]
  --add-linkh BOOLEAN             Add link hydrogens for severed bonds (carbon-
                                  only) in pockets.  [default: False]
  --selected-resn TEXT            Force-include residues (comma/space separated;
                                  chain/insertion codes allowed).  [default: ""]
  -l, --ligand-charge TEXT        Either a total charge (number) to distribute
                                  across unknown residues or a mapping like
                                  'GPP:-3,MMT:-1'.
  -q, --charge INTEGER            Force total system charge. Highest priority
                                  over derived charges.
  --parm FILE                     Pre-built AMBER parm7 topology file. When
                                  provided, mm_parm generation is skipped.
  --model-pdb FILE                Pre-built ML-region PDB (with B-factor layer
                                  info). When provided, ml_region generation is
                                  skipped.
  --auto-mm-ff-set [ff19sb|ff14sb]
                                  Force-field set forwarded to mm_parm (ff19SB
                                  uses OPC3; ff14SB uses TIP3P).  [default:
                                  ff19SB]
  --auto-mm-add-ter / --auto-mm-no-add-ter
                                  Control mm_parm TER insertion around
                                  ligand/water/ion blocks.  [default: auto-mm-
                                  add-ter]
  --auto-mm-keep-temp             Keep the mm_parm temporary working directory
                                  (for debugging).
  --auto-mm-ligand-mult TEXT      Spin multiplicity mapping forwarded to mm_parm
                                  (e.g., 'GPP:2,SAM:1'). If omitted, mm_parm
                                  defaults to 1 for all ligands.
  --verbose BOOLEAN               Enable INFO-level logging inside extractor.
                                  [default: True]
  -m, --multiplicity INTEGER      Multiplicity (2S+1).  [default: 1]
  --max-nodes INTEGER             Max internal nodes for **segment** GSM (String
                                  has max_nodes+2 images including endpoints).
                                  [default: 20]
  --max-cycles INTEGER            Maximum GSM optimization cycles.  [default:
                                  300]
  --climb BOOLEAN                 Enable transition-state climbing after growth
                                  for the **first** segment in each pair.
                                  [default: True]
  --opt-mode [grad|hess]          Optimizer mode forwarded to scan/path-search
                                  and used for single optimizations: grad
                                  (=LBFGS/Dimer) or hess (=RFO/RSIRFO).
                                  [default: grad]
  --opt-mode-post [grad|hess]     Optimizer mode for TSOPT and post-IRC endpoint
                                  optimizations. Takes precedence over --opt-
                                  mode for these stages.  [default: hess]
  --dump BOOLEAN                  Dump GSM / single-structure trajectories
                                  during the run, forwarding the same flag to
                                  scan/tsopt/freq.  [default: False]
  --refine-path / --no-refine-path
                                  If True, run recursive path_search on the full
                                  ordered series; if False, run a single-pass
                                  path-opt GSM between each adjacent pair and
                                  concatenate the segments (no path_search).
                                  [default: refine-path]
  --thresh TEXT                   Convergence preset (gau_loose|gau|gau_tight|ga
                                  u_vtight|baker|never). Defaults to 'gau' when
                                  not provided.
  --thresh-post TEXT              Convergence preset for post-IRC endpoint
                                  optimizations (gau_loose|gau|gau_tight|gau_vti
                                  ght|baker|never).  [default: baker]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running any stage.  [default: no-dry-
                                  run]
  --preopt BOOLEAN                If True, run initial single-structure
                                  optimizations of the pocket inputs.  [default:
                                  True]
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
  --tsopt BOOLEAN                 TS optimization + EulerPC IRC per reactive
                                  segment (or TSOPT-only mode for single-
                                  structure), and build energy diagrams.
                                  [default: False]
  --thermo BOOLEAN                Run freq on (R,TS,P) per reactive segment (or
                                  TSOPT-only mode) and build Gibbs free-energy
                                  diagram (MLIP).  [default: False]
  --dft BOOLEAN                   Run DFT single-point on (R,TS,P) and build DFT
                                  energy diagram. With --thermo True, also
                                  generate a DFT//MLIP Gibbs diagram.  [default:
                                  False]
  --tsopt-max-cycles INTEGER      Override tsopt --max-cycles value.
  --flatten / --no-flatten        Enable the extra-imaginary-mode flattening
                                  loop in tsopt (grad: dimer loop, hess: post-
                                  RSIRFO); --no-flatten forces
                                  flatten_max_iter=0.  [default: no-flatten]
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
  --dft-out-dir DIRECTORY         Override dft output base directory (relative
                                  paths resolved against the default).
  --dft-func-basis TEXT           Override dft --func-basis value.
  --dft-max-cycle INTEGER         Override dft --max-cycle value.
  --dft-conv-tol FLOAT            Override dft --conv-tol value.
  --dft-grid-level INTEGER        Override dft --grid-level value.
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
  --scan-one-based BOOLEAN        Override scan indexing interpretation (True =
                                  1-based, False = 0-based).
  --scan-max-step-size FLOAT      Override scan --max-step-size (Å).
  --scan-bias-k FLOAT             Override scan harmonic bias strength k
                                  (eV/Å^2).
  --scan-relax-max-cycles INTEGER
                                  Override scan relaxation max cycles per step.
  --scan-preopt BOOLEAN           Override scan --preopt flag.
  --scan-endopt BOOLEAN           Override scan --endopt flag.
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
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from ML region for MM
                                  point charges in xTB embedding. Default: 12.0
                                  Å. Only used when --embedcharge is enabled.
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  --mm-backend [hessian_ff|openmm]
                                  MM backend: hessian_ff (analytical Hessian,
                                  default) or openmm (finite-difference Hessian,
                                  slower).
  --cmap / --no-cmap              Enable CMAP (backbone cross-map) terms in
                                  model parm7. Default: disabled (Gaussian
                                  ONIOM-compatible).
  -h, --help                      Show this message and exit.
```
