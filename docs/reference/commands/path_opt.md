# `mlmm path-opt`

```text
Usage: mlmm path-opt [OPTIONS]

  MEP optimization via the Growing String method or Direct Max Flux.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE...             Two endpoint structures (reactant/product);
                                  both must be full-enzyme PDBs.  [required]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1). Defaults to 1 when
                                  omitted.
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --dmf-backend [cpu|gpu]         DMF compute backend (--mep-mode dmf only): gpu
                                  (dmf.torch / CUDA) or cpu (dmf / NumPy). On a
                                  GPU out-of-memory, retry with cpu.  [default:
                                  gpu]
  --max-nodes INTEGER             Number of internal nodes (for GSM: string has
                                  max_nodes+2 images including endpoints; for
                                  DMF: number of path waypoints).  [default: 20]
  --max-cycles INTEGER            Maximum optimization cycles.  [default: 300]
  --climb / --no-climb            Search for a transition state (climbing image)
                                  after path growth.  [default: climb]
  --preopt / --no-preopt          Pre-optimize the two endpoint structures with
                                  LBFGS before string growth.  [default: preopt]
  --preopt-max-cycles INTEGER     Maximum LBFGS cycles for endpoint pre-
                                  optimization when --preopt is enabled.
                                  [default: 10000]
  --fix-ends / --no-fix-ends      Fix endpoint structures during path growth.
                                  [default: fix-ends]
  --dump / --no-dump              Dump optimizer trajectory/restarts during the
                                  run.  [default: no-dump]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_opt/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for the string optimizer.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running path optimization.  [default:
                                  no-dry-run]
  --parm FILE                     Amber parm7 topology for the enzyme complex
                                  (MM layers).  [required]
  --model-pdb FILE                PDB defining the ML region (atom IDs used by
                                  the ML/MM calculator). Optional when --detect-
                                  layer is enabled.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  --freeze-atoms TEXT             Comma-separated 1-based indices to freeze
                                  (applied to every image).
  --hess-cutoff FLOAT             Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer.
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.
  --convert-files / --no-convert-files
                                  Convert XYZ/TRJ outputs into PDB companions
                                  based on the input format.  [default: convert-
                                  files]
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
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --ref-pdb FILE                  Full-size template PDBs in the same order as
                                  --input. Required when using XYZ inputs to
                                  provide topology and B-factor information.
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
  --backend-model TEXT            Model variant for the selected --backend (e.g.
                                  uma-s-1p1 / uma-m-1p1 for uma,
                                  orb_v3_conservative_omol for orb, MACE-OMOL-0
                                  / MACE-OFF23_small for mace). Default: the
                                  backend's built-in model.
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
