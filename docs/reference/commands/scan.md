# `mlmm scan`

```text

Usage: mlmm scan [OPTIONS]

  Bond-length driven scan with staged harmonic restraints and relaxation
  (ML/MM).

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Full-enzyme PDB used by the ML/MM calculator
                                  and as reference for conversions.  [required]
  --parm FILE                     Amber parm7 topology covering the entire
                                  enzyme complex.  [required]
  --model-pdb FILE                PDB defining the ML-region atoms for ML/MM.
                                  Optional when --detect-layer is enabled.
  --model-indices TEXT            Comma-separated atom indices for the ML region
                                  (ranges allowed like 1-5). Used when --model-
                                  pdb is omitted.
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (B=0/10/20). If disabled, you must provide
                                  --model-pdb or --model-indices.  [default:
                                  detect-layer]
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided.
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  GPP:-3,SAM:1) used to derive charge when -q is
                                  omitted (requires PDB input or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  Defaults to 1 when omitted.
  --freeze-atoms TEXT             Comma-separated 1-based atom indices to freeze
                                  (e.g., '1,3,5').
  --hess-cutoff FLOAT             Distance cutoff (Å) from ML region for MM
                                  atoms to include in Hessian calculation.
                                  Applied to movable MM atoms and can be
                                  combined with --detect-layer.
  --movable-cutoff FLOAT          Distance cutoff (Å) from ML region for movable
                                  MM atoms. MM atoms beyond this are frozen.
                                  Providing --movable-cutoff disables --detect-
                                  layer.
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
  --bias-k FLOAT                  Harmonic well strength k [eV/Å^2].  [default:
                                  300]
  --opt-mode [grad|hess|lbfgs|rfo|light|heavy]
                                  Compatibility option for mlmm all forwarding.
                                  scan relaxations currently use LBFGS
                                  regardless of this value.
  --max-cycles INTEGER            Maximum LBFGS cycles per biased step and per
                                  (pre|end)opt stage.  [default: 10000]
  --relax-max-cycles INTEGER      Compatibility alias of --max-cycles (overrides
                                  it when provided).
  --dump / --no-dump              Write per-step optimizer trajectory files.
                                  scan_trj.xyz is always written regardless.
                                  [default: no-dump]
  --out-dir TEXT                  Base output directory.  [default:
                                  ./result_scan/]
  --thresh [gau_loose|gau|gau_tight|gau_vtight|baker|never]
                                  Convergence preset for relaxations.
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
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  -h, --help                      Show this message and exit.
```
