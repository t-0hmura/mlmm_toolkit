# `mlmm dft`

```text

Usage: mlmm dft [OPTIONS]

  Single-point ML-region DFT with ML(dft)/MM energy recombination.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Full enzyme structure (PDB or XYZ). If XYZ,
                                  use --ref-pdb for topology.  [required]
  --ref-pdb FILE                  Reference PDB topology when input is XYZ. XYZ
                                  coordinates are used (higher precision) while
                                  PDB provides atom ordering and residue
                                  information.
  --parm FILE                     Amber parm7 topology for the full system.
                                  [required]
  --model-pdb FILE                PDB defining the ML region (atom IDs must
                                  match the enzyme PDB). Optional when --detect-
                                  layer is enabled.
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
  -q, --charge INTEGER            Charge of the ML region.  [required]
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  [default: 1]
  --freeze-atoms TEXT             Comma-separated 1-based indices to freeze
                                  (e.g., '1,3,5').
  --func-basis TEXT               Exchange-correlation functional and basis set
                                  as "FUNC/BASIS".  [default:
                                  wb97m-v/def2-tzvpd]
  --max-cycle INTEGER             Maximum SCF iterations.  [default: 100]
  --conv-tol FLOAT                SCF convergence tolerance (Hartree).
                                  [default: 1e-09]
  --grid-level INTEGER            DFT integration grid level (0=coarse,
                                  3=default, 9=ultrafine).  [default: 3]
  -o, --out-dir DIRECTORY         Output directory.  [default: result_dft]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --show-config / --no-show-config
                                  Print resolved configuration and continue
                                  execution.  [default: no-show-config]
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running DFT.  [default: no-dry-run]
  --convert-files / --no-convert-files
                                  Toggle XYZ/TRJ to PDB companions when a PDB
                                  template is available.  [default: convert-
                                  files]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable PySCF qmmm.mm_charge() point-charge
                                  embedding of MM partial charges into the DFT
                                  Hamiltonian.  [default: no-embedcharge]
  -h, --help                      Show this message and exit.
```
