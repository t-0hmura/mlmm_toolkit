# `mlmm dft`

```text
Usage: mlmm dft [OPTIONS]

  Single-point ML-region DFT with ML(dft)/MM energy recombination.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
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
  -q, --charge INTEGER            ML region charge. Required unless --ligand-
                                  charge is provided (PDB inputs or XYZ with
                                  --ref-pdb).
  -l, --ligand-charge TEXT        Total charge or per-resname mapping (e.g.,
                                  'SAM:1,GPP:-3') used to derive ML region
                                  charge when -q is omitted (requires PDB input
                                  or --ref-pdb).
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  [default: 1]
  --freeze-atoms TEXT             Comma-separated 1-based indices to freeze
                                  (e.g., '1,3,5').
  --func-basis TEXT               Exchange-correlation functional and basis set
                                  as "FUNC/BASIS".  [default:
                                  wb97m-v/def2-tzvpd]
  --max-cycle INTEGER             Maximum SCF iterations.  [default: 100]
  --conv-tol FLOAT                SCF energy convergence threshold (ΔE in
                                  Hartree between SCF cycles).  [default: 1e-09]
  --grid-level INTEGER            DFT integration grid level (0=coarse,
                                  3=default, 5=fine, 9=very fine).  [default: 3]
  --engine [gpu|cpu]              SCF backend: gpu (GPU4PySCF, raises error if
                                  unavailable) or cpu (PySCF).  [default: gpu]
  --lowmem / --no-lowmem          Use gpu4pyscf rks_lowmem.RKS for closed-shell
                                  GPU runs (memory-efficient direct JK; mlmm dft
                                  does not call density_fit() on either path).
                                  Open-shell or CPU engines fall back to
                                  standard RKS/UKS automatically.  [default:
                                  lowmem]
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
                                  Backend label recorded in output metadata; the
                                  ML region in dft is computed with DFT (PySCF),
                                  so this does not select a calculator.
  --embedcharge / --no-embedcharge
                                  Enable electrostatic embedding: MM point
                                  charges are added to the PySCF QM Hamiltonian
                                  via pyscf.qmmm.mm_charge().  [default: no-
                                  embedcharge]
  --embedcharge-cutoff FLOAT      Distance cutoff (Å) from ML region for MM
                                  point charges embedded in the PySCF QM
                                  Hamiltonian. Default: 12.0 Å. Only used when
                                  --embedcharge is enabled.
  --link-atom-method [scaled|fixed]
                                  Link-atom placement: 'scaled' (g-factor,
                                  Gaussian ONIOM standard, default) or 'fixed'
                                  (legacy 1.09 Å for C, 1.01 Å for N).
  --mm-backend [hessian_ff|openmm]
                                  MM backend for the low-level ONIOM evaluation:
                                  'hessian_ff' (default) or 'openmm'.
  --cmap / --no-cmap              Enable CMAP (backbone cross-map) terms in
                                  model parm7. Default: disabled (Gaussian
                                  ONIOM-compatible).
  --out-json / --no-out-json      Write machine-readable result.json to out_dir.
                                  [default: no-out-json]
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  --model-indices-one-based / --model-indices-zero-based
                                  Interpret --model-indices as 1-based (default)
                                  or 0-based.  [default: model-indices-one-
                                  based]
  -h, --help                      Show this message and exit.
```
