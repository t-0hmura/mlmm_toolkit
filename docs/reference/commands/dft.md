# mlmm dft

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli dft [OPTIONS]

  Single-point ML-region DFT with ML(dft)/MM energy recombination.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Full enzyme structure (PDB or XYZ). If XYZ,
                                  use --ref-pdb for topology.  [required]
  --parm FILE                     Amber parm7 topology for the full system.
                                  [required]
  --model-pdb FILE                PDB defining the ML region (atom IDs must
                                  match the enzyme PDB). Optional when --detect-
                                  layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            Charge of the ML region.  [required]
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1) for the ML region.
                                  [default: 1]
  --func-basis TEXT               Exchange-correlation functional and basis set
                                  as "FUNC/BASIS".  [default:
                                  wb97m-v/def2-tzvpd]
  -o, --out-dir DIRECTORY         Output directory.  [default: result_dft]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable electrostatic embedding: MM point
                                  charges are added to the PySCF QM Hamiltonian
                                  via pyscf.qmmm.mm_charge().  [default: no-
                                  embedcharge]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  -h, --help                      Show this message and exit.
```
