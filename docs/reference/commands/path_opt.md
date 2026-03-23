# mlmm path-opt

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli path-opt [OPTIONS]

  MEP optimization via the Growing String method or Direct Max Flux.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE...             Two endpoint structures (reactant/product);
                                  both must be full-enzyme PDBs.  [required]
  -q, --charge INTEGER            Total charge. Required unless --ligand-charge
                                  is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1). Defaults to 1 when
                                  omitted.
  --mep-mode [gsm|dmf]            MEP optimizer: Growing String Method (gsm) or
                                  Direct Max Flux (dmf).  [default: gsm]
  --max-nodes INTEGER             Number of internal nodes (for GSM: string has
                                  max_nodes+2 images including endpoints; for
                                  DMF: number of path waypoints).  [default: 20]
  --fix-ends / --no-fix-ends      Fix endpoint structures during path growth.
                                  [default: no-fix-ends]
  --out-dir TEXT                  Output directory.  [default:
                                  ./result_path_opt/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --parm FILE                     Amber parm7 topology for the enzyme complex
                                  (MM layers).  [required]
  --model-pdb FILE                PDB defining the ML region (atom IDs used by
                                  the ML/MM calculator). Optional when --detect-
                                  layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  --link-atom-method [scaled|fixed]
                                  Link-atom position mode: scaled (g-factor,
                                  default) or fixed (legacy 1.09/1.01 Å).
  -h, --help                      Show this message and exit.
```
