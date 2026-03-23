# mlmm path-search

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli path-search [OPTIONS]

  Multistep MEP search via recursive GSM segmentation.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Two or more structures in reaction order.
                                  Either repeat '-i' (e.g., '-i A -i B -i C') or
                                  use a single '-i' followed by multiple space-
                                  separated paths (e.g., '-i A B C').
                                  [required]
  --parm FILE                     Amber parm7 topology covering the full enzyme
                                  complex.  [required]
  --model-pdb FILE                PDB describing atoms that belong to the ML
                                  (high-level) region. Optional when --detect-
                                  layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            Total system charge. Required unless --ligand-
                                  charge is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1). Defaults to 1 when
                                  omitted.
  --mep-mode [gsm|dmf]            MEP method: gsm (GrowingString) or dmf (Direct
                                  Max Flux).  [default: gsm]
  --refine-mode [peak|minima]     Refinement seed around the highest-energy
                                  image: 'peak' uses HEI±1, 'minima' uses
                                  nearest local minima. Defaults to peak for gsm
                                  and minima for dmf.
  --max-nodes INTEGER             Number of internal nodes (string has
                                  max_nodes+2 images including endpoints). Used
                                  for *segment* GSM unless overridden by YAML
                                  search.max_nodes_segment.  [default: 10]
  -o, --out-dir TEXT              Output directory.  [default:
                                  ./result_path_search/]
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
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
