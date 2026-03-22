# `mlmm all`

```text
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
  -l, --ligand-charge TEXT        Either a total charge (number) to distribute
                                  across unknown residues or a mapping like
                                  'GPP:-3,MMT:-1'.
  -q, --charge INTEGER            Force total system charge. Highest priority
                                  over derived charges.
  --config FILE                   Base YAML configuration file applied before
                                  explicit CLI options.
  --dry-run / --no-dry-run        Validate options and print the execution plan
                                  without running any stage.  [default: no-dry-
                                  run]
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
  -s, --scan-lists TEXT           Scan targets: inline Python literal or a
                                  YAML/JSON spec file path. Multiple inline
                                  literals define sequential stages, e.g.
                                  "[(12,45,1.35)]"
                                  "[(10,55,2.20),(23,34,1.80)]". Indices refer
                                  to the original full PDB (1-based) or PDB atom
                                  selectors like "TYR,285,CA"; they are auto-
                                  mapped to the pocket after extraction.
  -b, --backend [uma|orb|mace|aimnet2]
                                  ML backend for the ONIOM high-level region
                                  (default: uma).
  --embedcharge / --no-embedcharge
                                  Enable xTB point-charge embedding correction
                                  for MM→ML environmental effects.  [default:
                                  no-embedcharge]
  -h, --help                      Show this message and exit.
```
