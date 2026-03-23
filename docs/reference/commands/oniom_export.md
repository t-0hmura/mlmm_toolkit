# mlmm oniom-export

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli oniom-export [OPTIONS]

  Export ONIOM input from Amber parm7 topology (Gaussian g16 or ORCA).

Options:
  --help-advanced             Show all options (including advanced settings) and
                              exit.
  --parm FILE                 Amber parm7 topology file.  [required]
  -i, --input FILE            Coordinate file (.pdb or .xyz) for the current
                              structure (atom order must match parm7).
  --model-pdb FILE            PDB file defining QM region atoms.
  -o, --output FILE           Output file path (.gjf/.com for g16, .inp for ORCA
                              when --mode is omitted).  [required]
  --mode [g16|orca]           Export mode. If omitted, inferred from -o suffix:
                              .gjf/.com -> g16, .inp -> orca.
  --method TEXT               QM method and basis set. Defaults depend on mode.
  -q, --charge INTEGER        Charge of QM region.  [required]
  -m, --multiplicity INTEGER  Multiplicity of QM region.  [default: 1]
  -h, --help                  Show this message and exit.
```
