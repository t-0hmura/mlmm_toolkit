# `mlmm oniom-export`

```text
Usage: mlmm oniom-export [OPTIONS]

  Export ONIOM input from Amber parm7 topology (Gaussian g16 or ORCA).

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+optimizer cycle tables,
                                  per-stage timing, VRAM, deliverable paths;
                                  3=everything (full config blocks, per-file
                                  paths, DEBUG logging).  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  --parm FILE                     Amber parm7 topology file.  [required]
  -i, --input FILE                Coordinate file (.pdb or .xyz) for the current
                                  structure (atom order must match parm7).
  --element-check / --no-element-check
                                  Validate that the element sequence in --input
                                  matches the parm7 topology.  [default:
                                  element-check]
  --model-pdb FILE                PDB file defining QM region atoms.
  -o, --output FILE               Output file path (.gjf/.com for g16, .inp for
                                  ORCA when --mode is omitted).  [required]
  --mode [g16|orca]               Export mode. If omitted, inferred from -o
                                  suffix: .gjf/.com -> g16, .inp -> orca.
  --method TEXT                   QM method and basis set. Defaults depend on
                                  mode.
  -q, --charge INTEGER            Charge of QM region.  [required]
  -m, --multiplicity INTEGER      Multiplicity of QM region.  [default: 1]
  --near FLOAT                    Distance cutoff for movable/active atoms
                                  (Angstrom).  [default: 6.0]
  --nproc INTEGER                 Number of processors.  [default: 8]
  --mem TEXT                      Memory allocation (g16 mode).  [default: 16GB]
  --total-charge INTEGER          Total charge of full QM+MM system for ORCA
                                  Charge_Total (orca mode).
  --total-mult INTEGER            Total multiplicity of full QM+MM system for
                                  ORCA Mult_Total (orca mode).
  --orcaff PATH                   Path to ORCAFF.prms (orca mode). If omitted,
                                  uses/creates <parm7_stem>.ORCAFF.prms in
                                  output directory.
  --convert-orcaff / --no-convert-orcaff
                                  If ORCAFF.prms is missing, try `orca_mm
                                  -convff -AMBER` automatically (orca mode).
                                  [default: convert-orcaff]
  --link-atom-method [scaled|fixed]
                                  Link-H placement rule. 'scaled' uses the
                                  Morokuma/Dapprich g-factor (matches MLMMCore
                                  runtime). 'fixed' uses the legacy 1.09/1.01 Å
                                  bond length.  [default: scaled]
  -h, --help                      Show this message and exit.
```
