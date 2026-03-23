# mlmm irc

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli irc [OPTIONS]

  Run an IRC calculation with EulerPC. Only the documented CLI options are
  accepted; all other settings come from YAML.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input FILE                Input structure file (.pdb, .xyz, _trj.xyz,
                                  etc.).  [required]
  --parm FILE                     Amber parm7 topology for the whole enzyme (MM
                                  region). If omitted, must be provided in YAML
                                  as calc.real_parm7.
  --model-pdb FILE                PDB defining atoms belonging to the ML region.
                                  Optional when --detect-layer is enabled.
  --detect-layer / --no-detect-layer
                                  Detect ML/MM layers from input PDB B-factors
                                  (ML=0, MovableMM=10, FrozenMM=20). If
                                  disabled, you must provide --model-pdb or
                                  --model-indices.  [default: detect-layer]
  -q, --charge INTEGER            Total charge; overrides calc.charge from YAML.
                                  Required unless --ligand-charge is provided.
  -m, --multiplicity INTEGER      Spin multiplicity (2S+1); overrides calc.spin
                                  from YAML.
  --max-cycles INTEGER            Maximum number of IRC steps; overrides
                                  irc.max_cycles from YAML.
  --step-size FLOAT               Step length in Bohr (unweighted Cartesian
                                  coordinates). Default: 0.10 Bohr. Overrides
                                  irc.step_length from YAML.
  --forward / --no-forward        Run the forward IRC; overrides irc.forward
                                  from YAML.
  --backward / --no-backward      Run the backward IRC; overrides irc.backward
                                  from YAML.
  -o, --out-dir TEXT              Output directory; overrides irc.out_dir from
                                  YAML.  [default: ./result_irc/]
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
