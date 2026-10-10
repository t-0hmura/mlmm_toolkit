# Glossary

Each abbreviation, method name, and unit used in the docs is defined here in one line, grouped by field. Options, output fields, and status values are on each command page and in [JSON Output](json-output.md).

## ML/MM & ONIOM

| Term | Full Name | Description |
|------|-----------|-------------|
| **ML/MM** | Machine Learning / Molecular Mechanics | A multi-scale method that couples a machine-learning interatomic potential (for the reactive region) with a classical force field (for the surrounding environment). Analogous to QM/MM but with ML replacing QM. |
| **ONIOM** | Our own N-layered Integrated molecular Orbital and molecular Mechanics | A multi-layer energy decomposition scheme. mlmm-toolkit uses an ONIOM-like subtraction: E_total = E_REAL_low + E_MODEL_high - E_MODEL_low. |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics | A multi-scale method coupling QM for the reactive region with MM for the environment. ML/MM replaces the QM layer with an MLIP backend. |
| **Real system** | — | The full set of atoms (all 3 layers). Evaluated at the MM (low) level in the ONIOM decomposition. Described by the parm7 topology; its MM energy is computed by the MM backend. |
| **Model system** | — | The ML region (Layer 1). Evaluated at both the MLIP (high) and MM (low) levels in the ONIOM decomposition. |
| **Link Hydrogen** | — | A hydrogen placed on each parm7 bond that crosses the ML/MM boundary. It lies along that bond, and its force is redistributed through a Jacobian. The cap hydrogens of `extract --add-linkh` are only for inspecting the pocket. |
| **Link atom** | — | See **Link Hydrogen**; in mlmm-toolkit the link atoms placed at severed ML/MM boundaries are hydrogens. |
| **hessian_ff** | — | A C++ native extension that evaluates Amber force field energies, forces, and analytical Hessians. Used as the MM engine in mlmm-toolkit. |
| **3-layer system** | — | mlmm-toolkit's B-factor partitioning scheme: ML (B=0.0), Movable-MM (B=10.0), Frozen-MM (B=20.0). |
| **B-factor encoding** | — | Convention of storing layer membership in the PDB B-factor (temperature factor) column: 0.0 = ML, 10.0 = Movable-MM, 20.0 = Frozen-MM. Hessian-target MM is controlled by cutoffs/explicit indices. See {ref}`The MM layers <mm-layers>`. |

## Amber & Force Field

| Term | Full Name | Description |
|------|-----------|-------------|
| **parm7** | Amber Parameter/Topology File | A file containing atom types, partial charges, bonding connectivity, and force field parameters for an Amber system. Also called .prmtop. |
| **rst7** | Amber Restart File | A file containing atomic coordinates (and optionally velocities and box dimensions) for an Amber system. Also called .inpcrd. |
| **AmberTools** | — | A free suite of tools for molecular dynamics preparation, including tleap, antechamber, and parmchk2. Required by `mlmm mm-parm`. |
| **tleap** | — | An AmberTools program that builds Amber topology/coordinate files from PDB structures and force field libraries. |
| **antechamber** | — | An AmberTools program that assigns GAFF2 atom types and AM1-BCC partial charges to small molecules. |
| **parmchk2** | — | An AmberTools program that checks and supplies missing force field parameters for GAFF2 typing. |
| **GAFF2** | General Amber Force Field 2 | A general-purpose force field for small organic molecules, used to parameterize non-standard residues (substrates, cofactors). |
| **ff19SB** | — | An Amber protein force field used for standard amino acid residues. The default in mlmm-toolkit. |
| **ff14SB** | — | An Amber protein force field (2014 version). Selectable with `--ff-set ff14SB`. |
| **AM1-BCC** | — | A charge model that combines AM1 (semi-empirical) Mulliken charges with bond charge corrections (BCC) to approximate HF/6-31G* RESP charges. |

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway connecting reactants to products through a transition state (on a potential energy surface). |
| **TS** | Transition State | A first-order saddle point on the potential energy surface — a stationary point with exactly one direction of negative curvature (one imaginary frequency) along the reaction coordinate. |
| **n_imag** | Number of imaginary modes | The number of vibrational modes below the imaginary-mode criterion (ν < −5.00 cm⁻¹ by default); a TS has n_imag = 1, and `result.json` records it as `n_imaginary_modes` (`tsopt`) or `n_imaginary` (`freq`). |
| **IRC** | Intrinsic Reaction Coordinate | A mass-weighted steepest-descent path from a TS toward reactants and products. Often used to validate TS connectivity. |
| **GSM** | Growing String Method | A string-based method that grows images from endpoints and optimizes them to approximate an MEP. |
| **DMF** | Direct Max Flux | A chain-of-states method for optimizing an MEP by maximizing flux along the pathway. In mlmm-toolkit it is selected with `--mep-mode dmf`. |
| **HEI** | Highest-Energy Image | The image along an MEP with maximum energy; often used as a TS guess. |
| **Image** | — | A single geometry (one "node") along a chain-of-states path. |
| **Segment** | — | An MEP between two adjacent endpoints (e.g., R → I1, I1 → I2, …). |
| **Kink** | — | A segment where only the conformation changes: the two structures that `path-search` optimizes on either side of the HEI (End1 and End2; [path-search → How it works](path-search.md#how-it-works), step 2) differ by no covalent bond change. `path-search` fills it with a few linearly interpolated nodes (`search.kink_max_nodes`, default 3) and optimizes each one instead of running a new GSM or DMF path. |

## Optimization Algorithms

| Term | Full Name | Description |
|------|-----------|-------------|
| **L-BFGS** | Limited-memory BFGS | A quasi-Newton optimization algorithm that approximates the Hessian using a limited history of gradients. Used in `opt --opt-mode grad`. |
| **RFO** | Rational Function Optimization | A trust-region optimization method that uses explicit Hessian information. Used in `--opt-mode hess`. |
| **RS-I-RFO** | Restricted-Step Image-RFO | A variant of RFO for saddle point (TS) optimization that follows one negative eigenvalue. |
| **Dimer** | Dimer Method | A TS optimization method that follows a low-curvature direction. mlmm-toolkit's Hessian-guided variant uses initial and periodic active-subspace Hessians, which is more robust than a random initial orientation for systems with many active degrees of freedom. Used in `--opt-mode grad` for TSOPT. |
| **PHVA** | Partial Hessian Vibrational Analysis | Computing vibrational frequencies using only the Hessian block for active (non-frozen) atoms. Default in `freq`. |

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Models for Atoms | Meta's family of pretrained MLIPs. The default ML backend in mlmm-toolkit (`--backend uma`). |
| **ORB** | ORB Models | A family of pretrained MLIPs from Orbital Materials. Supported as an alternative ML backend (`--backend orb`). Install with `pip install "mlmm-toolkit[orb]"`. |
| **MACE** | MACE (Message-passing Atomic Cluster Expansion) | A message-passing equivariant neural network MLIP. Supported as an alternative ML backend (`--backend mace`). Install in a **separate** conda env: `pip uninstall fairchem-core` (UMA pin clashes on `e3nn`), then `pip install mace-torch`. |
| **AIMNet2** | Atoms In Molecules Network 2 | A neural network potential for organic molecules. Supported as an alternative ML backend (`--backend aimnet2`). Install with `pip install "mlmm-toolkit[aimnet]"`. |
| **Analytical Hessian** | — | Computing second derivatives through the backend's differentiable/native Hessian path. Runtime and memory are backend- and system-dependent. Supported by UMA, ORB, MACE, and AIMNet2. |
| **Finite Difference** | — | Approximating derivatives from displaced-force evaluations. Runtime and memory are backend- and system-dependent. Available for every MLIP backend. |

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates; used for vibrational analysis and TS optimization. |
| **SP** | Single Point | A calculation at a fixed geometry (no optimization); often used for higher-level energy refinement. |
| **Spin Multiplicity** | — | 2S+1, where S is total spin. Singlet = 1, doublet = 2, triplet = 3, etc. |

## Structural Biology & Pocket Extraction

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. The calculation commands accept XYZ together with `--ref-pdb`. |
| **GJF** | Gaussian Job File | An input format for Gaussian; `oniom-export` writes it in `g16` mode, and `oniom-import` and `bond-summary` read it. |
| **Pocket** | Active-site Pocket | A truncated structure around the substrate(s), extracted by the `extract` subcommand. In the ML/MM workflow, this defines the ML region and surrounding MM environment. |
| **Extractor-only Link Hydrogen** | — | A cap hydrogen that `extract --add-linkh` adds so that you can inspect the pocket. ML/MM calculations do not use it: the link pairs come from the parm7 bonds at the ML/MM boundary. |
| **Backbone** | — | The main chain of a protein (N–Cα–C–O atoms). Can be excluded during pocket extraction with `--exclude-backbone`. |
| **B-factor** | Temperature Factor | The PDB temperature factor column. In mlmm-toolkit, used to encode 3-layer membership (0.0, 10.0, 20.0). |

## Thermochemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **ZPE** | Zero-Point Energy | The vibrational energy at 0 K; a quantum correction to the electronic energy. |
| **Gibbs Energy** | Gibbs Free Energy (G) | G = H − TS; includes thermal and entropic contributions. |
| **Enthalpy** | (H) | H = E + PV; total heat content at constant pressure. |
| **Entropy** | (S) | A measure of disorder; contributes −TS to Gibbs energy. |
| **QRRHO** | Quasi-Rigid-Rotor Harmonic Oscillator | A thermochemical approximation incorporating Grimme's correction for low-frequency vibrations. Automatically applied in `freq`. |

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ≈ 627.5 kcal/mol ≈ 27.21 eV. |
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **kJ/mol** | Kilojoules per mole; 1 kcal/mol ≈ 4.184 kJ/mol. |
| **eV** | Electron volt; 1 eV ≈ 23.06 kcal/mol. |
| **Bohr** | Atomic unit of length; 1 Bohr ≈ 0.529 Å. |
| **Angstrom (Å)** | 10⁻¹⁰ m; standard unit for interatomic distances. |
| **cm⁻¹** | Reciprocal centimeters (wavenumber); the standard unit for vibrational frequencies. Imaginary frequencies appear as negative values. |
| **Imaginary Frequency** | A vibrational frequency corresponding to a negative eigenvalue of the Hessian. A TS has exactly one (first-order saddle point). Reported as a negative cm⁻¹ value. |

(frequency-thresholds)=
### Imaginary-mode criterion and QRRHO rotor cutoff

Imaginary-mode classification and the QRRHO rotor cutoff serve different purposes:

| Threshold | Role | Source |
|-----------|------|--------|
| **ν < −5.00 cm⁻¹** | Default imaginary-mode criterion. | Configurable: `freq.zero_cutoff_cm`. |
| **100 cm⁻¹** | *QRRHO rotor cutoff* (Grimme). Positive low-frequency vibrations are damped between harmonic-oscillator and free-rotor entropy in `freq` thermochemistry; it changes only entropy / Gibbs free energy. | Fixed (not configurable in mlmm-toolkit). |

## CLI Conventions

Boolean options, residue selectors, and atom selectors are described in [Common options and selectors](cli-conventions.md).

## Notes

* **Imaginary modes and negative signs**: n_imag counts only the modes below the imaginary-mode criterion. Every frequency is still reported with its sign, and the count of all negative values (`n_negative_modes` in `result.json`) is a separate diagnostic that does not change convergence.

## See Also

- [Getting Started](getting-started.md) — the shortest run and which page to read next
- [Installation](installation.md) — setup and dependencies
- [all](all.md) — how pocket extraction, MEP search, and post-processing fit together
- [Troubleshooting](troubleshooting.md) — common errors and fixes
- [YAML Reference](yaml-reference.md) — configuration file format
- [MLIP Backends](backends.md) — machine learning potential details
- [ML/MM Calculator](mlmm-calc.md) — ONIOM coupling, link atoms, and the MM Hessian
