# Glossary

This page provides definitions for abbreviations and technical terms used throughout the mlmm_toolkit documentation.

---

## ML/MM & ONIOM

| Term | Full Name | Description |
|------|-----------|-------------|
| **ML/MM** | Machine Learning / Molecular Mechanics | A multi-scale method that couples a machine-learning interatomic potential (for the reactive region) with a classical force field (for the surrounding environment). Analogous to QM/MM but with ML replacing QM. |
| **ONIOM** | Our own N-layered Integrated molecular Orbital and molecular Mechanics | A multi-layer energy decomposition scheme. mlmm_toolkit uses an ONIOM-like subtraction: E_total = E_REAL_low + E_MODEL_high - E_MODEL_low. |
| **Real system** | — | The full set of atoms (all 3 layers). Evaluated at the MM (low) level in the ONIOM decomposition. |
| **Model system** | — | The ML region (Layer 1). Evaluated at both the UMA (high) and MM (low) levels in the ONIOM decomposition. |
| **hessian_ff** | — | A C++ native extension that evaluates Amber force field energies, forces, and analytical Hessians. Used as the MM engine in mlmm_toolkit. |
| **3-layer system** | — | mlmm_toolkit's B-factor partitioning scheme: ML (B=0.0), Movable-MM (B=10.0), Frozen (B=20.0). |
| **B-factor encoding** | — | Convention of storing layer membership in the PDB B-factor (temperature factor) column: 0.0 = ML, 10.0 = Movable-MM, 20.0 = Frozen. Hessian-target MM is controlled by cutoffs/explicit indices. |

---

## Amber & Force Field

| Term | Full Name | Description |
|------|-----------|-------------|
| **parm7** | Amber Parameter/Topology File | A file containing atom types, partial charges, bonding connectivity, and force field parameters for an Amber system. Also called.prmtop. |
| **rst7** | Amber Restart File | A file containing atomic coordinates (and optionally velocities and box dimensions) for an Amber system. Also called.inpcrd. |
| **AmberTools** | — | A free suite of tools for molecular dynamics preparation, including tleap, antechamber, and parmchk2. Required by `mlmm mm-parm`. |
| **tleap** | — | An AmberTools program that builds Amber topology/coordinate files from PDB structures and force field libraries. |
| **antechamber** | — | An AmberTools program that assigns GAFF2 atom types and AM1-BCC partial charges to small molecules. |
| **parmchk2** | — | An AmberTools program that checks and supplies missing force field parameters for GAFF2 typing. |
| **GAFF2** | General Amber Force Field 2 | A general-purpose force field for small organic molecules, used to parameterize non-standard residues (substrates, cofactors). |
| **ff19SB** | — | An Amber protein force field used for standard amino acid residues. |
| **AM1-BCC** | — | A charge model that combines AM1 (semi-empirical) Mulliken charges with bond charge corrections (BCC) to approximate HF/6-31G* RESP charges. |

---

## Reaction Path & Optimization

| Term | Full Name | Description |
|------|-----------|-------------|
| **MEP** | Minimum Energy Path | The lowest-energy pathway connecting reactants to products through a transition state (on a potential energy surface). |
| **TS** | Transition State | A first-order saddle point on the potential energy surface, typically the highest-energy point along the reaction coordinate. |
| **IRC** | Intrinsic Reaction Coordinate | A mass-weighted steepest-descent path from a TS toward reactants and products. Often used to validate TS connectivity. |
| **GSM** | Growing String Method | A string-based method that grows images from endpoints and optimizes them to approximate an MEP. |
| **DMF** | Direct Max Flux | A chain-of-states method for optimizing an MEP by maximizing flux along the pathway. In mlmm it is selected with `--mep-mode dmf`. |
| **NEB** | Nudged Elastic Band | A chain-of-states method that uses spring forces to maintain image spacing along a reaction path. |
| **HEI** | Highest-Energy Image | The image along an MEP with maximum energy; often used as a TS guess. |
| **Image** | — | A single geometry (one "node") along a chain-of-states path. |
| **Segment** | — | An MEP between two adjacent endpoints (e.g., R -> I1, I1 -> I2,...). |

---

## Optimization Algorithms

| Term | Full Name | Description |
|------|-----------|-------------|
| **L-BFGS** | Limited-memory BFGS | A quasi-Newton optimization algorithm that approximates the Hessian using a limited history of gradients. Used in `--opt-mode light`. |
| **RFO** | Rational Function Optimization | A trust-region optimization method that uses explicit Hessian information. Used in `--opt-mode heavy`. |
| **RS-I-RFO** | Restricted-Step Image-RFO | A variant of RFO for saddle point (TS) optimization that follows one negative eigenvalue. |
| **Dimer** | Dimer Method | A TS optimization method that estimates the lowest curvature mode without computing the full Hessian. Used in `--opt-mode light` for TSOPT. |
| **PHVA** | Partial Hessian Vibrational Analysis | Computing vibrational frequencies using only the Hessian block for active (non-frozen) atoms. Default in `freq`. |

---

## Machine Learning & Calculators

| Term | Full Name | Description |
|------|-----------|-------------|
| **MLIP** | Machine Learning Interatomic Potential | A model (often neural-network-based) that predicts energies and forces from atomic structures, trained on quantum-mechanical data. |
| **UMA** | Universal Machine-learning potential for Atoms | Meta's family of pretrained MLIPs used as the ML calculator backend in mlmm_toolkit. |
| **Analytical Hessian** | — | Computing the exact second derivatives of energy; faster but requires more VRAM. |
| **Finite Difference** | — | Approximating derivatives by small displacements; slower but more memory-efficient. |

---

## Quantum Chemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **QM** | Quantum Mechanics | First-principles electronic structure calculations (DFT, HF, post-HF, etc.). |
| **QM/MM** | Quantum Mechanics / Molecular Mechanics | A multi-scale method coupling QM for the reactive region with MM for the environment. mlmm_toolkit replaces the QM layer with ML (UMA). |
| **DFT** | Density Functional Theory | A quantum-mechanical method that models electronic structure via electron density functionals. |
| **Hessian** | — | The matrix of second derivatives of energy with respect to atomic coordinates; used for vibrational analysis and TS optimization. |
| **SP** | Single Point | A calculation at a fixed geometry (no optimization); often used for higher-level energy refinement. |
| **Spin Multiplicity** | — | 2S+1, where S is total spin. Singlet = 1, doublet = 2, triplet = 3, etc. |

---

## Structural Biology & Pocket Extraction

| Term | Full Name | Description |
|------|-----------|-------------|
| **PDB** | Protein Data Bank | A file format and database for macromolecular 3D structures. |
| **XYZ** | — | A simple text format listing atomic symbols and Cartesian coordinates. |
| **GJF** | Gaussian Job File | An input format for Gaussian; mlmm reads charge/multiplicity and coordinates from these files. |
| **Pocket** | Active-site Pocket | A truncated structure around the substrate(s) used to reduce system size for MEP/TS search. Also called "cluster model". |
| **Cluster Model** | — | Synonym for pocket; a computationally tractable subset of the full enzyme-substrate complex. |
| **Link Hydrogen** | — | A hydrogen atom added to cap severed bonds when extracting a pocket from a larger structure. |
| **Backbone** | — | The main chain of a protein (N-C_alpha-C-O atoms). Can be excluded during pocket extraction with `--exclude-backbone`. |

---

## Thermochemistry

| Term | Full Name | Description |
|------|-----------|-------------|
| **ZPE** | Zero-Point Energy | The vibrational energy at 0 K; a quantum correction to the electronic energy. |
| **Gibbs Energy** | Free Energy (G) | G = H - TS; includes thermal and entropic contributions. |
| **Enthalpy** | (H) | H = E + PV; total heat content at constant pressure. |
| **Entropy** | (S) | A measure of disorder; contributes -TS to Gibbs energy. |

---

## Units & Constants

| Term | Description |
|------|-------------|
| **Hartree** | Atomic unit of energy; 1 Hartree ~ 627.5 kcal/mol ~ 27.21 eV. |
| **kcal/mol** | Kilocalories per mole; a common unit for reaction energetics. |
| **kJ/mol** | Kilojoules per mole; 1 kcal/mol ~ 4.184 kJ/mol. |
| **eV** | Electron volt; 1 eV ~ 23.06 kcal/mol. |
| **Bohr** | Atomic unit of length; 1 Bohr ~ 0.529 Angstrom. |
| **Angstrom** | 10^-10 meters; standard unit for atomic distances. |

---

## CLI Conventions

| Term | Description |
|------|-------------|
| **Boolean option** | CLI flags that take `True` or `False` (capitalized). Example: `--tsopt`. |
| **Residue selector** | A specification like `'SAM,GPP'` (names) or `'A:123,B:456'` (chain:ID). |
| **Atom selector** | A specification like `'TYR,285,CA'` identifying a specific atom by residue name, number, and atom name. |

---

## See Also

- [Getting Started](getting_started.md) -- installation and a first run
- [Concepts & Workflow](concepts.md) -- how pocket extraction, ML/MM layers, MEP search, and post-processing fit together
- [Common Error Recipes](recipes_common_errors.md) -- symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- common errors and fixes
- [YAML Reference](yaml_reference.md) -- configuration file format
- [ML/MM Calculator](mlmm_calc.md) -- machine learning potential details
