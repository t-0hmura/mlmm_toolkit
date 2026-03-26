# Concepts & Workflow

This page explains the key terms in mlmm-toolkit -- the ML/MM 3-layer system, ONIOM decomposition, segments, images, and templates -- and how the `all` command ties together the subcommands.

---

## Workflow at a glance

Most workflows follow this flow:

```text
Full system(s) (PDB/XYZ)
 |
 +- ML region definition [extract] <- requires PDB when you use --center/-c
 | |
 | ML region structure(s) (PDB)
 | |
 | +- Amber topology [mm-parm] <- generates parm7/rst7 via AmberTools
 | | |
 | +- 3-layer assignment [define-layer] <- B-factor encoding for ML/MM layers
 | | |
 | | v
 | | ML/MM system (PDB with B-factors)
 | | |
 | +- (optional) staged scan [scan] <- single-structure workflows
 | | |
 | | v
 | | Ordered intermediates
 | |
 | +- MEP search [path-search] or [path-opt]
 | |
 | MEP trajectory (mep_trj.xyz) + energy diagrams
 |
 +- (optional) TS optimization + IRC [tsopt] -> [irc]
 +- (optional) thermochemical analysis [freq] to obtain ΔG
 +- (optional) single-point DFT [dft]
```

Each stage is available as an individual subcommand. The `mlmm all` command runs many stages end-to-end.

```{important}
Transition states: treat HEI / [`tsopt`](tsopt.md) outputs as **TS candidates** until validated via [`freq`](freq.md) (a single imaginary mode) and [`irc`](irc.md) (endpoints reach intended minima).
```

---

## ML/MM 3-layer system

A central concept in mlmm-toolkit is the **3-layer ML/MM partitioning** of the system. Each atom belongs to one of three layers, encoded in the PDB B-factor column:

| Layer | B-factor | Description |
|-------|----------|-------------|
| **ML** (Layer 1) | 0.0 | The reactive region. Full MLIP energy, forces, and Hessian. |
| **Movable-MM** (Layer 2) | 10.0 | MM atoms allowed to move during optimization. |
| **Frozen** (Layer 3) | 20.0 | Coordinates are fixed; no optimization. |

B-factor values are encoded in PDB columns 61-66 (the temperature factor column).

The [`define-layer`](define_layer.md) subcommand assigns these B-factors based on distance from the ML region:

- Atoms/residues within `--radius-freeze` (default 8.0 Å) are assigned to Movable-MM.
- Atoms/residues beyond `--radius-freeze` are Frozen.

Hessian-target MM atoms are controlled by calculator options (`hess_cutoff`, explicit `hess_mm_atoms`, etc.), not by a dedicated B-factor layer.

```{tip}
The B-factor encoding allows you to visually inspect layer assignments in any molecular viewer that can color by B-factor.
```

---

## ONIOM-like energy decomposition

mlmm-toolkit uses an **ONIOM-like** scheme to combine ML and MM energies:

```
E_total = E_REAL_low + E_MODEL_high - E_MODEL_low
```

where:
- **REAL** = the full system (all atoms)
- **MODEL** = the ML region (subset of atoms)
- **high** = the selected MLIP backend (default: UMA; ORB, MACE, AIMNet2 also supported)
- **low** = hessian_ff (Amber-based classical force field)

This means:
1. The **full system** is evaluated at the MM level (hessian_ff).
2. The **ML region** is evaluated at both the MLIP level and the MM level.
3. The MM contribution of the ML region is subtracted to avoid double-counting.

The same decomposition applies to forces and (where applicable) Hessians. Link-hydrogen contributions are redistributed to the ML and MM host atoms via a Jacobian.

The MLIP backend is selected via `-b/--backend` (default: `uma`). Alternative backends (`orb`, `mace`, `aimnet2`) are installed as optional dependencies (e.g., `pip install "mlmm-toolkit[orb]"`).

When `--embedcharge` is enabled, an xTB point-charge embedding correction is applied to account for the electrostatic influence of the MM environment on the ML region.

---

## Link atoms

When the ML/MM boundary cuts a covalent bond, a **link hydrogen atom** is inserted to cap the dangling bond in the model (ML) system. Two placement methods are available via `--link-atom-method`:

| Method | Placement | Recommended |
|--------|-----------|-------------|
| **scaled** (g-factor, default) | `r_L = r_QM + g·(r_MM − r_QM)` where `g = (CR_QM + CR_H)/(CR_QM + CR_MM)` | Yes (smooth PES, constant Jacobian) |
| **fixed** | `r_L = r_QM + d·û` where `d` = 1.09 Å (C) / 1.01 Å (N), `û` = unit vector toward MM | No (geometry-dependent Jacobian) |

The default is **scaled** (Morokuma–Dapprich g-factor), the same method used in Gaussian ONIOM. The link atom position scales linearly with the QM–MM distance, producing a smooth PES and a constant Jacobian (no second-derivative correction needed). The **fixed** method places the link hydrogen at a fixed distance along the bond axis regardless of QM–MM distance.

Forces on the link atom are redistributed back to the QM and MM host atoms via the Jacobian:

```
F_QM += (1−g) · F_link    (scaled)
F_MM += g · F_link
```

The same transformation applies to the Hessian: `H_redistributed = Jᵀ H_link J`.

---

## Microiteration

For systems with many movable MM atoms, simultaneous optimization of all coordinates is expensive because the high-level (MLIP) gradient must be evaluated at every step — even when only the MM environment is relaxing.

**Microiteration** (Gaussian 16-style) splits the optimization into alternating macro and micro steps:

```
repeat until converged:
    MACRO step  — 1 RFO step on ML atoms + link-atom MM parents (full ONIOM force)
    MICRO step  — L-BFGS relaxation of remaining MM atoms (MM-only force)
```

| | Macro step | Micro step |
|---|---|---|
| **Calculator** | Full ONIOM (`E_MM_real + E_ML − E_MM_model`) | MM force field only (`E_MM_real`) |
| **Coordinates optimized** | ML atoms + link-atom MM parents | Movable MM (excluding link-atom MM parents) |
| **Optimizer** | RFO (explicit Hessian, BFGS-updated) | L-BFGS (Hessian-free, from scratch each cycle) |
| **Convergence** | `--thresh` (default: `gau`) | `--micro-thresh` (default: same as `--thresh`) |

```{note}
**Why link-atom MM parents move in the macro step:**
Link-atom MM parent atoms are included in the macro optimization set to maintain consistency at the ML/MM boundary. When using the scaled (g-factor) link atom method, `r_L = (1−g)·r_QM + g·r_MM` couples the link atom position to **both** the QM and MM parents. If the MM parent moved during the micro step (under MM-only forces with no ML contribution), the link atom position would shift between cycles, creating energy oscillation in the macro step. Freezing the MM parents during the micro step and moving them with the ML atoms in the macro step eliminates this coupling mismatch.
```

### Convergence thresholds

pysisyphus provides several preset thresholds (units: Hartree/Bohr for forces, Bohr for steps):

| Preset | max(force) | rms(force) | max(step) | rms(step) |
|--------|-----------|-----------|----------|----------|
| `gau_loose` | 2.5×10⁻³ | 1.7×10⁻³ | 1.0×10⁻² | 6.7×10⁻³ |
| `gau` | 4.5×10⁻⁴ | 3.0×10⁻⁴ | 1.8×10⁻³ | 1.2×10⁻³ |
| `gau_tight` | 1.5×10⁻⁵ | 1.0×10⁻⁵ | 6.0×10⁻⁵ | 4.0×10⁻⁵ |
| `baker` | 3.0×10⁻⁴ | 2.0×10⁻⁴ | 3.0×10⁻⁴ | 2.0×10⁻⁴ |

`overachieve_factor` is a convergence shortcut: when `max(force)` and `rms(force)` are both below `threshold / overachieve_factor`, convergence is declared even if the step-size criteria have not yet been met. In the default configuration, `overachieve_factor` is set to **0.0** (disabled) for the main optimizer (opt/tsopt macro steps). It is only active inside the **microiteration** MM relaxation loop (`LayerOpt`, where it is set to 3). Users can enable it for the main optimizer via YAML (`overachieve_factor: 3`).

Enable microiteration with `--microiter` (default for `--opt-mode hess`):

```bash
mlmm opt -i layered.pdb --parm system.parm7 -q 0 --opt-mode hess --microiter
mlmm opt -i layered.pdb --parm system.parm7 -q 0 --opt-mode hess --no-microiter  # disable
```

---

## hessian_ff: the MM engine

`hessian_ff` is a C++ native extension that evaluates Amber force field energies, forces, and especially **analytical Hessians**. It reads Amber parm7/rst7 topology files and supports:

- Bond, angle, dihedral, and improper terms
- Van der Waals (Lennard-Jones) interactions
- Electrostatic interactions
- Analytical second derivatives (Hessian)
- CPU execution (GPU memory is reserved for MLIP inference)
- CMAP torsion corrections (implemented but disabled by default, as in Gaussian)

Unlike OpenMM, `hessian_ff` is designed specifically for providing the **MM Hessian** needed by the ONIOM-like coupling and vibrational analysis.

---

## Amber parm7/rst7 topology

The internal MM calculation requires Amber topology files:

- **parm7** (parameter/topology file): Contains atom types, charges, bonding connectivity, and force field parameters.
- **rst7** (restart/coordinate file): Contains atomic coordinates. (Not directly specified by users; mlmm reads coordinates from the input PDB.)

These are generated by the [`mm-parm`](mm_parm.md) subcommand using **AmberTools** (tleap, antechamber, parmchk2). The command automatically:

- Identifies non-standard residues (substrates, cofactors)
- Parameterizes them with **GAFF2** (General Amber Force Field 2)
- Assigns AM1-BCC partial charges
- Builds the full topology with ff19SB for protein residues

---

## Key objects and terms

### Full system vs. ML region
- **Full system**: your original structure(s). In enzyme use-cases this is typically a protein-ligand complex.
- **ML region**: the reactive region treated with the MLIP backend. Defined by the `extract` subcommand based on distance from specified substrates.

ML region definition is controlled by:
- `-c/--center`: how to locate the substrate (residue IDs, residue names, or a substrate-only PDB).
- `-r/--radius`, `--radius-het2het`, `--include-H2O`, `--exclude-backbone`, `--add-linkH`, `--selected-resn`.

### Real system vs. Model system (ONIOM terminology)
- **Real system**: the entire set of atoms (all 3 layers). Evaluated at the MM (low) level.
- **Model system**: the ML region (Layer 1 only). Evaluated at both the MLIP (high) and MM (low) levels.

### Images and segments
- **Image**: a single geometry (one "node") along a chain-of-states path in Minimum Energy Path (MEP) search methods such as the Growing String Method (GSM).
- **Segment**: an MEP between two adjacent endpoints (e.g., R -> I1, I1 -> I2,...). A multi-structure run is decomposed into segments.

### Templates and file conversion (`--convert-files`)
`mlmm-toolkit` often writes a **trajectory** (e.g., `mep_trj.xyz`, `irc_trj.xyz`). When you supply topology-aware inputs (PDB templates), it can optionally write companion files:
- `.pdb` companions when a PDB template exists

This behavior is controlled globally by `--convert-files/--no-convert-files` (default: `True`).

---

## Three common workflow modes

### 1) Multi-structure MEP search (R ->... -> P)
Use this when you already have **two or more** full structures along a reaction coordinate.

Typical command:

```bash
mlmm -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

### 2) Single-structure staged scan -> MEP
Use this when you prefer to define reaction coordinates yourself, rather than providing multiple endpoint structures.

Typical command:

```bash
mlmm -i holo.pdb -c '308,309' \
 --scan-lists '[("TYR,285,CA","MMT,309,C10",2.20)]'
```

### 3) TSOPT-only mode (ML/MM TS optimization)
Use this when you already have a TS candidate (or want a quick TS optimization on one structure).

Typical command:

```bash
mlmm -i ts_guess.pdb -c 'SAM,GPP' --tsopt
```

---

## When to use `all` vs individual subcommands

### Prefer `mlmm all` when...
- You want an **end-to-end** run (ML/MM model setup -> MEP search -> TSOPT/IRC -> freq/DFT).
- You are still exploring the workflow and want a single command to manage outputs.

### Prefer subcommands when...
- You want to run each stage step by step, verifying results at each point. For complex reactions, this approach often works better than running everything at once.
- You want to mix-and-match a custom workflow (e.g., your own endpoint preparation).
- You already have parm7/rst7 and layer-assigned PDB files from a previous run.
- You want to generate Gaussian/ORCA ONIOM input files via `oniom-export --mode g16|orca`.

---

## A few CLI conventions worth knowing

```{important}
- Boolean options accept both `--flag` / `--no-flag` and value style `--flag True/False` (`yes/no`, `1/0` are also accepted). Prefer toggle style.
- With multiple PDB inputs, all files should have the **same atoms in the same order** (only coordinates differ).
- For enzyme use-cases, you usually want hydrogens present in the input PDB.
- Most subcommands require `--parm` and `--model-pdb` for ML/MM calculations.
```

---

## Next steps

### Getting started
- [Getting Started](getting_started.md) -- installation and first run
- [Common Error Recipes](recipes_common_errors.md) -- symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- common errors and fixes

### Core subcommands
| Subcommand | Purpose | Documentation |
|------------|---------|---------------|
| `all` | End-to-end workflow | [all.md](all.md) |
| `extract` | ML region definition | [extract.md](extract.md) |
| `mm-parm` | Amber topology generation | [mm_parm.md](mm_parm.md) |
| `define-layer` | 3-layer assignment | [define_layer.md](define_layer.md) |
| `path-search` | Recursive MEP search | [path_search.md](path_search.md) |
| `tsopt` | TS optimization | [tsopt.md](tsopt.md) |
| `freq` | Vibrational analysis | [freq.md](freq.md) |
| `dft` | Single-point DFT | [dft.md](dft.md) |
| `oniom-export` | Gaussian ONIOM / ORCA QM/MM input generation (`--mode g16\|orca`) | [oniom_export.md](oniom_export.md) |

### Reference
- [YAML Reference](yaml_reference.md) -- complete YAML configuration options
- [Glossary](glossary.md) -- terminology reference
