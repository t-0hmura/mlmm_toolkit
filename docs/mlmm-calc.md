#  ML/MM Calculator

## Overview

ONIOM-like ML/MM calculator for PySisyphus, coupling an MLIP backend (high-level ML) and hessian_ff (low-level MM) to compute energies, forces, and especially **analytical Hessians** for enzyme active-site models. Multiple MLIP backends are supported via `-b/--backend`.

`mlmm_calc.mlmm` implements a subtractive ONIOM-style ML/MM calculator that combines a machine learning interatomic potential (MLIP) with a molecular-mechanics force field (Amber prmtop-based `hessian_ff`). It serves as the core calculator for all ML/MM optimization, path search, scan, frequency, and IRC workflows in `mlmm`.

### Multi-backend architecture

The ML (high-level) component is provided by one of several MLIP backends, selected via the `-b/--backend` CLI option or the `mlmm.backend` YAML key:

| Backend | Value | Package | Install |
| --- | --- | --- | --- |
| FAIR-Chem UMA | `uma` (default) | `fairchem-core` | `pip install mlmm-toolkit` |
| ORB | `orb` | `orb-models` | `pip install "mlmm-toolkit[orb]"` |
| MACE | `mace` | `mace-torch` | dedicated env: `pip uninstall -y fairchem-core && pip install mace-torch` |
| AIMNet2 | `aimnet2` | `aimnet` | `pip install "mlmm-toolkit[aimnet]"` |

See [MLIP Backends](backends.md) for per-backend kwargs, model identifiers, precision options, and how to add a backend.

The calculator automatically generates link hydrogen atoms at covalent ML/MM boundaries. The ML region is defined by a model PDB (`model.pdb`), the MM topology comes from an Amber prmtop (`real.parm7`), and coordinates are taken from the input PDB (`input.pdb`). An internal `real.rst7` is generated via ParmEd by combining `real.parm7` with coordinates from `input.pdb` -- no external `real.rst7` or `real.pdb` is required.

## Three-layer scheme (energy / force / Hessian)

The calculator combines three evaluations using the ONIOM subtraction:

| Layer | System | Method | Description |
| --- | --- | --- | --- |
| **REAL-low** | Full system | MM (hessian_ff) | Full system evaluated with Amber prmtop-based MM |
| **MODEL-low** | ML subset | MM (hessian_ff) | ML region evaluated with MM |
| **MODEL-high** | ML subset + link-H | ML (MLIP) | ML region evaluated with the selected MLIP backend (default: UMA) |

The combined energy is:

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

Forces and Hessians follow the same subtraction pattern.

## Layering for Hessian / optimization

The implementation uses 3-layer B-factor encoding and optional Hessian-target MM selection:

- **ML region** (B-factor = 0.0): Treated with the selected MLIP backend (default: UMA)
- **Movable-MM** (B-factor = 10.0): MM atoms that move during optimization
- **Frozen** (B-factor = 20.0): Fixed MM atoms
- **Hessian-target MM** (not a dedicated B-factor): selected by `hess_cutoff` and/or explicit `hess_mm_atoms`

Layer assignment is controlled by `hess_cutoff`, `movable_cutoff`, `use_bfactor_layers`, and explicit `*_mm_atoms` lists.

## Features

### Link-atom redistribution
Forces and Hessian contributions from link atoms are redistributed to the ML/MM parent atoms via a Jacobian. The redistribution adds:
- The self term `J^T H J`
- The geometry-dependent second term `sum (dJ^T/dx * f_L)` in-place to the parent atoms

### MM Hessian

The MM backend can be selected via the `mm_backend` parameter:

- **`"hessian_ff"`** (default backend): CPU-only MM engine with an analytical-Hessian capability. The effective default remains finite difference (`mm_fd: true`); set `mm_fd: false` to use its analytical Hessian. Active blocks can optionally be expanded to full Cartesian shape with frozen rows/columns zero-filled.
  - CMAP torsion corrections, preserved when present in the parm7
- **`"openmm"`**: Finite-difference (FD) Hessian via OpenMM. Supports both CPU and CUDA platforms. Covers force fields not supported by `hessian_ff`, or cases where OpenMM is already in your workflow. See [Device Configuration & HPC Setup](device-hpc.md) for mm_backend/mm_device YAML examples and VRAM trade-offs.

ML Hessian: `Analytical` (backend autograd/native Hessian for UMA, ORB, MACE, and AIMNet2) or `FiniteDifference` (central differences of forces for any backend). An explicit analytical request fails if the installed backend lacks its required API; it is never silently downgraded. See [YAML Reference](yaml-reference.md) for VRAM guidance.

### CMAP in the two MM layers

CMAP (Cross-Map backbone dihedral correction) is a 5-atom torsion correction term used by force fields such as ff19SB. The REAL and MODEL MM calculations must use the same CMAP policy in the subtractive expression.

| Region | E_MM(real) | E_MM(model) | ONIOM net effect |
|--------|-----------|------------|-----------------|
| `use_cmap: true` (default) | CMAP included when present | CMAP included when present | Complete model-internal CMAP cancels; boundary terms remain in the low-level coupling |
| `use_cmap: false` | CMAP excluded | CMAP excluded | Explicit modified-force-field calculation without CMAP |

For ff19SB, CMAP replaces the corresponding zeroed backbone cosine terms ([Tian et al., 2020](https://doi.org/10.1021/acs.jctc.9b00591)). Preserving it is therefore the force-field-faithful default. `use_cmap: false` removes CMAP from both MM layers; it is not an ff19SB-compatible default.

**Example YAML configuration:**
```yaml
mlmm:
 use_cmap: false  # Explicitly remove CMAP from both MM layers
```

## Inputs

| Input | Description |
| --- | --- |
| `input.pdb` | Input structure (residue/atom names are read from here) |
| `real.parm7` | Amber prmtop (topology of the full REAL system) |
| `model.pdb` | PDB defining the ML region (used to determine atom IDs) |

## Units

| Quantity | Internal unit | PySisyphus interface |
| --- | --- | --- |
| Energy | eV | Hartree |
| Forces | eV/Å | Hartree/Bohr |
| Hessian | eV/Å² | Hartree/Bohr² |

The PySisyphus interface returns values converted to atomic units (Hartree/Bohr).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [opt](opt.md) — Single-structure geometry optimization using the ML/MM calculator
- [tsopt](tsopt.md) — Transition state optimization
- [freq](freq.md) — Vibrational frequency analysis
- [YAML Reference](yaml-reference.md) — `calc`/`mlmm` configuration keys
- [MLIP Backends](backends.md) — Backend selection, install, precision, and the add-a-backend recipe (canonical backend reference)
- [Device Configuration & HPC Setup](device-hpc.md) — ML/MM device settings and HPC submission
