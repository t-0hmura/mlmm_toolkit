# ML/MM Calculator

## Overview

> **Summary:** ONIOM-like ML/MM calculator for PySisyphus, coupling FAIR-Chem UMA (high-level ML) and hessian_ff (low-level MM) to compute energies, forces, and Hessians for enzyme active-site models.

`mlmm_calc.mlmm` implements a subtractive ONIOM-style ML/MM calculator that combines a machine-learning interatomic potential (FAIR-Chem UMA) with a molecular-mechanics force field (Amber prmtop-based `hessian_ff`). It serves as the core calculator for all ML/MM optimization, path search, scan, frequency, and IRC workflows in `mlmm_toolkit`.

The calculator operates without link atoms. The ML region is defined by a model PDB (`model.pdb`), the MM topology comes from an Amber prmtop (`real.parm7`), and coordinates are taken from the input PDB (`input.pdb`). An internal `real.rst7` is generated via ParmEd by combining `real.parm7` with coordinates from `input.pdb` -- no external `real.rst7` or `real.pdb` is required.

## Three-layer scheme (energy / force / Hessian)

The calculator combines three evaluations using the ONIOM subtraction:

| Layer | System | Method | Description |
| --- | --- | --- | --- |
| **REAL-low** | Full system | MM (hessian_ff) | Full system evaluated with Amber prmtop-based MM |
| **MODEL-low** | ML subset | MM (hessian_ff) | ML region evaluated with MM |
| **MODEL-high** | ML subset + link-H | ML (UMA) | ML region evaluated with FAIR-Chem UMA |

The combined energy is:

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

Forces and Hessians follow the same subtraction pattern.

## Layering for Hessian / optimization

The implementation uses 3-layer B-factor encoding and optional Hessian-target MM selection:

- **ML region** (B-factor = 0.0): Treated with UMA machine-learning potential
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

- **`"hessian_ff"`** (default): Analytical Hessian via `hessian_ff` (active atoms only, then optionally expanded to full Cartesian shape with frozen rows/cols zero-filled). CPU-only backend.
- **`"openmm"`**: Finite-difference (FD) Hessian via OpenMM. Supports both CPU and CUDA platforms. Useful for force fields not supported by `hessian_ff` or when OpenMM is already in your workflow.

**Example YAML configuration:**
```yaml
mlmm:
  mm_backend: openmm  # Use OpenMM for MM calculations
  mm_device: cuda     # Use CUDA (or "cpu")
```

### UMA Hessian modes
- `"Analytical"`: Second-order autograd on the selected device.
- `"FiniteDifference"`: Central differences of forces.

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
| Forces | eV/A | Hartree/Bohr |
| Hessian | eV/A^2 | Hartree/Bohr^2 |

The PySisyphus interface returns values converted to atomic units (Hartree/Bohr).

## Implementation notes
- To reduce memory usage, Hessian assembly uses in-place `add_`/`sub_`/`mul_`.
- Atom ID mapping between REAL/MODEL is based solely on `input.pdb`.
- `real.rst7` is created internally from `input.pdb` coordinates using ParmEd.

---

## See Also

- [opt](opt.md) -- Single-structure geometry optimization using the ML/MM calculator
- [tsopt](tsopt.md) -- Transition state optimization
- [freq](freq.md) -- Vibrational frequency analysis
- [YAML Reference](yaml-reference.md) -- `calc`/`mlmm` configuration keys
