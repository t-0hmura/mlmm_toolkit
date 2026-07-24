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

### Electronic embedding

Electronic embedding is unavailable in v0.3.3. `--embedcharge`,
`--embedcharge-cutoff`, and `calc.embedcharge: true` are retained only so older
commands fail with an explicit diagnostic before calculation. The previous
experimental correction added an electronic ML--MM interaction on top of the
Amber interaction already retained by the subtractive ONIOM expression and
used an uncapped model inconsistent with the link-H high-level system. Use the
default mechanical embedding (`--no-embedcharge`) and do not reuse results
generated with the experimental path.

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
  - CMAP torsion corrections (implemented but disabled by default, as in Gaussian)
- **`"openmm"`**: Finite-difference (FD) Hessian via OpenMM. Supports both CPU and CUDA platforms. Covers force fields not supported by `hessian_ff`, or cases where OpenMM is already in your workflow. See [Device Configuration & HPC Setup](device-hpc.md) for mm_backend/mm_device YAML examples and VRAM trade-offs.

ML Hessian: `Analytical` (backend autograd/native Hessian for UMA, ORB, MACE, and AIMNet2) or `FiniteDifference` (central differences of forces for any backend). An explicit analytical request fails if the installed backend lacks its required API; it is never silently downgraded. See [YAML Reference](yaml-reference.md) for VRAM guidance.

### CMAP in the model system

CMAP (Cross-Map backbone dihedral correction) is a 5-atom torsion correction term used in protein AMBER force fields to improve backbone conformational sampling. In ONIOM, the model system parm7 is generated by slicing the real topology to the ML region.

By default (`use_cmap: false`), CMAP terms are **excluded** from the model parm7:

| Region | E_MM(real) | E_MM(model) | ONIOM net effect |
|--------|-----------|------------|-----------------|
| `use_cmap: false` (default) | CMAP included | CMAP **excluded** | Model backbone CMAP remains in E_total |
| `use_cmap: true` | CMAP included | CMAP included | Model backbone CMAP cancels in subtraction |

This default behavior is consistent with Gaussian ONIOM, which also omits CMAP from model MM parameters. For typical active-site models (ligand + functional residues, no backbone atoms in ML region), CMAP in the model is zero in either case.

**Example YAML configuration:**
```yaml
mlmm:
 use_cmap: true  # Enable CMAP in model parm7 (non-Gaussian-compatible behavior)
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
