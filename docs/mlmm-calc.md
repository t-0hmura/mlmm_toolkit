# ML/MM Calculator

This page explains how the ML/MM calculator works: it computes the ML region with an MLIP, computes the surrounding protein and solvent with an Amber MM force field, and combines the two by ONIOM subtraction. It also shows how to call the calculator from Python.

All ML/MM optimization, path search, scan, frequency, and IRC commands use this calculator. The ML region is computed by the backend selected with `-b/--backend`: `uma` (default), `orb`, `mace`, `aimnet2`, or `dft` for DFT/MM. Installation, model names, and options for each backend are in [MLIP Backends](backends.md).

## ONIOM energy decomposition

The calculator combines three evaluations by ONIOM subtraction:

| Evaluation | System | Method | Description |
| --- | --- | --- | --- |
| **REAL-low** | Full system | MM (`hessian_ff` or OpenMM) | Full system evaluated with the Amber parm7 force field |
| **MODEL-low** | ML region | MM (same engine) | ML region evaluated with MM |
| **MODEL-high** | ML region + link H | MLIP or DFT | ML region evaluated with the selected backend (default: UMA) |

The combined energy is:

```
E_ONIOM = E(REAL-low) - E(MODEL-low) + E(MODEL-high)
```

The full system is evaluated with MM, the ML region is evaluated at both the high level and the MM level, and the MM energy of the ML region is subtracted so that it is not counted twice. Forces and Hessians follow the same subtraction.

### Comparison with conventional QM/MM

| Aspect | Conventional QM/MM | mlmm-toolkit ML/MM |
| --- | --- | --- |
| High-level method | DFT, HF, post-HF | MLIP (UMA, ORB, MACE, AIMNet2) or DFT (`-b dft`) |
| Low-level method | OpenMM / Amber | `hessian_ff` (default) / OpenMM |
| Link atoms | Usually required | Added automatically for every parm7 bond that crosses the ML/MM boundary |
| Embedding | Electrostatic embedding is common | Mechanical by default; `--embedcharge` adds the MM point charges (an xTB correction for MLIP backends, the PySCF Hamiltonian for DFT) |
| Speed | Slow (QM is the bottleneck) | Fast with an MLIP (GPU inference); `-b dft` runs at DFT cost |

## Layering for Hessian / optimization

Each atom belongs to one of three layers (ML, movable MM, or frozen MM), read from the B-factors of the input PDB ({ref}`Building the ML region and layers › The MM layers <mm-layers>`). Two further settings decide which MM atoms enter the Hessian and which MM atoms move:

- **Hessian-target MM** (no B-factor of its own): the movable MM atoms whose Hessian rows and columns are computed. By default every movable MM atom is included. `--hessian-cutoff` keeps only the movable MM atoms within that distance (Å) of the ML region.
- **Movable MM by distance**: `--movable-cutoff` makes the MM atoms within that distance (Å) of the ML region movable and freezes the rest, in place of the B-factor layers.

## Features

### Link-atom redistribution

When the ML/MM boundary cuts a covalent bond, a link hydrogen caps the ML side of that bond in the MODEL-high calculation. The calculator adds one for every parm7 bond with exactly one end in the ML region (or for the pairs given in `link_mlmm`), so `model.pdb` contains no link hydrogens. `--link-atom-method` selects how the link hydrogen is placed:

| Method | Placement | Recommended |
| --- | --- | --- |
| **scaled** (g-factor, default) | `r_L = r_QM + g·(r_MM − r_QM)` with `g = (CR_QM + CR_H)/(CR_QM + CR_MM)` (covalent radii) | Yes (smooth PES, constant Jacobian) |
| **fixed** | `r_L = r_QM + d·û` with `d` = 1.09 Å (C parent) / 1.01 Å (N parent) and `û` the unit vector toward the MM atom | No (coordinate-dependent Jacobian) |

The scaled placement is the Morokuma–Dapprich g-factor method used in Gaussian ONIOM: the link hydrogen moves linearly with the QM–MM distance. The fixed placement keeps the link hydrogen at a constant distance along the bond axis.

Forces on the link hydrogen are passed to its ML parent and MM parent through the Jacobian `J` of the link position:

```
F_QM += (1−g) · F_link    (scaled)
F_MM += g · F_link
```

For the Hessian, the self term `Jᵀ H_link J` is added to the blocks of the two parent atoms. With the fixed placement, `J` depends on the coordinates, so the force-weighted second-derivative term `Σ (∂Jᵀ/∂x) f_L` is also added. The scaled placement is linear in the parent coordinates and needs no such term.

(microiteration)=
### Microiteration

When many MM atoms can move, optimizing all coordinates together evaluates the MLIP gradient at every step, even while only the MM environment relaxes. Microiteration, as in Gaussian 16, alternates two kinds of step:

```
relax the MM environment once, then repeat until converged:
    MACRO step  — one optimizer step on the ML atoms + the MM parents of link atoms (full ONIOM force)
    MICRO step  — L-BFGS relaxation of the other movable MM atoms (MM force only)
```

| | Macro step | Micro step |
|---|---|---|
| **Calculator** | Full ONIOM (`E_MM(real) + E_ML(model) − E_MM(model)`) | MM force field only (`E_MM(real)`) |
| **Coordinates optimized** | ML atoms + MM parents of link atoms | Movable MM atoms other than the MM parents of link atoms |
| **Optimizer** | `opt`: RFO from an exact Hessian, updated by TS-BFGS; `tsopt`: the selected Hessian TS optimizer | L-BFGS (no Hessian, started fresh at every micro step) |
| **Convergence** | `--thresh` (default `gau` in `opt`, `baker` in `tsopt`) | `microiter.micro_thresh` (default: same as `--thresh`) |

Microiteration is controlled by `--microiter/--no-microiter` (default: on). It runs in `opt --opt-mode hess` and in the Hessian TS modes of `tsopt` (`hess`, `rsirfo`, `rsprfo`, `trim`); `tsopt` uses `hess` by default. The micro-step keys are listed under [`microiter`](yaml-reference.md#microiter) in the YAML Reference.

```{note}
**Why the MM parents of link atoms move in the macro step:**
With the scaled (g-factor) link atom, `r_L = (1−g)·r_QM + g·r_MM` ties the link-atom position to **both** the QM and the MM parent. If the MM parent moved during the micro step (under MM forces with no ML contribution), the link atom would shift between cycles and the macro-step energy would oscillate. Moving the MM parents together with the ML atoms in the macro step, and holding them fixed in the micro step, removes this mismatch.
```

### MM Hessian

`--mm-backend` (YAML `calc.mm_backend`) selects the MM engine:

- **`hessian_ff`** (default): a CPU-only MM engine for Amber parm7 force fields, bundled with mlmm-toolkit. It evaluates the bond, angle, dihedral, improper, Lennard-Jones, electrostatic, and CMAP terms, and it can compute the MM Hessian analytically. The MM Hessian uses finite differences by default (`calc.mm_fd: true`); `calc.mm_fd: false` switches to the analytical `hessian_ff` Hessian. Its C++ kernels are built automatically on first use and need a C++20 compiler ([Installation](installation.md)). Running MM on the CPU leaves the GPU memory to the ML region.
- **`openmm`**: OpenMM on the CPU or CUDA, with a finite-difference Hessian. Use it for force fields that `hessian_ff` does not cover, or when OpenMM is already part of your workflow. YAML examples for `mm_backend` and `mm_device` and the VRAM trade-offs are in [Device Configuration & HPC Setup](device-hpc.md).

When only part of the system is active, the Hessian blocks of the active atoms can be expanded to the full Cartesian shape with the frozen rows and columns filled with zeros (`return_partial_hessian`).

### ML Hessian

`--hessian-calc-mode` (YAML `calc.hessian_calc_mode`) selects how the Hessian of the ML region is built:

- `FiniteDifference` (default): central differences of the forces. Works with every backend.
- `Analytical`: the autograd or native Hessian of the backend. UMA (with one worker), ORB, MACE, and AIMNet2 provide it when the installed version supports it, and the DFT backend provides it without `--embedcharge`. The per-backend details and the memory trade-offs are in [MLIP Backends › Hessian evaluation mode](backends.md#hessian-evaluation-mode).

The combined Hessian is assembled in float64 by default (`calc.H_double: true`).

### CMAP in the two MM layers

CMAP (Cross-Map backbone dihedral correction) is a 5-atom torsion correction term used by force fields such as ff19SB. The REAL and MODEL MM calculations must use the same CMAP policy in the subtractive expression.

| Region | E_MM(real) | E_MM(model) | ONIOM net effect |
|--------|-----------|------------|-----------------|
| `use_cmap: true` (default) | CMAP included when present | CMAP included when present | Complete model-internal CMAP cancels; boundary terms remain in the low-level coupling |
| `use_cmap: false` | CMAP excluded | CMAP excluded | Explicit modified-force-field calculation without CMAP |

For ff19SB, CMAP replaces the corresponding zeroed backbone cosine terms ([Tian et al., 2020](https://doi.org/10.1021/acs.jctc.9b00591)), so preserving it is the force-field-faithful default. `use_cmap: false` (CLI `--no-cmap`) removes CMAP from both MM layers.

**Example YAML configuration:**
```yaml
calc:
 use_cmap: false  # Explicitly remove CMAP from both MM layers
```

## Inputs

| Input | CLI | Description |
| --- | --- | --- |
| `input.pdb` | `-i` | Input structure; residue and atom names and the B-factor layers are read from it |
| `real.parm7` | `--parm7` | Amber topology of the full (REAL) system |
| `model.pdb` | `--model-pdb` | PDB that defines the ML region (used to identify the ML atoms) |

The atom order of `input.pdb` must match `real.parm7`. The calculator writes its own `real.rst7` from `real.parm7` and the coordinates of `input.pdb` with ParmEd, so no separate `real.rst7` or `real.pdb` is needed. On the command line, the ML region can also come from `--model-indices` or from the ML layer of the B-factors (`--detect-layer`).

## Units

| Quantity | Internal unit | PySisyphus interface |
| --- | --- | --- |
| Energy | eV | Hartree |
| Forces | eV/Å | Hartree/Bohr |
| Hessian | eV/Å² | Hartree/Bohr² |

## Python API

The calculator can also be used from Python without the CLI. The `mlmm` package exports `MLMMCore` (the engine), `MLMMASECalculator` (ASE interface), and `mlmm` (pysisyphus Calculator). The constructor arguments correspond to the CLI options: `-i` → `input_pdb`, `--parm7` → `real_parm7`, `--model-pdb` → `model_pdb`, `-q`/`-m` → `model_charge`/`model_mult`, `-b` → `backend`, and `--mm-backend` → `mm_backend`. Most other arguments share their names with the keys of the YAML [`calc` section](yaml-reference.md#calc-section).

### Quick start

```python
# cd examples/methyltransferase
from ase.io import read
from mlmm import MLMMCore

# Base engine: returns energy (eV), forces (eV/Å), Hessian (eV/Å²)
core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,   # net charge of the ML region in this example
)

coords = read("r_layered.pdb").get_positions()   # shape (N, 3), Å
result = core.compute(coords, return_forces=True, return_hessian=False)
print(result["energy"], result["forces"].shape)
```

### API levels

mlmm-toolkit provides three API levels:

| Level | Class | Input units | Output units | Use case |
|-------|-------|-------------|--------------|----------|
| Base engine | `MLMMCore` | Å | eV, eV/Å, eV/Å² | Direct Python scripting |
| ASE | `MLMMASECalculator` | Å (via `Atoms`) | eV, eV/Å | ASE-based workflows (DMF, MD) |
| pysisyphus | `mlmm` (Calculator) | Bohr (via `Geometry`) | Hartree, Hartree/Bohr | pysisyphus optimization, IRC, freq |

### MLMMCore

The core ML/MM engine. It sets up the topology, the force field, and the MLIP backend once; each `compute()` call then only updates the coordinates.

```python
from mlmm import MLMMCore

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
    model_mult=1,
    backend="uma",               # uma | orb | mace | aimnet2 | dft (dft also needs dft_settings)
    return_partial_hessian=True, # partial Hessian for the Hessian-target atoms
)
```

#### Key parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_pdb` | `str` | *required* | Input PDB (full system with B-factor layers) |
| `real_parm7` | `str` | *required* | Amber prmtop for the full system |
| `model_pdb` | `str` | *required* | PDB defining the ML region |
| `model_charge` | `int` | `0` | Net charge of the ML region. The constructor logs it at INFO level, so you can check it for charged systems. |
| `model_mult` | `int` | `1` | Spin multiplicity of the ML region |
| `backend` | `str` | `"uma"` | MLIP backend |
| `mm_backend` | `str` | `"hessian_ff"` | MM engine (`hessian_ff` or `openmm`) |
| `return_partial_hessian` | `bool` | `True` | If `True`, `compute()` returns a 4D `(n_active, 3, n_active, 3)` sub-Hessian plus a `within_partial_hessian` metadata dict (active-atom indices, DOF maps). If `False`, it returns the expanded 4D `(N, 3, N, 3)` full-system Hessian. |
| `link_mlmm` | `List[Tuple[str, str]]` | `None` | Manual link-atom pairs as `[("RESN RESID ATOMNAME", "RESN RESID ATOMNAME"), ...]` (first = ML side, second = MM side, e.g. `[("SAM 359 CA", "SAM 359 N")]`). `None` adds a link for every bond in the parm7 topology that crosses the ML/MM boundary (not by distance). |

#### compute()

```python
result = core.compute(
    coord_ang,                   # numpy (N, 3), Angstrom
    return_forces=True,
    return_hessian=False,
)
# result["energy"]   : float (eV)
# result["forces"]   : numpy (N, 3) (eV/Å)
# result["hessian"]  : torch 4D (eV/Å²), only if return_hessian=True.
#   Shape (n_active, 3, n_active, 3) when return_partial_hessian=True (default),
#   else expanded to (N, 3, N, 3). With the partial Hessian, the accompanying
#   key `within_partial_hessian` (dict with active_atoms / active_dofs /
#   full_to_active mappings) is also returned.
```

### MLMMASECalculator

An ASE `Calculator` that wraps `MLMMCore` and returns energy and forces. It works with ASE optimizers, MD, and DMF.

```python
from mlmm import MLMMCore, MLMMASECalculator
from ase.io import read

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
)
calc = MLMMASECalculator(core)

atoms = read("r_layered.pdb")
atoms.calc = calc
print(atoms.get_potential_energy())   # eV
print(atoms.get_forces().shape)       # (N, 3), eV/Å
```

### pysisyphus Calculator (`mlmm`)

For pysisyphus optimization, IRC, and frequency analysis. It takes the same arguments as `MLMMCore`.

```python
from mlmm import mlmm as MLMMCalc
from pysisyphus.io import geom_from_pdb

calc = MLMMCalc(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=-1,
)
geom = geom_from_pdb("r_layered.pdb")
geom.set_calculator(calc)
energy = geom.energy            # Hartree
forces = geom.forces            # Hartree/Bohr (flat)
```

## Notes

- `--hessian-calc-mode Analytical` stops with an error when the selected backend, or its installed version, has no analytical Hessian. The calculator does not switch to finite differences on its own.
- `use_cmap: false` (`--no-cmap`) is not an ff19SB-compatible setting; it removes CMAP from both MM layers.
- Microiteration is turned off when `--embedcharge` is on, because the MM-only micro steps would leave out the embedding forces; the optimizer then moves the ML and MM atoms together.
- The fixed link-atom placement supports only C and N parents on the ML side.
- With the default `return_partial_hessian=True`, `compute()` returns the Hessian of the Hessian-target atoms only, as a 4D array `(n_active, 3, n_active, 3)`. Use `within_partial_hessian` to map it back to the full system.
- `backend="dft"` also needs `dft_settings`. On the command line, `-b dft` builds them from the YAML `calc.dft` block and the DFT options ([Refine an MLIP TS with DFT](dft-backend.md)).

## See Also

- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [opt](opt.md) — Single-structure geometry optimization using the ML/MM calculator
- [tsopt](tsopt.md) — Transition state optimization
- [freq](freq.md) — Vibrational frequency analysis
- [YAML Reference](yaml-reference.md) — `calc` and `microiter` configuration keys
- [MLIP Backends](backends.md) — Backend selection, install, precision, and the add-a-backend recipe
- [Device Configuration & HPC Setup](device-hpc.md) — ML/MM device settings and HPC submission
