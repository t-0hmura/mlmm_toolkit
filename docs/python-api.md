# Python API

Use mlmm-toolkit as a Python library — `MLMMCore` (base engine), `MLMMASECalculator` (ASE interface), and `mlmm` (pysisyphus Calculator).

## Quick Start

```python
# cd examples/methyltransferase
from mlmm import MLMMCore, MLMMASECalculator, mlmm

# Base engine — returns energy (eV), forces (eV/Å), Hessian (eV/Å²)
core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=0,
)

import numpy as np
coords = np.loadtxt(...)  # shape (N, 3), Angstrom
result = core.compute(coords, return_forces=True, return_hessian=False)
print(result["energy"], result["forces"].shape)
```

## API Levels

mlmm-toolkit provides three API levels depending on the context:

| Level | Class | Input units | Output units | Use case |
|-------|-------|-------------|--------------|----------|
| Base engine | `MLMMCore` | Å | eV, eV/Å, eV/Å² | Direct Python scripting |
| ASE | `MLMMASECalculator` | Å (via `Atoms`) | eV, eV/Å | ASE-based workflows (DMF, MD) |
| pysisyphus | `mlmm` (Calculator) | Bohr (via `Geometry`) | Hartree, Hartree/Bohr | pysisyphus optimization, IRC, freq |

## MLMMCore

The core ML/MM engine. Initializes topology, force field, and MLIP backend once; subsequent `compute()` calls only update coordinates.

```python
from mlmm import MLMMCore

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=0,
    model_mult=1,
    backend="uma",               # uma | orb | mace | aimnet2 | dft (dft also needs dft_settings)
    return_partial_hessian=True, # partial Hessian for the active Hessian atoms
)
```

### Key parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_pdb` | `str` | *required* | Input PDB (full system with B-factor layers) |
| `real_parm7` | `str` | *required* | Amber prmtop for the full system |
| `model_pdb` | `str` | *required* | PDB defining the ML region |
| `model_charge` | `int` | `0` | Net charge of the ML region. The constructor logs the resolved ML-region net charge at INFO level so you can verify it matches your expectation for charged systems. |
| `model_mult` | `int` | `1` | Spin multiplicity of the ML region |
| `backend` | `str` | `"uma"` | MLIP backend |
| `mm_backend` | `str` | `"hessian_ff"` | MM engine (`hessian_ff` or `openmm`) |
| `return_partial_hessian` | `bool` | `True` | If `True`, `compute()` returns a 4D `(n_active, 3, n_active, 3)` sub-Hessian plus a `within_partial_hessian` metadata dict (active-atom indices, DOF maps). If `False`, returns the expanded 4D `(N, 3, N, 3)` full-system Hessian. |
| `link_mlmm` | `List[Tuple[str, str]]` | `None` | Manual link-atom pairs as `[("RESN RESID ATOMNAME", "RESN RESID ATOMNAME"), ...]` (first = ML-side, second = MM-side, e.g. `[("SAM 359 CA", "SAM 359 N")]`). `None` = derive every crossing bond from the supplied parm7 topology (not from distance). |

### compute()

```python
result = core.compute(
    coord_ang,                   # numpy (N, 3), Angstrom
    return_forces=True,
    return_hessian=False,
)
# result["energy"]   : float (eV)
# result["forces"]   : numpy (N, 3) (eV/Å)
# result["hessian"]  : torch 4D (eV/Å²) — only if return_hessian=True.
#   Shape (n_active, 3, n_active, 3) when return_partial_hessian=True (default),
#   else expanded to (N, 3, N, 3). Companion key `within_partial_hessian`
#   (dict with active_atoms / active_dofs / full_to_active mappings) accompanies
#   the partial-Hessian result.
```

## MLMMASECalculator

ASE `Calculator` wrapping `MLMMCore`. Compatible with ASE optimizers, MD, and DMF.

```python
from mlmm import MLMMCore, MLMMASECalculator
from ase.io import read

core = MLMMCore(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
)
calc = MLMMASECalculator(core)

atoms = read("r_layered.pdb")
atoms.calc = calc
print(atoms.get_potential_energy())   # eV
print(atoms.get_forces().shape)       # (N, 3), eV/Å
```

## pysisyphus Calculator (`mlmm`)

For use with pysisyphus optimization, IRC, and frequency analysis.

```python
from mlmm import mlmm as MLMMCalc
from pysisyphus.helpers import geom_loader

calc = MLMMCalc(
    input_pdb="r_layered.pdb",
    real_parm7="complex.parm7",
    model_pdb="pocket.pdb",
    model_charge=0,
)
geom = geom_loader("r_layered.pdb")
geom.set_calculator(calc)
energy = geom.energy            # Hartree
forces = geom.forces            # Hartree/Bohr (flat)
```

## See Also

- [ML/MM Calculator](mlmm-calc.md) — Architecture and internal details
- [YAML Reference](yaml-reference.md) — Configuration keys for `--config` YAML
