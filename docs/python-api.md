# Python API

> **Summary:** Use mlmm-toolkit as a Python library — `MLMMCore` (base engine), `MLMMASECalculator` (ASE interface), `mlmm` (pysisyphus Calculator), and the `mlmm_ase()` compatibility wrapper.

## Quick Start

```python
from mlmm import MLMMCore, MLMMASECalculator, mlmm

# Base engine — returns energy (eV), forces (eV/Å), Hessian (eV/Å²)
core = MLMMCore(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
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
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
    model_charge=0,
    model_mult=1,
    backend="uma",               # uma | orb | mace | aimnet2
    embedcharge=False,           # xTB point-charge correction
    return_partial_hessian=True, # partial Hessian (ML region only)
)
```

### Key parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_pdb` | `str` | *required* | Input PDB (full system with B-factor layers) |
| `real_parm7` | `str` | *required* | Amber prmtop for the full system |
| `model_pdb` | `str` | *required* | PDB defining the ML region |
| `model_charge` | `int` | `0` | Net charge of the ML region |
| `model_mult` | `int` | `1` | Spin multiplicity of the ML region |
| `backend` | `str` | `"uma"` | MLIP backend |
| `embedcharge` | `bool` | `False` | Enable xTB point-charge embedding |
| `mm_backend` | `str` | `"hessian_ff"` | MM engine (`hessian_ff` or `openmm`) |
| `return_partial_hessian` | `bool` | `True` | Return partial Hessian (ML + boundary) |
| `link_mlmm` | `list` | `None` | Manual link atom specification |

### compute()

```python
result = core.compute(
    coord_ang,                   # numpy (N, 3), Angstrom
    return_forces=True,
    return_hessian=False,
)
# result["energy"]   : float (eV)
# result["forces"]   : numpy (N, 3) (eV/Å)
# result["hessian"]  : torch (3N, 3N) (eV/Å²)  — only if return_hessian=True
```

## MLMMASECalculator

ASE `Calculator` wrapping `MLMMCore`. Compatible with ASE optimizers, MD, and DMF.

```python
from mlmm import MLMMCore, MLMMASECalculator
from ase.io import read

core = MLMMCore(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
)
calc = MLMMASECalculator(core)

atoms = read("complex_layered.pdb")
atoms.calc = calc
print(atoms.get_potential_energy())   # eV
print(atoms.get_forces().shape)       # (N, 3), eV/Å
```

## pysisyphus Calculator (`mlmm`)

For use with pysisyphus optimization, IRC, frequency analysis.

```python
from mlmm import mlmm as MLMMCalc
from pysisyphus.helpers import geom_loader

calc = MLMMCalc(
    input_pdb="complex_layered.pdb",
    real_parm7="real.parm7",
    model_pdb="ml_region.pdb",
    model_charge=0,
)
geom = geom_loader("complex_layered.pdb")
geom.set_calculator(calc)
energy = geom.energy            # Hartree
forces = geom.forces            # Hartree/Bohr (flat)
```

The `mlmm` calculator is also registered in pysisyphus's `CALC_DICT`, so it can be used in YAML workflows:

```yaml
calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
```

See [mlmm pysis](pysis.md) for running YAML workflows.

## v0.1.x Compatibility

### Parameter aliases

The following v0.1.x parameter names are accepted with a `DeprecationWarning`:

| v0.1.x name | v0.2.x name | Notes |
|-------------|-------------|-------|
| `real_pdb` | `input_pdb` | Alias — maps to `input_pdb` |
| `real_rst7` | *(removed)* | Ignored with warning (auto-generated internally) |
| `vib_run` | *(removed)* | Ignored with warning |
| `vib_dir` | *(removed)* | Ignored with warning |

```python
# v0.1.x style — still works, emits DeprecationWarning
from mlmm import MLMMCore
core = MLMMCore(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
```

### mlmm_ase() factory

The v0.1.x `mlmm_ase(real_pdb=..., ...)` convenience function is preserved:

```python
from mlmm import mlmm_ase

# v0.1.x style — emits DeprecationWarning
calc = mlmm_ase(real_pdb="complex.pdb", real_parm7="real.parm7", model_pdb="ml.pdb")
# Equivalent to: MLMMASECalculator(MLMMCore(input_pdb="complex.pdb", ...))
```

### import path

v0.2.x uses the same import path as v0.1.x:

```python
import mlmm
from mlmm import MLMMCore
```

## See Also

- [ML/MM Calculator](mlmm-calc.md) — Architecture and internal details
- [mlmm pysis](pysis.md) — Running pysisyphus YAML workflows
- [YAML Reference](yaml-reference.md) — Configuration keys for YAML mode
