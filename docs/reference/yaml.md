# YAML Schema

- Source template: `scripts/generate_reference.py::_ALL_TEMPLATE`
- Template digest: `63bc0343cb6d`

## Top-level Keys

| Key |
|---|
| `calc` |
| `extract` |
| `path_search` |
| `scan` |
| `tsopt` |
| `freq` |
| `dft` |

## Starter Template

```yaml
# Starter config for `mlmm all`

calc:
  backend: uma              # ML backend: uma, orb, mace, aimnet2
  orb_model: null           # ORB model name (when backend=orb)
  orb_precision: null       # ORB precision setting (when backend=orb)
  mace_model: null          # MACE model path or name (when backend=mace)
  mace_dtype: null          # MACE dtype, e.g. float64 (when backend=mace)
  aimnet2_model: null       # AIMNet2 model name (when backend=aimnet2)
  embedcharge: false        # Enable xTB point-charge embedding correction
  embedcharge_step: 1.0e-3   # Numerical Hessian step for embedding correction (Å)
  xtb_cmd: xtb              # Path or command for the xTB executable
  xtb_acc: 0.2              # xTB SCF accuracy parameter
  xtb_workdir: tmp          # Working directory for xTB scratch files
  xtb_keep_files: false     # Keep xTB intermediate files after completion
  xtb_ncores: 4             # Number of CPU cores for xTB

extract:
  radius: 2.6
  radius_het2het: 0.0

path_search:
  max_nodes: 10
  max_cycles: 300

scan:
  max_step_size: 0.2
  bias_k: 300.0
  relax_max_cycles: 10000

tsopt:
  max_cycles: 10000

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value
  temperature: 298.15
  pressure_atm: 1.0

dft:
  func_basis: wb97m-v/def2-tzvpd
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## Scalar Defaults

| Key | Type | Default |
|---|---|---|
| `calc.backend` | `str` | `'uma'` |
| `calc.orb_model` | `NoneType` | `None` |
| `calc.orb_precision` | `NoneType` | `None` |
| `calc.mace_model` | `NoneType` | `None` |
| `calc.mace_dtype` | `NoneType` | `None` |
| `calc.aimnet2_model` | `NoneType` | `None` |
| `calc.embedcharge` | `bool` | `False` |
| `calc.embedcharge_step` | `float` | `0.001` |
| `calc.xtb_cmd` | `str` | `'xtb'` |
| `calc.xtb_acc` | `float` | `0.2` |
| `calc.xtb_workdir` | `str` | `'tmp'` |
| `calc.xtb_keep_files` | `bool` | `False` |
| `calc.xtb_ncores` | `int` | `4` |
| `extract.radius` | `float` | `2.6` |
| `extract.radius_het2het` | `float` | `0.0` |
| `path_search.max_nodes` | `int` | `10` |
| `path_search.max_cycles` | `int` | `300` |
| `scan.max_step_size` | `float` | `0.2` |
| `scan.bias_k` | `float` | `300.0` |
| `scan.relax_max_cycles` | `int` | `10000` |
| `tsopt.max_cycles` | `int` | `10000` |
| `freq.max_write` | `int` | `10` |
| `freq.amplitude_ang` | `float` | `0.8` |
| `freq.n_frames` | `int` | `20` |
| `freq.sort` | `str` | `'value'` |
| `freq.temperature` | `float` | `298.15` |
| `freq.pressure_atm` | `float` | `1.0` |
| `dft.func_basis` | `str` | `'wb97m-v/def2-tzvpd'` |
| `dft.max_cycle` | `int` | `100` |
| `dft.conv_tol` | `float` | `1e-09` |
| `dft.grid_level` | `int` | `3` |

## Scan Spec Shapes

Accepted by `scan`, `scan2d`, and `scan3d` with `-s/--scan-lists`.

```yaml
# scan (1D staged)
one_based: false
stages:
  - - [1, 2, 1.65]
  - - [2, 3, 2.30]

# scan2d / scan3d
one_based: false
pairs:
  - [1, 2, 1.40, 2.20]
  - [2, 3, 1.20, 2.00]
  - [3, 4, 1.00, 1.80]  # required only for scan3d
```
