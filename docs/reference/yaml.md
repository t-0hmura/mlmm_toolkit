# Curated `mlmm all` Starter Snapshot

This page is a **curated, non-exhaustive** starter snapshot for `mlmm all`. It shows a common subset of keys with their current default values; it is **not** the full configuration schema. For every configurable section and option, see the [YAML Reference](../yaml-reference.md).

## Included Sections

| Section |
|---|
| `calc` |
| `freq` |
| `thermo` |
| `dft` |

## Starter Template

```yaml
# Starter config for `mlmm all`

calc:
  backend: uma              # High-level backend: uma, orb, mace, aimnet2, or dft
  orb_model: orbmol-v2      # ORB model name (when backend=orb)
  orb_precision: float64    # ORB precision default (when backend=orb; "float32-high" = TF32 matmul, also via --precision tf32; "float32-highest" = IEEE FP32, also via --precision fp32 or "float32")
  mace_model: MACE-OMOL-0   # MACE model path or name (when backend=mace)
  mace_dtype: float64       # MACE dtype, e.g. float32 / float64 (when backend=mace)
  aimnet2_model: aimnet2    # AIMNet2 model name (when backend=aimnet2)

freq:
  max_write: 10
  amplitude_ang: 0.8
  n_frames: 20
  sort: value

thermo:
  temperature: 298.15
  pressure_atm: 1.0
  symmetry_number: null

dft:
  func_basis: wb97m-v/def2-svp
  max_cycle: 100
  conv_tol: 1.0e-9
  grid_level: 3
```

## Scalar Defaults

| Key | Type | Default |
|---|---|---|
| `calc.backend` | `str` | `'uma'` |
| `calc.orb_model` | `str` | `'orbmol-v2'` |
| `calc.orb_precision` | `str` | `'float64'` |
| `calc.mace_model` | `str` | `'MACE-OMOL-0'` |
| `calc.mace_dtype` | `str` | `'float64'` |
| `calc.aimnet2_model` | `str` | `'aimnet2'` |
| `freq.max_write` | `int` | `10` |
| `freq.amplitude_ang` | `float` | `0.8` |
| `freq.n_frames` | `int` | `20` |
| `freq.sort` | `str` | `'value'` |
| `thermo.temperature` | `float` | `298.15` |
| `thermo.pressure_atm` | `float` | `1.0` |
| `thermo.symmetry_number` | `NoneType` | `None` |
| `dft.func_basis` | `str` | `'wb97m-v/def2-svp'` |
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

# scan2d
one_based: false
pairs:
  - [1, 2, 1.40, 2.20]
  - [2, 3, 1.20, 2.00]

# scan3d
one_based: false
pairs:
  - [1, 2, 1.40, 2.20]
  - [2, 3, 1.20, 2.00]
  - [3, 4, 1.00, 1.80]
```
