# Curated `mlmm all` Starter Snapshot

This page is a **curated, non-exhaustive** starter snapshot for `mlmm all`. It shows a common subset of keys whose values are pinned to (and equal) their runtime owners; it is **not** the full configuration schema. For every configurable section and option, see the [YAML Reference](../yaml-reference.md).

- Source template: `.github/scripts/generate_reference.py::_ALL_TEMPLATE`
- Template digest: `10aefa8eaab4`

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
  backend: uma              # ML backend: uma, orb, mace, aimnet2
  orb_model: orb_v3_conservative_omol  # ORB model name (when backend=orb)
  orb_precision: float64    # ORB precision default (when backend=orb; "float32-high" = TF32 matmul, also via --precision fp32; legacy "float32" alias accepted)
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

Each scalar is pinned to (and equals) the runtime owner shown.

| Key | Type | Default | Runtime owner |
|---|---|---|---|
| `calc.backend` | `str` | `'uma'` | `MLMM_CALC_KW["backend"]` |
| `calc.orb_model` | `str` | `'orb_v3_conservative_omol'` | `MLMM_CALC_KW["orb_model"]` |
| `calc.orb_precision` | `str` | `'float64'` | `MLMM_CALC_KW["orb_precision"]` |
| `calc.mace_model` | `str` | `'MACE-OMOL-0'` | `MLMM_CALC_KW["mace_model"]` |
| `calc.mace_dtype` | `str` | `'float64'` | `MLMM_CALC_KW["mace_dtype"]` |
| `calc.aimnet2_model` | `str` | `'aimnet2'` | `MLMM_CALC_KW["aimnet2_model"]` |
| `freq.max_write` | `int` | `10` | `FREQ_KW["max_write"]` |
| `freq.amplitude_ang` | `float` | `0.8` | `FREQ_KW["amplitude_ang"]` |
| `freq.n_frames` | `int` | `20` | `FREQ_KW["n_frames"]` |
| `freq.sort` | `str` | `'value'` | `FREQ_KW["sort"]` |
| `thermo.temperature` | `float` | `298.15` | `THERMO_KW["temperature"]` |
| `thermo.pressure_atm` | `float` | `1.0` | `THERMO_KW["pressure_atm"]` |
| `thermo.symmetry_number` | `NoneType` | `None` | `THERMO_KW["symmetry_number"]` |
| `dft.func_basis` | `str` | `'wb97m-v/def2-svp'` | `DFT_KW["func_basis"]` |
| `dft.max_cycle` | `int` | `100` | `DFT_KW["max_cycle"]` |
| `dft.conv_tol` | `float` | `1e-09` | `DFT_KW["conv_tol"]` |
| `dft.grid_level` | `int` | `3` | `DFT_KW["grid_level"]` |

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
