# mlmm pysis — pysisyphus YAML Workflows

> **Summary:** Run pysisyphus YAML workflow files with the mlmm calculator pre-registered. Compatible with v0.1.x YAML-based workflows.

## Overview

The `mlmm pysis` subcommand executes pysisyphus YAML workflow files with the `mlmm` calculator type automatically registered. This provides backward compatibility with v0.1.x YAML-based workflows where `mlmm opt.yaml` was the primary interface.

```bash
mlmm pysis opt.yaml
mlmm pysis tsopt.yaml
mlmm pysis irc.yaml
```

## YAML Format

pysisyphus YAML files have three main sections: `geom`, `calc`, and a workflow section (`opt`, `tsopt`, `irc`, etc.).

### Geometry section

```yaml
geom:
  type: cart          # Coordinate system (cart, redund, dlc)
  fn: complex_layered.pdb
```

For frozen atoms, use `freeze_atoms` (1-based indices):

```yaml
geom:
  type: cart
  fn: complex_layered.pdb
  freeze_atoms: [101, 102, 103, 200, 201]  # 1-based
```

### Calculator section

```yaml
calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
  model_mult: 1
  backend: uma
  embedcharge: false
```

All `MLMMCore` parameters are accepted as YAML keys. v0.1.x parameter names (`real_pdb`, etc.) are also accepted with a deprecation warning.

### Workflow sections

#### Optimization

```yaml
opt:
  type: lbfgs
  max_cycles: 300
  thresh: gau_loose    # gau, gau_tight, gau_vtight, gau_loose, baker
  dump: true
```

#### TS optimization

```yaml
tsopt:
  type: dimer
  max_cycles: 100
  thresh: gau_loose
  dump: true
```

#### IRC

```yaml
irc:
  type: eulerpc
  max_cycles: 75
  step_length: 0.15
  dump: true
```

## Full YAML Examples

### Geometry optimization

```yaml
geom:
  type: cart
  fn: complex_layered.pdb
  freeze_atoms: [101, 102, 200, 201]

calc:
  type: mlmm
  input_pdb: complex_layered.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0
  backend: uma

opt:
  type: lbfgs
  max_cycles: 300
  thresh: gau_loose
  dump: true
```

Run: `mlmm pysis opt.yaml`

### TS optimization

```yaml
geom:
  type: cart
  fn: ts_candidate.pdb
  freeze_atoms: [101, 102, 200, 201]

calc:
  type: mlmm
  input_pdb: ts_candidate.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0

tsopt:
  type: dimer
  max_cycles: 100
  thresh: gau_loose
  dump: true
```

Run: `mlmm pysis tsopt.yaml`

### IRC from TS

```yaml
geom:
  type: cart
  fn: ts_optimized.pdb

calc:
  type: mlmm
  input_pdb: ts_optimized.pdb
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  model_charge: 0

irc:
  type: eulerpc
  max_cycles: 75
  step_length: 0.15
  dump: true
```

Run: `mlmm pysis irc.yaml`

## Migration from v0.1.x

| v0.1.x | v0.2.x |
|--------|--------|
| `mlmm opt.yaml` | `mlmm pysis opt.yaml` |
| `real_pdb: complex.pdb` | `input_pdb: complex.pdb` (or keep `real_pdb` with warning) |
| `real_rst7: real.rst7` | No longer needed (ignored with warning) |

Existing v0.1.x YAML files work without modification — deprecated parameter names are accepted with warnings.

## See Also

- [Python API](python_api.md) — Using mlmm-toolkit as a Python library
- [ML/MM Calculator](mlmm_calc.md) — Calculator architecture
- [YAML Reference](yaml_reference.md) — Full YAML configuration reference
