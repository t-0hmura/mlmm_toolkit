# `define-layer`

## Overview

> **Summary:** Define a 3-layer ML/MM system based on distance from the ML region and encode the layer assignments as B-factors in the output PDB.

`mlmm define-layer` partitions an enzyme system into three layers around the ML region and writes the assignments as PDB B-factors. The ML region can be specified via a model PDB, explicit atom indices, or a combination of both.

### 3-Layer System
| Layer | Name | B-factor | Description |
| --- | --- | --- | --- |
| 1 | ML | 0.0 | Atoms in the ML region |
| 2 | Movable-MM | 10.0 | MM atoms/residues within `--radius-freeze` of ML |
| 3 | Frozen | 20.0 | MM atoms/residues beyond `--radius-freeze` |

### Layer assignment strategy
- **Residues without ML atoms:** the entire residue is assigned to a single layer based on the minimum distance from any ML atom to any atom in the residue.
- **Residues with ML atoms:** non-ML atoms in the same residue are classified individually by distance.

## Usage
```bash
mlmm define-layer -i INPUT.pdb --model-pdb MODEL.pdb [--model-indices "0,1,2,..."] \
 [--radius-partial-hessian FLOAT] [--radius-freeze FLOAT] \
 [-o OUTPUT.pdb] [--one-based|--zero-based]
```

### Examples
```bash
# Using model PDB to define ML region
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb

# Using explicit atom indices (0-based)
mlmm define-layer -i system.pdb --model-indices "0,1,2,3,4" --zero-based -o labeled.pdb

# Custom movable/frozen cutoff
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb \
 --radius-freeze 10.0 -o labeled.pdb
```

## Description
1. **ML region identification** -- The ML region is defined by `--model-pdb` (atom matching against the input PDB) or `--model-indices` (explicit atom indices). If `--model-indices` is provided, it takes precedence over `--model-pdb`.
2. **Distance computation** -- For each non-ML atom (or residue), the minimum distance from any ML atom is computed.
3. **Layer assignment** -- Non-ML atoms/residues are assigned to Movable-MM or Frozen by `--radius-freeze`.
4. **Output** -- The output PDB has B-factors set to layer values (0, 10, 20). A summary of layer assignments is printed to the console.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file containing the full system. | Required |
| `--model-pdb PATH` | PDB file defining atoms in the ML region. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (e.g. `"0,1,2,3"` or `"1-10,15,20-25"`). Takes precedence over `--model-pdb`. | _None_ |
| `--radius-partial-hessian FLOAT` | Deprecated in 3-layer mode (kept ). | `0.0` |
| `--radius-freeze FLOAT` | Distance cutoff (A) from ML region for Movable-MM. Atoms beyond this are Frozen. | `8.0` |
| `-o, --output PATH` | Output PDB file with B-factors set to layer values. | `<input>_layered.pdb` |
| `--one-based / --zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |

## Outputs
```
<output>.pdb # PDB with B-factors set to 0 / 10 / 20
(stdout) # Summary table of layer assignments and atom counts
```

## Notes
- Distance is calculated as the minimum distance from any ML atom to any atom in the residue.
- `--radius-partial-hessian` is accepted but ignored in 3-layer mode.
- The output PDB preserves all original atom records; only the B-factor column is modified.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [mm_parm](mm_parm.md) -- Build AMBER topology (parm7/rst7) before layer definition
- [opt](opt.md) -- Single-structure optimization using the layered system
- [all](all.md) -- End-to-end workflow that includes automatic layer definition
