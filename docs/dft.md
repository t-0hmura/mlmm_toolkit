# `dft`

## Overview

> **Summary:** Run a single-point DFT calculation on the ML region using PySCF/GPU4PySCF, then recombine with MM energies to obtain the ML(dft)/MM total energy.

`mlmm dft` extracts the ML region from the full enzyme PDB, appends link hydrogens, and runs a single-point PySCF (or GPU4PySCF) calculation. After the DFT evaluation, the script recomputes the **ML(dft)/MM total energy** by combining the PySCF high-level energy with MM evaluations of the full system (REAL-low) and the ML subset (MODEL-low):

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

The GPU4PySCF backend is activated automatically when available; otherwise PySCF CPU is used. The default functional/basis is `wb97m-v/6-31g**`.

## Usage
```bash
mlmm dft -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m SPIN] [--freeze-atoms "1,3,5"] [--func-basis "FUNC/BASIS"] \
    [--max-cycle N] [--conv-tol Eh] [--grid-level L] [--out-dir DIR] [--args-yaml FILE]
```

### Examples
```bash
# Basic single-point DFT
mlmm dft -i enzyme.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1

# Custom functional/basis, frozen atoms, tighter SCF
mlmm dft -i enzyme.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q -1 -m 2 \
    --func-basis "wb97m-v/def2-tzvpd" --freeze-atoms "1,3,5" --max-cycle 150 --conv-tol 1e-9
```

## Workflow
1. **Input handling** -- The full enzyme PDB (`-i`), Amber topology (`--real-parm7`), and ML-region definition (`--model-pdb` or `--model-indices` or B-factor detection via `--detect-layer`) are loaded. Link hydrogens are appended automatically (C/N parents within 1.7 A) unless explicit `link_mlmm` pairs are provided via YAML.
2. **Configuration merge** -- Defaults -> CLI -> YAML (`geom`, `calc`/`mlmm`, `dft` blocks). YAML overrides follow the precedence **CLI > YAML > defaults**.
3. **SCF build** -- `--func-basis` is parsed into functional and basis. Density fitting and nonlocal corrections follow PySCF/GPU4PySCF defaults.
4. **ML(dft)/MM recombination** -- After the DFT converges, MM evaluations of the full system (REAL-low) and the ML subset (MODEL-low) are computed. The combined energy is reported in Hartree and kcal/mol.
5. **Population analysis & outputs** -- Mulliken, meta-Lowdin, and IAO charges and spin densities (UKS only) are written alongside the combined energy block.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme PDB file (must be `.pdb`). | Required |
| `--real-parm7 PATH` | Amber parm7 topology for the full system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs must match the enzyme PDB). Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20). | `True` |
| `-q, --charge INT` | Charge of the ML region. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (e.g. `"1,3,5"`). Merged with YAML `geom.freeze_atoms`. | _None_ |
| `--func-basis TEXT` | Functional/basis pair as `"FUNC/BASIS"`. | `wb97m-v/6-31g**` |
| `--max-cycle INT` | Maximum SCF iterations. | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance (Hartree). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level. | `3` |
| `--out-dir DIR` | Output directory. | `./result_dft/` |
| `--args-yaml FILE` | YAML overrides under `geom`, `calc`/`mlmm`, `dft`. | _None_ |

## Outputs
```
out_dir/  (default: ./result_dft/)
├── ml_region_with_linkH.xyz    # ML-region coordinates (with link-H) used for DFT
├── result.yaml                 # DFT + ML(dft)/MM energy summary, charges, spin densities
└── (stdout)                    # Pretty-printed configuration blocks and energies
```

## YAML configuration (`--args-yaml`)

Accepts a mapping root; the `dft` section (and optional `geom`, `calc`/`mlmm`) is applied when present. YAML values override CLI values.

`dft` keys (defaults in parentheses):
- `func_basis` (`"wb97m-v/6-31g**"`): Combined `FUNC/BASIS` string.
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`4`): PySCF verbosity (0-9).
- `out_dir` (`"./result_dft/"`): Output directory root.

```yaml
geom:
  coord_type: cart
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
dft:
  func_basis: wb97m-v/6-31g**
  conv_tol: 1.0e-09
  max_cycle: 100
  grid_level: 3
  verbose: 4
  out_dir: ./result_dft/
```

## Notes
- Link hydrogens are detected automatically (C/N parents within 1.7 A) unless explicit `link_mlmm` pairs are provided via YAML. Unsupported parent elements raise an error.
- DFT options (functional/basis, SCF controls) remain YAML-overridable under the `dft` key.
- Exit codes: `0` when the SCF converges, `3` when it does not, `130` for keyboard interrupt, `1` on other errors.

---

## See Also

- [freq](freq.md) -- Vibrational frequency analysis (often precedes DFT refinement)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow with `--dft True`
