# `dft`

## Overview

> **Summary:** Run a single-point DFT calculation on the ML region using PySCF/GPU4PySCF, then recombine with MM energies to obtain the ML(dft)/MM total energy. Results include energy and population analysis (Mulliken, meta-Lowdin, IAO charges).

### At a glance
- **Use when:** You want a higher-level single-point energy for your ML/MM system using DFT on the ML region with ONIOM-style recombination.
- **Method:** PySCF (CPU) or GPU4PySCF (GPU) single-point DFT on the ML region + link hydrogens, combined with MM evaluations for the ONIOM total energy.
- **Outputs:** `ml_region_with_linkH.xyz`, `result.yaml` with ML(dft)/MM combined energy.
- **Defaults:** `--func-basis wb97m-v/def2-tzvpd`, `--max-cycle 100`, `--conv-tol 1e-9`.
- **Next step:** Compare DFT//UMA energies across R/TS/P states, or use within [all](all.md) `--dft` for automated diagrams.
- **Prerequisites:** DFT dependencies (PySCF, GPU4PySCF) are **not** included in the default install. Install them with `pip install mlmm_toolkit[dft]`.
- **System size limit:** DFT single-point calculations are practical only for ML regions up to **~500 atoms**. Larger ML regions will require prohibitive compute time and memory.

`mlmm dft` extracts the ML region from the full enzyme PDB, appends link hydrogens, and runs a single-point PySCF (or GPU4PySCF) calculation. After the DFT evaluation, the script recomputes the **ML(dft)/MM total energy** by combining the PySCF high-level energy with MM evaluations of the full system (REAL-low) and the ML subset (MODEL-low):

```
E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low
```

The GPU4PySCF backend is activated automatically when available; otherwise PySCF CPU is used. The default functional/basis is `wb97m-v/def2-tzvpd`.

## Minimal example

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_dft
```

## Output checklist

- `result_dft/ml_region_with_linkH.xyz`
- `result_dft/result.yaml`
- Standard output block with ML(dft)/MM combined energy

## Common examples

1. Change functional/basis for a higher-level single point.

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --func-basis "wb97m-v/def2-tzvpd" --out-dir ./result_dft_tz
```

2. Freeze selected atoms in the ML/MM setup before DFT.

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q -1 -m 2 --freeze-atoms "1,3,5" --out-dir ./result_dft_freeze
```

3. Tighten SCF convergence and allow more cycles.

```bash
mlmm dft -i enzyme.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --conv-tol 1e-10 --max-cycle 200 --out-dir ./result_dft_tight
```

## Workflow

1. **Input handling** -- The full enzyme PDB (`-i`), Amber topology (`--parm`), and ML-region definition (`--model-pdb` or `--model-indices` or B-factor detection via `--detect-layer`) are loaded. Link hydrogens are appended automatically (C/N parents within 1.7 A) unless explicit `link_mlmm` pairs are provided via YAML.
2. **SCF build** -- `--func-basis` is parsed into functional and basis. Density fitting is enabled automatically with PySCF defaults. The GPU4PySCF backend is used when available; otherwise CPU PySCF is used.
3. **ML(dft)/MM recombination** -- After the DFT converges, MM evaluations of the full system (REAL-low) and the ML subset (MODEL-low) are computed. The combined energy is reported in Hartree and kcal/mol.
4. **Population analysis & outputs** -- Mulliken, meta-Lowdin, and IAO charges and spin densities (UKS only) are written alongside the combined energy block in `result.yaml`.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme PDB file (must be `.pdb`). | Required |
| `--parm PATH` | Amber parm7 topology for the full system. | Required |
| `--model-pdb PATH` | PDB defining the ML region (atom IDs must match the enzyme PDB). Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors (B=0/10/20). | `True` |
| `-q, --charge INT` | Charge of the ML region. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (e.g. `"1,3,5"`). Merged with YAML `geom.freeze_atoms`. | _None_ |
| `--func-basis TEXT` | Functional/basis pair as `"FUNC/BASIS"`. | `wb97m-v/def2-tzvpd` |
| `--max-cycle INT` | Maximum SCF iterations. | `100` |
| `--conv-tol FLOAT` | SCF convergence tolerance (Hartree). | `1e-9` |
| `--grid-level INT` | PySCF numerical integration grid level. | `3` |
| `--out-dir DIR` | Output directory. | `./result_dft/` |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved configuration and continue execution. | `False` |
| `--backend CHOICE` | MLIP backend used for the low-level ONIOM recombination: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running DFT. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |

## Outputs

```
out_dir/ (default: ./result_dft/)
├── ml_region_with_linkH.xyz    # ML-region coordinates (with link-H) used for DFT
├── result.yaml                 # DFT + ML(dft)/MM energy summary, charges, spin densities
└── (stdout)                    # Pretty-printed configuration blocks and energies
```

- `result.yaml` expands to:
  - `energy`: Hartree/kcal/mol values, convergence flag, wall time, backend info (gpu4pyscf vs pyscf(cpu)).
  - `charges`: Mulliken, meta-Lowdin, and IAO atomic charges (`null` when a method fails).
  - `spin_densities`: Mulliken, meta-Lowdin, and IAO spin densities (UKS-only for spins).
- It also summarizes charge, multiplicity, spin (2S), functional, basis, convergence knobs, and resolved output directory.

## YAML configuration

Accepts a mapping root; the `dft` section (and optional `geom`, `calc`/`mlmm`) is applied when present. Merge order is:
- defaults
- `--config`
- explicit CLI options

`dft` keys (defaults in parentheses):
- `func_basis` (`"wb97m-v/def2-tzvpd"`): Combined `FUNC/BASIS` string.
- `conv_tol` (`1e-9`): SCF convergence threshold (Hartree).
- `max_cycle` (`100`): Maximum SCF iterations.
- `grid_level` (`3`): PySCF `grids.level`.
- `verbose` (`4`): PySCF verbosity (0-9).
- `out_dir` (`"./result_dft/"`): Output directory root.

```yaml
geom:
 coord_type: cart                  # optional geom_loader settings
calc:
 charge: 0                         # ML region charge
 spin: 1                           # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
dft:
 func_basis: wb97m-v/def2-tzvpd      # exchange-correlation functional / basis set
 conv_tol: 1.0e-09                # SCF convergence tolerance (Hartree)
 max_cycle: 100                    # maximum SCF iterations
 grid_level: 3                     # PySCF grid level
 verbose: 4                        # PySCF verbosity (0-9)
 out_dir: ./result_dft/            # output directory root
```

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Link hydrogens are detected automatically (C/N parents within 1.7 A) unless explicit `link_mlmm` pairs are provided via YAML. Unsupported parent elements raise an error.
- The GPU4PySCF backend is activated automatically when available; otherwise PySCF CPU is used. If **Blackwell architecture** GPUs are detected, a warning is emitted because current GPU4PySCF may be unsupported.
- Density fitting is always attempted with PySCF defaults (no auxiliary basis guessing is implemented).
- DFT options (functional/basis, SCF controls) remain YAML-overridable under the `dft` key.
- IAO spin/charge analysis may fail for challenging systems; corresponding columns in `result.yaml` become `null` and a warning is printed.
- Exit codes: `0` when the SCF converges, `3` when it does not, `130` for keyboard interrupt, `1` on other errors.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [freq](freq.md) -- Vibrational frequency analysis (often precedes DFT refinement)
- [opt](opt.md) -- Single-structure geometry optimization
- [all](all.md) -- End-to-end workflow with `--dft`
- [YAML Reference](yaml_reference.md) -- Full `dft` configuration options
- [Glossary](glossary.md) -- Definitions of DFT, SP (Single Point)
