# `freq`

## Overview

> **Summary:** Compute ML/MM vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) with PHVA support. When VRAM permits, `--hessian-calc-mode Analytical` speeds Hessian evaluation. Imaginary frequencies appear as negative values.

### At a glance
- **Use when:** You want to validate a minimum/TS candidate and/or compute thermo corrections from ML/MM.
- **Method:** Full or partial Hessian vibrational analysis (PHVA) with the ML/MM calculator (MLIP backend + hessian_ff). Backend selected via `--backend` (default: `uma`).
- **Outputs:** `frequencies_cm-1.txt`, per-mode `_trj.xyz` and `.pdb` animations, plus `thermoanalysis.yaml` when enabled.
- **Next step:** Use results to confirm a minimum (no imaginary frequencies) or a TS (exactly one imaginary frequency).

`mlmm freq` performs vibrational analysis with the ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`), honoring frozen atoms via PHVA. It exports normal-mode animations as `_trj.xyz` and `.pdb` (mapped back onto the enzyme ordering), and prints a Gaussian-style thermochemistry summary when the optional `thermoanalysis` package is installed.

## Minimal example

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

## Output checklist

- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*_trj.xyz`
- `result_freq/mode_*.pdb` (for PDB inputs)

## Common examples

1. Limit the number of exported modes for quick inspection.

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. Run PHVA with explicit frozen atoms and dump thermo payload.

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

3. Use analytical Hessian mode on VRAM-rich nodes.

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Workflow

1. **ML/MM calculator setup** -- The ML region is supplied via `--model-pdb`; Amber parameters are read from `--parm`. `--hessian-calc-mode` selects analytical or finite-difference Hessians. The calculator may return either the full 3N x 3N Hessian or an active-DOF sub-block.
   - When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
2. **PHVA & TR projection** -- With frozen atoms, eigenanalysis occurs inside the active subspace with translation/rotation modes projected there. Both 3N x 3N and active-block Hessians are accepted, and frequencies are reported in cm^-1 (negatives = imaginary).
3. **Active DOF mode** -- `--active-dof-mode` controls which atoms are included in the frequency analysis: `all` (all atoms), `ml-only` (ML layer, B=0), `partial` (ML + MovableMM; default), `unfrozen` (non-frozen layers, typically B=0/10).
4. **Mode export** -- `--max-write` limits how many modes are animated. Modes are sorted by value (or absolute value with `--sort abs`). Each exported mode writes `_trj.xyz` (XYZ-like trajectory) and `.pdb` files (PDB animation mapped back onto the enzyme ordering). The sinusoidal animation amplitude (`--amplitude-ang`) and frame count (`--n-frames`) match the YAML defaults.
5. **Thermochemistry** -- If `thermoanalysis` is installed, a QRRHO-like summary (EE, ZPE, E/H/G corrections, heat capacities, entropies) is printed using PHVA frequencies. CLI pressure in atm is converted internally to Pa. When `--dump`, a `thermoanalysis.yaml` snapshot is also written.
6. **Device selection** -- `ml_device="auto"` triggers CUDA when available, otherwise CPU. The internal TR projection/mode assembly runs on the same device to minimize transfers.
7. **Exit behavior** -- Keyboard interrupts exit with code 130; other failures print a traceback and exit with code 1.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme PDB (no link atoms). | Required |
| `--parm PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Explicit ML-region atom indices (alternative to `--model-pdb`). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Indexing convention for `--model-indices`. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors. | `True` |
| `-q, --charge INT` | ML region charge. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | 1-based comma-separated frozen atom indices. | _None_ |
| `--hess-cutoff FLOAT` | Cutoff distance for MovableMM atoms. | _None_ |
| `--movable-cutoff FLOAT` | Cutoff distance for movable-MM layer. | _None_ |
| `--hessian-calc-mode CHOICE` | Hessian mode (`Analytical` or `FiniteDifference`). | _None_ |
| `--max-write INT` | Number of modes to export. | `10` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (angstrom). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm^-1) or `abs`. | `value` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |
| `-o, --out-dir TEXT` | Output directory. | `./result_freq/` |
| `--active-dof-mode CHOICE` | Active DOF selection: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--hess-device CHOICE` | Device for Hessian assembly/diagonalization: `auto`, `cuda`, `cpu`. Use `cpu` to avoid VRAM issues with large systems. | `auto` |
| `--ref-pdb FILE` | Reference PDB topology for non-PDB inputs. | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. Shown in `--help-advanced`. | `False` |

## Outputs

```
out_dir/ (default: ./result_freq/)
├─ mode_XXXX_±freqcm-1_trj.xyz   # Per-mode animations (XYZ-like trajectory)
├─ mode_XXXX_±freqcm-1.pdb       # PDB animation mapped back onto the enzyme ordering
├─ frequencies_cm-1.txt           # Full frequency list using the selected sort order
└─ thermoanalysis.yaml            # Present when thermoanalysis is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## YAML configuration

Provide mappings with merge order **defaults < config < explicit CLI < override**.
Shared sections reuse [YAML Reference](yaml_reference.md).
An additional `thermo` section is supported for thermochemistry controls.

```yaml
geom:
 coord_type: cart                  # coordinate type: cartesian vs dlc internals
 freeze_atoms: []                  # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0                         # total charge (CLI override)
 spin: 1                           # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
 backend: uma                      # MLIP backend: uma | orb | mace | aimnet2
 embedcharge: false                # xTB point-charge embedding correction
 uma_model: uma-s-1p2              # uma-s-1p1 | uma-s-1p2 | uma-m-1p1
 uma_task_name: omol                # UMA task name (UMA backend only)
 ml_device: auto                   # ML backend device selection
 ml_hessian_mode: Analytical         # Hessian mode selection
 out_hess_torch: true              # request torch-form Hessian
 mm_fd: true                       # MM finite-difference toggle
 return_partial_hessian: true      # allow partial Hessians (PHVA default)
freq:
 amplitude_ang: 0.8                # displacement amplitude for modes (Å)
 n_frames: 20                      # number of frames per mode
 max_write: 10                     # maximum number of modes to write
 sort: value                       # sort order: value vs abs
thermo:
 temperature: 298.15               # thermochemistry temperature (K)
 pressure_atm: 1.0                 # thermochemistry pressure (atm)
 dump: false                       # write thermoanalysis.yaml when true
```

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The ML/MM calculator returns Hessians in Hartree/Bohr^2. Conversions to cm^-1 follow the PySisyphus/ASE conventions used elsewhere in the toolkit.
- Thermochemistry relies on the optional `thermoanalysis` package; when absent, only a warning is printed and execution continues.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI < override).
- `freq` is partial-first: unless YAML explicitly sets `calc.return_partial_hessian`, the calculation uses partial Hessian (PHVA-oriented path) by default.
- Imaginary modes are reported as negative frequencies. `freq` prints how many were detected and dumps details when `--dump`.
- **PHVA details:** With frozen atoms, the active subspace is formed by excluding frozen DOF. Translation/rotation modes are projected within this subspace, ensuring correct mode shapes. Both full 3N x 3N and active-block Hessians are accepted. The mass-weighted eigendecomposition uses `UPLO="U"` to minimize GPU memory.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [tsopt](tsopt.md) -- Optimize TS candidates (validate with freq/IRC; expected: one imaginary frequency)
- [opt](opt.md) -- Geometry optimization (often precedes freq)
- [dft](dft.md) -- Single-point DFT for higher-level energy refinement
- [all](all.md) -- End-to-end workflow with `--thermo`
- [YAML Reference](yaml_reference.md) -- Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) -- Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
