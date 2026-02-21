# `freq`

## Overview

> **Summary:** Compute ML/MM vibrational frequencies and thermochemistry (ZPE, Gibbs energy, etc.) with PHVA support. Imaginary frequencies appear as negative values; a properly converged TS should have exactly one.

### Quick reference
- **Use when:** You want to validate a minimum/TS candidate and/or compute thermo corrections from ML/MM.
- **Frozen atoms:** Supported via PHVA (partial Hessian vibrational analysis).
- **Outputs:** `frequencies_cm-1.txt`, per-mode `.trj` and `.pdb` animations, plus `thermoanalysis.yaml` when enabled.
- **TS check:** A properly converged TS is expected to have **exactly one** imaginary frequency (negative cm^-1).

`mlmm freq` performs vibrational analysis with the ML/MM calculator (`mlmm_toolkit.mlmm_calc.mlmm`), honoring frozen atoms via PHVA. It exports normal-mode animations as `.trj` and `.pdb` (mapped back onto the enzyme ordering), and prints a Gaussian-style thermochemistry summary when the optional `thermoanalysis` package is installed.

Configuration follows **defaults < `--config` < explicit CLI < `--override-yaml`** (`geom`, `calc`/`mlmm`, `freq`, `thermo`). Legacy `--args-yaml` is still accepted as an alias of `--override-yaml`. The ML region is supplied via `--model-pdb`, the full enzyme topology via `--real-parm7`, and the input coordinates must be a full enzyme PDB (no link atoms).

## Minimal example

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_freq
```

## Output checklist

- `result_freq/summary.md`
- `result_freq/key_frequencies.txt`
- `result_freq/key_mode_1.trj`
- `result_freq/frequencies_cm-1.txt`
- `result_freq/mode_*.trj`
- `result_freq/mode_*.pdb` (for PDB inputs)

## Common examples

1. Limit the number of exported modes for quick inspection.

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

2. Run PHVA with explicit frozen atoms and dump thermo payload.

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

3. Use analytical Hessian mode on VRAM-rich nodes.

```bash
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Usage
```bash
mlmm freq -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m MULT] [--freeze-atoms "1,3,5"] \
    [--max-write N] [--amplitude-ang FLOAT] [--n-frames N] [--sort {value|abs}] \
    [--temperature K] [--pressure atm] [--dump/--no-dump] \
    [--hessian-calc-mode {Analytical|FiniteDifference}] \
    [--active-dof-mode {all|ml-only|partial|unfrozen}] \
    [--out-dir DIR] [--config FILE] [--override-yaml FILE|--args-yaml FILE] \
    [--show-config] [--dry-run]
```

### Examples
```bash
# Minimal run
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# PHVA with custom options
mlmm freq -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --freeze-atoms "1,3,5,7" --max-write 10 --sort abs --dump --config ./freq.yaml
```

## Workflow
- **Geometry loading & freeze handling**: structures are read via
  `pysisyphus.helpers.geom_loader`. `--freeze-atoms "1,3,5"` accepts 1-based indices and
  merges them with YAML `geom.freeze_atoms`; the merged list is passed to both the
  geometry echo and the ML/MM calculator, enabling PHVA.
- **ML/MM calculator**: the ML region is supplied via `--model-pdb`; Amber parameters are
  read from `--real-parm7`. `--hessian-calc-mode` selects analytical or finite-difference
  Hessians. The calculator may return either the full 3N x 3N Hessian or an active-DOF
  sub-block.
- **PHVA & TR projection**: with frozen atoms, eigenanalysis occurs inside the active
  subspace with translation/rotation modes projected there. Both 3N x 3N and active-block
  Hessians are accepted, and frequencies are reported in cm^-1 (negatives = imaginary).
- **Active DOF mode**: `--active-dof-mode` controls which atoms are included in the
  frequency analysis: `all` (all atoms), `ml-only` (ML layer, B=0),
  `partial` (ML + Hessian-target MM; default), `unfrozen` (non-frozen layers, typically B=0/10).
- **Mode export**: `--max-write` limits how many modes are animated. Modes are sorted by
  value (or absolute value with `--sort abs`). Each exported mode writes `.trj` (XYZ-like
  trajectory) and `.pdb` files (PDB animation mapped back onto the enzyme ordering).
  The sinusoidal animation amplitude (`--amplitude-ang`) and frame count (`--n-frames`)
  match the YAML defaults.
- **Thermochemistry**: if `thermoanalysis` is installed, a QRRHO-like summary (EE, ZPE,
  E/H/G corrections, heat capacities, entropies) is printed using PHVA frequencies.
  CLI pressure in atm is converted internally to Pa. When `--dump`, a
  `thermoanalysis.yaml` snapshot is also written.
- **Device selection**: `ml_device="auto"` triggers CUDA when available, otherwise CPU.
  The internal TR projection/mode assembly runs on the same device to minimise transfers.
- **Exit behavior**: keyboard interrupts exit with code 130; other failures print a
  traceback and exit with code 1.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Full enzyme PDB (no link atoms). | Required |
| `--real-parm7 PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region. | Required |
| `--model-indices TEXT` | Explicit ML-region atom indices (alternative to `--model-pdb`). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Indexing convention for `--model-indices`. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors. | `True` |
| `-q, --charge INT` | ML region charge. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | 1-based comma-separated frozen atom indices. | _None_ |
| `--hess-cutoff FLOAT` | Cutoff distance for Hessian-target MM atoms. | _None_ |
| `--movable-cutoff FLOAT` | Cutoff distance for movable-MM layer. | _None_ |
| `--hessian-calc-mode CHOICE` | Hessian mode (`Analytical` or `FiniteDifference`). | _None_ |
| `--max-write INT` | Number of modes to export. | `20` |
| `--amplitude-ang FLOAT` | Mode animation amplitude (angstrom). | `0.8` |
| `--n-frames INT` | Frames per mode animation. | `20` |
| `--sort CHOICE` | Mode ordering: `value` (cm^-1) or `abs`. | `value` |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_freq/` |
| `--active-dof-mode CHOICE` | Active DOF selection: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--override-yaml FILE` | Final YAML override (highest-priority YAML layer). | _None_ |
| `--args-yaml FILE` | Legacy alias of `--override-yaml`. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. | `False` |
| `--ref-pdb FILE` | Reference PDB topology for non-PDB inputs. | _None_ |

## Outputs
```
out_dir/ (default: ./result_freq/)
  summary.md                      # Quick index of key outputs
  key_frequencies.txt             # Shortcut to frequencies_cm-1.txt
  key_mode_1.trj                  # Shortcut to a representative mode trajectory
  key_mode_1.pdb                  # Shortcut to representative mode PDB (when available)
  key_thermo.yaml                 # Shortcut to thermoanalysis.yaml (when available)
  mode_XXXX_{+/-freq}cm-1.trj    # XYZ-like trajectory, sinusoidal animation per mode
  mode_XXXX_{+/-freq}cm-1.pdb    # PDB animation mapped back onto the enzyme ordering
  frequencies_cm-1.txt            # All computed frequencies (cm^-1) sorted by the chosen key
  thermoanalysis.yaml             # Optional thermochemistry payload when --dump
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- The ML/MM calculator returns Hessians in Hartree/Bohr^2. Conversions to cm^-1 follow the
  PySisyphus/ASE conventions used elsewhere in the toolkit.
- Thermochemistry relies on the optional `thermoanalysis` package; when absent, only a
  warning is printed and execution continues.
- `--hessian-calc-mode` follows the standard precedence (defaults < config < explicit CLI < override).

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing

- [tsopt](tsopt.md) -- Optimize TS candidates (validate with freq/IRC; expected: one imaginary frequency)
- [opt](opt.md) -- Geometry optimization (often precedes freq)
- [dft](dft.md) -- Single-point DFT for higher-level energy refinement
- [all](all.md) -- End-to-end workflow with `--thermo`
- [YAML Reference](yaml-reference.md) -- Full `freq` and `thermo` configuration options
