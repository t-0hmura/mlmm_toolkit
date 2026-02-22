# `tsopt`

## Overview

> **Summary:** Transition-state optimization using a partial-Hessian guided Dimer (`--opt-mode light`) or RS-I-RFO (`--opt-mode heavy`, default) workflow with the ML/MM calculator (FAIR-Chem UMA + hessian_ff).

`mlmm tsopt` carries out transition-state optimization tailored to the ML/MM calculator. The optimizer starts from a TS guess (e.g. a highest-energy image from `path-opt`/`path-search`, or your own structure) and refines it to a first-order saddle point.

### Key characteristics
- **Partial Hessian guided Dimer:** During the loose/final Dimer loops, the hessian_ff finite-difference Hessian is disabled (`mm_fd=False`). The UMA Hessian is embedded into the full 3N x 3N space with MM atoms zero-padded, providing a partial Hessian that still guides the Dimer direction updates.
- **Flatten loop with full Hessian:** Once the search enters the flatten loop, a full ML/MM Hessian (including the MM finite-difference block) is computed exactly once and then updated by Bofill steps in the active subspace between Dimer segments.
- **PHVA + TR projection:** Active-DOF projection and mass-weighted translation/rotation removal mirror `freq.py`, ensuring consistent imaginary-mode analysis and mode writing.

### Choosing `--opt-mode`
- **`heavy` (RS-I-RFO, default):** Conservative optimizer with full Hessian work. Generally more robust.
- **`light` (Dimer):** Lighter-weight search using partial Hessian guided Dimer. Often cheaper per step.

## Minimal example

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_tsopt
```

## Output checklist

- `result_tsopt/summary.md`
- `result_tsopt/key_ts.xyz` (or `key_ts.pdb`)
- `result_tsopt/key_imag_mode.trj`
- `result_tsopt/final_geometry.pdb` (or `final_geometry.xyz`)
- `result_tsopt/vib/final_imag_mode_*.trj`
- `result_tsopt/vib/final_imag_mode_*.pdb`

## Common examples

1. Use light mode with analytical Hessian when VRAM is sufficient.

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. Keep a full optimization trajectory for inspection.

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. Run heavy mode with YAML overrides.

```bash
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

## Usage
```bash
mlmm tsopt -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb \
    -q CHARGE [-m SPIN] [--freeze-atoms "1,3,5"] [--max-cycles N] \
    [--dump/--no-dump] [--out-dir DIR] [--thresh PRESET] \
    [--opt-mode light|heavy] [--hessian-calc-mode Analytical|FiniteDifference] \
    [--config FILE] [--override-yaml FILE] [--show-config] [--dry-run] [--args-yaml FILE]
```

### Examples
```bash
# Default heavy mode (RS-I-RFO)
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --max-cycles 8000 --dump --out-dir ./result_tsopt/

# Light mode (Dimer) with analytical Hessian
mlmm tsopt -i ts_guess.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt/
```

## Workflow
1. **Input handling** -- Load the enzyme PDB, Amber topology, and ML-region definition. Resolve charge/spin. Freeze atoms from CLI and YAML are merged.
2. **ML/MM calculator setup** -- Build the ML/MM calculator (FAIR-Chem UMA + hessian_ff). The `--hessian-calc-mode` controls whether UMA evaluates Hessians analytically or via finite difference.
3. **Light mode (Dimer):**
   - The Hessian Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian (active subspace, TR-projected).
   - When the flatten loop is enabled, the stored active Hessian is updated via Bofill using displacements and gradient differences. Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, and runs a dimer + LBFGS micro-segment.
4. **Heavy mode (RS-I-RFO):**
   - Runs the RS-I-RFO optimizer with optional Hessian reference files and micro-cycle controls defined in the `rsirfo` YAML section.
5. **Mode export** -- The converged imaginary mode is written to `vib/` as `.trj`/`.pdb` pairs. Final geometry and optional trajectory are also saved.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Starting geometry (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--real-parm7 PATH` | Amber parm7 topology for the whole enzyme. | Required |
| `--model-pdb PATH` | PDB containing the ML-region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | Total charge of the ML region. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | _None_ |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (A) for MM Hessian atoms. Providing cutoffs disables `--detect-layer`. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (A) for movable MM atoms. | _None_ |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode: `Analytical` or `FiniteDifference`. | _None_ |
| `--max-cycles INT` | Maximum total optimizer cycles. | `10000` |
| `--dump/--no-dump` | Write concatenated trajectory `optimization_all.trj`. | `False` |
| `--out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Convergence preset (`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`). | _None_ |
| `--opt-mode CHOICE` | TS optimizer mode: `light` (Dimer) or `heavy` (RS-I-RFO). | `heavy` |
| `--partial-hessian-flatten / --full-hessian-flatten` | Use partial Hessian (ML only) for imaginary mode detection in flatten loop. | `True` (partial) |
| `--active-dof-mode CHOICE` | Active DOF for final frequency analysis: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--override-yaml FILE` | Final YAML override file (highest-priority YAML layer). | _None_ |
| `--args-yaml FILE` | Legacy alias of `--override-yaml`. | _None_ |
| `--show-config/--no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--dry-run/--no-dry-run` | Validate inputs/config and print the execution plan without running TS optimization. | `False` |

## Outputs
```
out_dir/  (default: ./result_tsopt/)
├── summary.md                   # Quick index of key outputs
├── key_ts.xyz                   # Shortcut to final TS geometry (or key_ts.pdb)
├── key_imag_mode.trj            # Shortcut to a representative imaginary mode
├── key_opt.trj                  # Shortcut to optimization trajectory (when available)
├── final_geometry.xyz            # Always written
├── final_geometry.pdb            # When the input was PDB
├── optimization_all.trj          # Concatenated Dimer segments (when --dump)
├── optimization_all.pdb          # PDB companion (when --dump and input was PDB)
├── vib/
│   ├── final_imag_mode_±XXXX.Xcm-1.trj   # Imaginary mode trajectory
│   └── final_imag_mode_±XXXX.Xcm-1.pdb   # Imaginary mode PDB companion
└── .dimer_mode.dat               # Dimer orientation seed (light mode)
```

## YAML configuration (`--config` / `--override-yaml` / `--args-yaml`)

Use `--config` for the base mapping and `--override-yaml` for the final YAML override (`--args-yaml` is a legacy alias of `--override-yaml`). Merge precedence is:

`defaults < config < explicit CLI < override`.

```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
opt:
  thresh: baker
  max_cycles: 10000
  dump: false
  out_dir: ./result_tsopt/
hessian_dimer:
  thresh_loose: gau_loose
  thresh: baker
  update_interval_hessian: 500
  neg_freq_thresh_cm: 5.0
  flatten_amp_ang: 0.1
  flatten_max_iter: 50
  root: 0
  dimer:
    length: 0.0189
    rotation_max_cycles: 15
  lbfgs:
    thresh: baker
    max_cycles: 10000
    max_step: 0.3
rsirfo:
  thresh: baker
  max_cycles: 10000
  roots: [0]
  hessian_update: bofill
```

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Imaginary-mode detection uses a default threshold of ~5 cm^-1 (configurable via `hessian_dimer.neg_freq_thresh_cm`). The selected `root` determines which imaginary mode is exported.
- `--freeze-atoms` accepts 1-based indices and is merged with YAML `geom.freeze_atoms`.
- Convergence presets propagate to both the outer bookkeeping (`opt`) and the inner LBFGS segments (`hessian_dimer.lbfgs`).
- PHVA translation/rotation projection mirrors the implementation in `freq`, reducing GPU memory consumption while preserving correct eigenvectors in the active space.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing

- [opt](opt.md) -- Single-structure geometry optimization
- [freq](freq.md) -- Confirm a single imaginary frequency for the validated TS
- [irc](irc.md) -- Trace the reaction path from an optimized TS
- [all](all.md) -- End-to-end workflow that chains extraction -> MEP -> tsopt -> IRC -> freq
