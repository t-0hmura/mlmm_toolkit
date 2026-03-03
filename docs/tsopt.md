# `tsopt`

## Overview

> **Summary:** Optimize a transition-state *candidate* using Dimer (`--opt-mode grad`) or RS-I-RFO (`--opt-mode hess`, default). Microiteration (`--microiter`, default on) alternates ML 1-step RS-I-RFO and MM relaxation in `hess` mode. A validated TS should show **exactly one** imaginary frequency; always confirm the mode/connectivity with freq/IRC.

### At a glance
- **Use when:** You have a TS guess (HEI from `path-opt`/`path-search`, or your own structure) and want to refine it to a first-order saddle point with ML/MM.
- **Method:** `hess` = RS-I-RFO (default, generally more robust). `grad` = Hessian Guided Dimer (often cheaper per step). Aliases `heavy`/`light` and `rsirfo`/`dimer` are also accepted.
- **Outputs:** `final_geometry.xyz`/`.pdb`, imaginary-mode animations in `vib/`.
- **Next step:** Run [freq](freq.md) to confirm exactly one imaginary frequency, then [irc](irc.md) to verify connectivity.

### Choosing `--opt-mode`
- Use **`--opt-mode hess` (RS-I-RFO)** when you want the default, conservative optimizer and you can afford Hessian work. With `--microiter` (default on), ML and MM regions are optimized alternately.
- Use **`--opt-mode grad` (Dimer)** when you want a lighter-weight search, or when you plan to iterate quickly from several TS guesses. `--ml-only-hessian-dimer` uses only the ML-region Hessian for dimer orientation (faster but less accurate).

`mlmm tsopt` carries out transition-state optimization tailored to the ML/MM calculator. The optimizer starts from a TS guess and refines it to a first-order saddle point.

### Key characteristics
- **Partial Hessian guided Dimer:** During the loose/final Dimer loops, the hessian_ff finite-difference Hessian is disabled (`mm_fd=False`). The UMA Hessian is embedded into the full 3N x 3N space with MM atoms zero-padded, providing a partial Hessian that still guides the Dimer direction updates.
- **Flatten loop with full Hessian:** Once the search enters the flatten loop, a full ML/MM Hessian (including the MM finite-difference block) is computed exactly once and then updated by Bofill steps in the active subspace between Dimer segments.
- **PHVA + TR projection:** Active-DOF projection and mass-weighted translation/rotation removal mirror `freq.py`, ensuring consistent imaginary-mode analysis and mode writing.
- **Output conversion:** With `--convert-files` (default), PDB inputs can be mirrored to `.pdb` (when `--dump`), and the imaginary mode is exported as `.pdb` alongside `_trj.xyz`.

## Minimal example

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_tsopt
```

## Output checklist

- `result_tsopt/final_geometry.pdb` (or `final_geometry.xyz`)
- `result_tsopt/vib/final_imag_mode_*_trj.xyz`
- `result_tsopt/vib/final_imag_mode_*.pdb`

## Common examples

1. Use light mode with analytical Hessian when VRAM is sufficient.

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode light --hessian-calc-mode Analytical --out-dir ./result_tsopt_light
```

2. Keep a full optimization trajectory for inspection.

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --dump --out-dir ./result_tsopt_dump
```

3. Run heavy mode with YAML overrides.

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --opt-mode heavy --config tsopt.yaml --out-dir ./result_tsopt_heavy
```

## Workflow

1. **Input handling** -- Load the enzyme PDB, Amber topology, and ML-region definition. Resolve charge/spin. Freeze atoms from CLI and YAML are merged.
2. **ML/MM calculator setup** -- Build the ML/MM calculator (FAIR-Chem UMA + hessian_ff). The `--hessian-calc-mode` controls whether UMA evaluates Hessians analytically or via finite difference.
4. **Light mode (Dimer):**
   - The Hessian Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian (active subspace, TR-projected).
   - When the flatten loop is enabled (`--flatten`), the stored active Hessian is updated via Bofill using displacements and gradient differences. Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, and runs a dimer + LBFGS micro-segment.
5. **Heavy mode (RS-I-RFO):**
   - Runs the RS-I-RFO optimizer with optional Hessian reference files and micro-cycle controls defined in the `rsirfo` YAML section.
   - When `--flatten` is enabled and more than one imaginary mode remains after convergence, the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten iteration cap is reached.
6. **Mode export & conversion** -- The converged imaginary mode is always written to `vib/final_imag_mode_*_trj.xyz` and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimization trajectory and final geometry are also converted to PDB via the input template when `--dump`.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Starting geometry (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the whole enzyme. | Required |
| `--model-pdb PATH` | PDB containing the ML-region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | Total charge of the ML region. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--hess-cutoff FLOAT` | Distance cutoff (A) for MM Hessian atoms. Providing cutoffs disables `--detect-layer`. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (A) for movable MM atoms. | _None_ |
| `--hessian-calc-mode CHOICE` | UMA Hessian mode: `Analytical` or `FiniteDifference`. | _None_ |
| `--max-cycles INT` | Maximum total optimizer cycles. | `10000` |
| `--opt-mode CHOICE` | TS optimizer mode: `grad` (Dimer) or `hess` (RS-I-RFO). Aliases `light`/`heavy` and `dimer`/`rsirfo` accepted. | `hess` |
| `--microiter/--no-microiter` | Microiteration: alternate ML 1-step (RS-I-RFO) + MM relaxation (LBFGS). Only effective in `hess` mode. | `True` |
| `--ml-only-hessian-dimer/--no-ml-only-hessian-dimer` | Use ML-region-only Hessian for dimer orientation in `grad` mode. Faster but less accurate. | `False` |
| `--flatten/--no-flatten` | Enable the extra-imaginary-mode flattening loop. Applies to both light and heavy modes. | `False` |
| `--dump/--no-dump` | Write concatenated trajectory `optimization_all_trj.xyz`. | `False` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions for PDB inputs. | `True` |
| `--out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Convergence preset (`gau_loose\|gau\|gau_tight\|gau_vtight\|baker\|never`). | _None_ |
| `--partial-hessian-flatten / --full-hessian-flatten` | Use partial Hessian (ML only) for imaginary mode detection in flatten loop. | `True` (partial) |
| `--active-dof-mode CHOICE` | Active DOF for final frequency analysis: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--config FILE` | Base YAML configuration file applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--dry-run/--no-dry-run` | Validate inputs/config and print the execution plan without running TS optimization. | `False` |

## Outputs

```
out_dir/ (default: ./result_tsopt/)
├── final_geometry.xyz             # Always written
├── final_geometry.pdb             # When the input was PDB
├── optimization_all_trj.xyz       # Concatenated Dimer segments (when --dump)
├── optimization_all.pdb           # PDB companion (when --dump and input was PDB)
├── vib/
│   ├── final_imag_mode_±XXXX.Xcm-1_trj.xyz  # Imaginary mode trajectory
│   └── final_imag_mode_±XXXX.Xcm-1.pdb      # Imaginary mode PDB companion
└── .dimer_mode.dat                # Dimer orientation seed (light mode)
```

## YAML configuration

Settings are applied with **defaults < config < explicit CLI < override**.
Shared sections reuse [YAML Reference](yaml_reference.md). Keep the full block below intact if it already matches your workflow -- adjust only the values you need to change.

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
 uma_model: uma-s-1p1              # UMA model tag
 uma_task_name: omol                # UMA task name
 ml_device: auto                   # UMA device selection
 ml_hessian_mode: FiniteDifference  # Hessian mode selection
opt:
 thresh: baker                     # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000                 # optimizer cycle cap
 print_every: 100                  # logging stride
 min_step_norm: 1.0e-08            # minimum norm for step acceptance
 assert_min_step: true             # stop if steps fall below threshold
 rms_force: null                   # explicit RMS force target
 rms_force_only: false             # rely only on RMS force convergence
 max_force_only: false             # rely only on max force convergence
 force_only: false                 # skip displacement checks
 converge_to_geom_rms_thresh: 0.05  # geom RMS threshold when converging to ref
 overachieve_factor: 0.0           # factor to tighten thresholds
 check_eigval_structure: false     # validate Hessian eigenstructure
 line_search: true                 # enable line search
 dump: false                       # dump trajectory/restart data
 dump_restart: false               # dump restart checkpoints
 prefix: ""                        # filename prefix
 out_dir: ./result_tsopt/          # output directory
hessian_dimer:
 thresh_loose: gau_loose           # loose convergence preset
 thresh: baker                     # main convergence preset
 update_interval_hessian: 500      # Hessian rebuild cadence
 neg_freq_thresh_cm: 5.0           # negative frequency threshold (cm^-1)
 flatten_amp_ang: 0.1              # flattening amplitude (A)
 flatten_max_iter: 50              # flattening iteration cap (disabled when --no-flatten)
 flatten_sep_cutoff: 0.0           # minimum distance between representative atoms (A)
 flatten_k: 10                     # representative atoms sampled per mode
 flatten_loop_bofill: false        # Bofill update for flatten displacements
 mem: 100000                       # memory limit for solver
 device: auto                      # device selection for eigensolver
 root: 0                           # targeted TS root index
 dimer:
  length: 0.0189                   # dimer separation (Bohr)
  rotation_max_cycles: 15          # max rotation iterations
  rotation_method: fourier         # rotation optimizer method
  rotation_thresh: 0.0001          # rotation convergence threshold
  rotation_tol: 1                  # rotation tolerance factor
  rotation_max_element: 0.001      # max rotation matrix element
  rotation_interpolate: true       # interpolate rotation steps
  rotation_disable: false          # disable rotations entirely
  rotation_disable_pos_curv: true  # disable when positive curvature detected
  rotation_remove_trans: true      # remove translational components
  trans_force_f_perp: true         # project forces perpendicular to translation
  bonds: null                      # bond list for constraints
  N_hessian: null                  # Hessian size override
  bias_rotation: false             # bias rotational search
  bias_translation: false          # bias translational search
  bias_gaussian_dot: 0.1           # Gaussian bias dot product
  seed: null                       # RNG seed for rotations
  write_orientations: true         # write rotation orientations
  forward_hessian: true            # propagate Hessian forward
 lbfgs:
  thresh: baker                    # LBFGS convergence preset
  max_cycles: 10000                # iteration limit
  print_every: 100                 # logging stride
  min_step_norm: 1.0e-08           # minimum accepted step norm
  assert_min_step: true            # assert when steps stagnate
  max_step: 0.3                    # maximum step length
  control_step: true               # control step length adaptively
  double_damp: true                # double damping safeguard
  keep_last: 7                     # history size for LBFGS buffers
  beta: 1.0                        # initial damping beta
  mu_reg: null                     # regularization strength
  max_mu_reg_adaptions: 10         # cap on mu adaptations
rsirfo:
 thresh: baker                     # RS-IRFO convergence preset
 max_cycles: 10000                 # iteration cap
 print_every: 100                  # logging stride
 min_step_norm: 1.0e-08            # minimum accepted step norm
 assert_min_step: true             # assert when steps stagnate
 roots: [0]                        # target root indices
 hessian_ref: null                 # reference Hessian
 hessian_update: bofill            # Hessian update scheme override
 hessian_recalc_reset: true        # reset recalc counter after exact Hessian
 max_micro_cycles: 50              # micro-iterations per macro cycle
 augment_bonds: false              # augment reaction path based on bond analysis
 min_line_search: true             # enforce minimum line-search step
 max_line_search: true             # enforce maximum line-search step
 assert_neg_eigval: false          # require a negative eigenvalue at convergence
```

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Imaginary-mode detection uses a default threshold of ~5 cm^-1 (configurable via `hessian_dimer.neg_freq_thresh_cm`). The selected `root` determines which imaginary mode is exported.
- `--freeze-atoms` accepts 1-based indices and is merged with YAML `geom.freeze_atoms`.
- Convergence presets propagate to both the outer bookkeeping (`opt`) and the inner LBFGS segments (`hessian_dimer.lbfgs`).
- PHVA translation/rotation projection mirrors the implementation in `freq`, reducing GPU memory consumption while preserving correct eigenvectors in the active space.
- `return_partial_hessian` is partial-first in `tsopt`: when YAML does not explicitly set `calc.return_partial_hessian`, partial Hessian is used. Set it to `false` explicitly to force full Hessian.
- Config merge precedence is `defaults < config < explicit CLI < override`.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [opt](opt.md) -- Single-structure geometry optimization
- [freq](freq.md) -- Confirm a single imaginary frequency for the validated TS
- [irc](irc.md) -- Trace the reaction path from an optimized TS
- [all](all.md) -- End-to-end workflow that chains extraction -> MEP -> tsopt -> IRC -> freq
- [YAML Reference](yaml_reference.md) -- Full `hessian_dimer` (Hessian Guided Dimer) and `rsirfo` configuration options
- [Glossary](glossary.md) -- Definitions of TS, Dimer, RS-I-RFO, Hessian
