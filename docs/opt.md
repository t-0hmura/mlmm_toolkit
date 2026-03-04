# `opt`

## Overview

> **Summary:** Optimizes a single structure to a local minimum using L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). Optional imaginary-mode flattening can be enabled with `--flatten`. Microiteration (`--microiter`, default on) alternates ML 1-step and MM relaxation in `hess` mode.

`mlmm opt` optimizes a single structure to a local minimum using L-BFGS (`--opt-mode grad`, default) or RFO (`--opt-mode hess`). Aliases `light`/`heavy` and `lbfgs`/`rfo` are also accepted. The ML/MM calculator (MLIP backend + hessian_ff) provides energies, gradients, and Hessians. The MLIP backend is selected via `--backend` (default: `uma`; choices: `uma`, `orb`, `mace`, `aimnet2`). Input structures can be `.pdb`, `.xyz`, `_trj.xyz`, or any format supported by `geom_loader`. Settings follow precedence: **defaults < config < explicit CLI < override**.

When the starting structure is a PDB, the command also writes `.pdb` companions, controlled by `--convert-files/--no-convert-files` (enabled by default). PDB-specific conveniences include:
- Output conversion produces `final_geometry.pdb` (and `optimization.pdb` when dumping trajectories) using the input PDB as the topology reference.
- B-factors are annotated as: ML-region atoms = 100.00, frozen atoms = 50.00, atoms that are both ML and frozen = 150.00.

### At a glance
- **Use when:** You want to minimize a single enzyme structure to a local energy minimum with ML/MM.
- **Method:** L-BFGS (grad, default) or RFO (hess). Aliases `light`/`heavy` and `lbfgs`/`rfo` are also accepted. In `hess` mode, microiteration (default on) alternates ML 1-step RFO with MM LBFGS relaxation.
- **Outputs:** `final_geometry.xyz`, `final_geometry.pdb` (PDB inputs), optional trajectory.
- **Next step:** Run [freq](freq.md) to confirm the structure is a true minimum (no imaginary frequencies).

## Minimal example

```bash
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_opt
```

## Output checklist

- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb` (when the input is PDB and conversion is enabled)
- `result_opt/optimization_trj.xyz` (when `--dump` is enabled)

## Common examples

1. Tighten convergence and keep an optimization trajectory.

```bash
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --thresh gau_tight --dump --out-dir ./result_opt_tight
```

2. Apply one harmonic distance restraint during optimization.

```bash
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --dist-freeze "[(12,45,2.20)]" --bias-k 20.0 --out-dir ./result_opt_rest
```

3. Switch to heavy mode (RFO).

```bash
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --opt-mode heavy --out-dir ./result_opt_rfo
```

4. Use the ORB backend instead of the default UMA.

```bash
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --backend orb --out-dir ./result_opt_orb
```

## Workflow

1. **Input handling** -- The tool requires `-i/--input` to be a PDB file (the enzyme complex). The optimizer reads coordinates from this PDB via `pysisyphus.helpers.geom_loader`. ML/MM layer definitions come from `--model-pdb`, `--model-indices`, or `--detect-layer` (B-factor encoding: B=0 ML, B=10 Hessian-target MM, B=20 frozen MM).
2. **ML/MM calculator setup** -- Build the ML/MM calculator (MLIP backend + hessian_ff). The `--backend` option selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). `--parm` provides Amber MM topology; `--model-pdb` defines the ML region. When `--embedcharge` is enabled, xTB point-charge embedding is applied to correct for MM-to-ML environmental electrostatic effects.
4. **Optimization** -- `--opt-mode light` runs L-BFGS and `--opt-mode heavy` runs RFOptimizer (RFO).
   - `--flatten` enables post-optimization flattening of imaginary modes. All detected imaginary modes are flattened each iteration until none remain or the internal loop cap is reached.
5. **Restraints** -- `--dist-freeze` consumes Python-literal tuples `(i, j, target_A)` where `target_A` is the target distance in angstrom; omitting the third element restrains the starting distance. `--bias-k` sets a global harmonic strength (eV/A^2). Indices default to 1-based but can be flipped to 0-based with `--zero-based`.
6. **Dumping & conversion** -- `--dump` writes `optimization_trj.xyz`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs (with B-factor annotations). `opt.dump_restart` can emit restart YAML snapshots.
7. **Exit codes** -- `0` success, `2` zero step (step norm < `min_step_norm`), `3` optimizer failure, `130` keyboard interrupt, `1` unexpected error.

## CLI options

> **Note:** Default values shown are used when the option is not specified.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader`. | Required |
| `--ref-pdb PATH` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Alternative to `--model-pdb`. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Index convention for `--model-indices`. | 1-based |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors (0/10/20). | Enabled |
| `-q, --charge INT` | Charge of the ML region. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze. | _None_ |
| `--radius-partial-hessian FLOAT` | Distance cutoff (A) from ML region for Hessian-MM atoms. Can be combined with `--detect-layer`. | _None_ |
| `--radius-freeze FLOAT` | Distance cutoff (A) from ML region for movable MM atoms. Atoms beyond this are frozen. | _None_ |
| `--dist-freeze TEXT` | Python-literal `(i, j, target_A)` tuples for harmonic restraints. | _None_ |
| `--one-based / --zero-based` | Index convention for `--dist-freeze`. | 1-based |
| `--bias-k FLOAT` | Harmonic bias strength (eV/A^2). | `10.0` |
| `--max-cycles INT` | Hard limit on optimization iterations. | `10000` |
| `--opt-mode [grad\|hess\|light\|heavy\|lbfgs\|rfo]` | Optimizer mode: `grad` (LBFGS) or `hess` (RFO). Aliases `light`/`heavy` and `lbfgs`/`rfo` accepted. | `grad` |
| `--microiter/--no-microiter` | Microiteration: alternate ML 1-step (RFO) + MM relaxation (LBFGS). Only effective in `hess` mode. | `True` |
| `--flatten/--no-flatten` | Enable/disable the post-optimization imaginary-mode flatten loop. | `False` |
| `--dump/--no-dump` | Emit trajectory dumps (`optimization_trj.xyz`). | `False` |
| `--convert-files/--no-convert-files` | Enable or disable XYZ/TRJ to PDB companions for PDB inputs. | `True` |
| `--out-dir TEXT` | Output directory for all files. | `./result_opt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | `gau` |
| `--config FILE` | Base YAML configuration file. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layer information before execution. | `False` |
| `--backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running optimization. | `False` |

### Convergence threshold presets

Forces in Hartree/bohr, steps in bohr.

| Preset | Purpose | max\|F\| | RMS(F) | max\|step\| | RMS(step) |
| --- | --- | --- | --- | --- | --- |
| `gau_loose` | Loose/quick pre-optimization; rough path searches | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | Standard Gaussian-like tightness for routine work | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | Tighter; better structures / freq / TS refinement | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | Very tight; benchmarking/high-precision final structures | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | Baker-style rule (converged if `max\|F\| < 3e-4` **and** `\|dE\| < 1e-6 or max\|step\| < 3e-4`) | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

## Outputs

```
out_dir/ (default: ./result_opt/)
├─ final_geometry.xyz          # Always written
├─ final_geometry.pdb          # Only when the input was a PDB and conversion is enabled (B-factors annotated)
├─ optimization_trj.xyz        # Only if dumping is enabled
├─ optimization.pdb            # PDB conversion of the trajectory (PDB inputs, conversion enabled)
└─ restart*.yml                # Optional restarts when opt.dump_restart is set
```

Console output prints resolved configuration blocks (`geom`, `calc`, `opt`, `lbfgs`), progress every `print_every` cycles, and a final wall-clock time summary.

## YAML configuration

Settings are applied with **defaults < config < explicit CLI < override**. Accepted sections:

### `geom`

- `coord_type` (`"cart"` default): Cartesian vs. `"dlc"` delocalized internal coordinates.
- `freeze_atoms` (`[]`): 0-based indices to freeze during optimization.

### `calc` / `mlmm`

- `input_pdb`, `real_parm7`, `model_pdb`: required file paths (strings).
- `model_charge` (`-q/--charge`, required) and `model_mult` (`-m/--multiplicity`, default 1).
- `link_mlmm`: optional list of `(ML_atom_id, MM_atom_id)` strings to pin ML/MM link pairs (no link atoms created).
- ML backend controls: `backend` (default `"uma"`; choices `uma`, `orb`, `mace`, `aimnet2`), `embedcharge` (default `false`). UMA-specific: `uma_model` (default `"uma-s-1p1"`), `uma_task_name` (default `"omol"`). Shared: `ml_hessian_mode` (`"Analytical"` or `"FiniteDifference"`), `out_hess_torch` (bool), `H_double` (bool).
- Device selection: `ml_device` (`"auto"`/`"cuda"`/`"cpu"`), `ml_cuda_idx`, `mm_device`, `mm_cuda_idx`, `mm_threads`.
- MM finite difference: `mm_fd` (bool), `mm_fd_dir` (output dir for FD info), and whether to `return_partial_hessian`.
- `return_partial_hessian`: for `opt`, partial Hessian is used by default when this key is not explicitly set in YAML. Set `calc.return_partial_hessian: false` to force full Hessian output.
- `freeze_atoms`: propagated from `geom.freeze_atoms` so ML/MM and optimizer share the same frozen atoms.

### `opt`

Shared optimizer controls:
- `thresh` presets (see convergence table above).
- Common controls: `max_cycles` (default 10000), `print_every` (100), `min_step_norm` (1e-8) with `assert_min_step` True.
- Convergence toggles: `rms_force`, `rms_force_only`, `max_force_only`, `force_only`.
- Extras: `converge_to_geom_rms_thresh`, `overachieve_factor`, `check_eigval_structure`.
- Line search: `line_search` (bool, default True).
- Bookkeeping: `dump`, `dump_restart`, `prefix`, `out_dir` (default `./result_opt/`).

### `lbfgs`

Extends `opt` with L-BFGS specifics: `keep_last`, `beta`, `gamma_mult`, `max_step`, `control_step`, `double_damp`, `mu_reg`, `max_mu_reg_adaptions`.

### `rfo`

Extends `opt` with RFOptimizer fields: trust-region sizing (`trust_radius`, `trust_min`, `trust_max`, `trust_update`), `max_energy_incr`, Hessian management (`hessian_update`, `hessian_init`, `hessian_recalc`, `hessian_recalc_adapt`, `small_eigval_thresh`), micro-iteration controls (`alpha0`, `max_micro_cycles`, `rfo_overlaps`), DIIS helpers (`gdiis`, `gediis`, thresholds, `gdiis_test_direction`), and `adapt_step_func`.

### Example YAML
```yaml
geom:
 coord_type: cart               # coordinate type: cartesian vs dlc internals
 freeze_atoms: []               # 0-based frozen atoms merged with CLI/link detection
calc:
 charge: 0                      # total charge (CLI override)
 spin: 1                        # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7         # Amber parm7 topology for the full enzyme
 model_pdb: ml_region.pdb       # PDB defining the ML region
 backend: uma                   # MLIP backend: uma | orb | mace | aimnet2
 embedcharge: false             # xTB point-charge embedding correction
 uma_model: uma-s-1p1           # UMA model tag (UMA backend only)
 uma_task_name: omol             # UMA task name (UMA backend only)
 ml_device: auto                # ML backend device selection
 ml_hessian_mode: FiniteDifference  # Hessian mode selection
 out_hess_torch: true           # request torch-form Hessian
 mm_fd: false                   # MM finite-difference toggle
 return_partial_hessian: true   # allow partial Hessians (default for opt)
opt:
 thresh: gau                    # convergence preset (Gaussian/Baker-style)
 max_cycles: 10000              # optimizer cycle cap
 print_every: 100               # logging stride
 min_step_norm: 1.0e-08         # minimum norm for step acceptance
 assert_min_step: true          # stop if steps fall below threshold
 rms_force: null                # explicit RMS force target
 rms_force_only: false          # rely only on RMS force convergence
 max_force_only: false          # rely only on max force convergence
 force_only: false              # skip displacement checks
 converge_to_geom_rms_thresh: 0.05  # geom RMS threshold when converging to ref
 overachieve_factor: 0.0        # factor to tighten thresholds
 check_eigval_structure: false  # validate Hessian eigenstructure
 line_search: true              # enable line search
 dump: false                    # dump trajectory/restart data
 dump_restart: false            # dump restart checkpoints
 prefix: ""                     # filename prefix
 out_dir: ./result_opt/         # output directory
lbfgs:
 thresh: gau                    # LBFGS convergence preset
 max_cycles: 10000              # iteration limit
 print_every: 100               # logging stride
 min_step_norm: 1.0e-08         # minimum accepted step norm
 assert_min_step: true          # assert when steps stagnate
 rms_force: null                # explicit RMS force target
 rms_force_only: false          # rely only on RMS force convergence
 max_force_only: false          # rely only on max force convergence
 force_only: false              # skip displacement checks
 converge_to_geom_rms_thresh: 0.05  # RMS threshold when targeting geometry
 overachieve_factor: 0.0        # tighten thresholds
 check_eigval_structure: false  # validate Hessian eigenstructure
 line_search: true              # enable line search
 dump: false                    # dump trajectory/restart data
 dump_restart: false            # dump restart checkpoints
 prefix: ""                     # filename prefix
 out_dir: ./result_opt/         # output directory
 keep_last: 7                   # history size for LBFGS buffers
 beta: 1.0                      # initial damping beta
 gamma_mult: false              # multiplicative gamma update toggle
 max_step: 0.3                  # maximum step length
 control_step: true             # control step length adaptively
 double_damp: true              # double damping safeguard
 mu_reg: null                   # regularization strength
 max_mu_reg_adaptions: 10       # cap on mu adaptations
rfo:
 thresh: gau                    # RFOptimizer convergence preset
 max_cycles: 10000              # iteration cap
 print_every: 100               # logging stride (matches shared opt defaults)
 min_step_norm: 1.0e-08         # minimum accepted step norm
 assert_min_step: true          # assert when steps stagnate
 rms_force: null                # explicit RMS force target
 rms_force_only: false          # rely only on RMS force convergence
 max_force_only: false          # rely only on max force convergence
 force_only: false              # skip displacement checks
 converge_to_geom_rms_thresh: 0.05  # RMS threshold when targeting geometry
 overachieve_factor: 0.0        # tighten thresholds
 check_eigval_structure: false  # validate Hessian eigenstructure
 line_search: true              # enable line search
 dump: false                    # dump trajectory/restart data
 dump_restart: false            # dump restart checkpoints
 prefix: ""                     # filename prefix
 out_dir: ./result_opt/         # output directory
 trust_radius: 0.1              # trust-region radius
 trust_update: true             # enable trust-region updates
 trust_min: 0.0                 # minimum trust radius
 trust_max: 0.1                 # maximum trust radius
 max_energy_incr: null          # allowed energy increase per step
 hessian_update: bfgs           # Hessian update scheme
 hessian_init: calc             # Hessian initialization source
 hessian_recalc: 200            # rebuild Hessian every N steps
 hessian_recalc_adapt: null     # adaptive Hessian rebuild limit
 small_eigval_thresh: 1.0e-08   # eigenvalue threshold for stability
 alpha0: 1.0                    # initial micro step
 max_micro_cycles: 50           # micro-iteration limit
 rfo_overlaps: false            # enable RFO overlaps
 gediis: false                  # enable GEDIIS
 gdiis: true                    # enable GDIIS
 gdiis_thresh: 0.0025           # GDIIS acceptance threshold
 gediis_thresh: 0.01            # GEDIIS acceptance threshold
 gdiis_test_direction: true     # test descent direction before DIIS
 adapt_step_func: true          # adaptive step scaling
```

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli_conventions.md).
- **Devices:** `ml_device="auto"` selects CUDA when available; `mm_device` controls MM backend device placement.
- **Hessians:** `calc.out_hess_torch=True` returns PyTorch tensors (optionally double precision via `calc.H_double`).
- **Freeze atoms:** CLI freeze-link logic is merged with YAML `geom.freeze_atoms`, then propagated to the ML/MM calculator (`calc.freeze_atoms`).
- **Exit codes:** `ZeroStepLength` -> exit code **2**; `OptimizationError` -> **3**; `KeyboardInterrupt` -> **130**; any other unhandled exception -> **1**.
- **Precedence:** settings are applied with **defaults < config < explicit CLI < override**.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [tsopt](tsopt.md) -- Optimize transition states (saddle points) instead of minima
- [freq](freq.md) -- Vibrational analysis to confirm a minimum was reached
- [all](all.md) -- End-to-end workflow that pre-optimizes endpoints
- [YAML Reference](yaml_reference.md) -- Full `opt`, `lbfgs`, `rfo` configuration options
- [Glossary](glossary.md) -- Definitions of L-BFGS, RFO
