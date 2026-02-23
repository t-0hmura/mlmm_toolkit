# `opt`

## Overview

> **Summary:** Optimizes a single structure to a local minimum using L-BFGS (default) or RFO with the ML/MM calculator. Input must be a PDB file; output includes XYZ and PDB geometries with B-factor annotations.


## Minimal example

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_opt
```

## Output checklist

- `result_opt/summary.md`
- `result_opt/key_opt.xyz` (or `key_opt.pdb`)
- `result_opt/key_opt_trj.xyz` (when trajectory is available)
- `result_opt/final_geometry.xyz`
- `result_opt/final_geometry.pdb` (when the input is PDB)
- `result_opt/optimization_trj.xyz` (when `--dump` is enabled)

## Common examples

1. Tighten convergence and keep an optimization trajectory.

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --thresh gau_tight --dump --out-dir ./result_opt_tight
```

2. Apply one harmonic distance restraint during optimization.

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --dist-freeze "[(12,45,2.20)]" --bias-k 20.0 --out-dir ./result_opt_rest
```

3. Switch to heavy mode (RFO).

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 -q 0 --opt-mode heavy --out-dir ./result_opt_rfo
```

## Usage

```bash
mlmm opt -i INPUT.pdb --real-parm7 real.parm7 --model-pdb model.pdb -q CHARGE [-m MULT]
 [--dist-freeze "[(I,J,TARGET_A),...]"] [--one-based|--zero-based] [--bias-k FLOAT]
 [--freeze-atoms "1,3,5"] [--max-cycles N] [--thresh PRESET]
```

### Examples

```bash
mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0

mlmm opt -i pocket.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
 --freeze-atoms "1,3,5,7" --thresh gau_tight --dump --out-dir ./result_opt/ \
```

## Workflow

- **Input handling:** the tool requires `-i/--input` to be a PDB file (the enzyme complex). The optimizer reads coordinates from this PDB via `pysisyphus.helpers.geom_loader`.
- **Harmonic distance restraints** are available via `--dist-freeze` with strength set by `--bias-k` (eV/Angstrom^2).
- **PDB conversion note:** after converting XYZ/TRJ to PDB, B-factors are annotated as: ML-region atoms = 100.00, frozen atoms = 50.00, atoms that are both ML and frozen = 150.00.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure file (PDB or XYZ). | Required |
| `--ref-pdb PATH` | Reference PDB topology when input is XYZ. | _None_ |
| `--real-parm7 PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Alternative to `--model-pdb`. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Index convention for `--model-indices`. | 1-based |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors (0/10/20). | Enabled |
| `-q, --charge INT` | Charge of the ML region. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze. | _None_ |
| `--radius-partial-hessian FLOAT` | Distance cutoff (A) from ML region for Hessian-MM atoms. Can be combined with `--detect-layer`. | _None_ |
| `--radius-freeze FLOAT` | Distance cutoff (A) from ML region for movable MM atoms. Atoms beyond this are frozen. | _None_ |
| `--dist-freeze TEXT` | Python-literal `(i, j, target_A)` tuples for harmonic restraints. | _None_ |
| `--one-based / --zero-based` | Index convention for `--dist-freeze`. | 1-based |
| `--bias-k FLOAT` | Harmonic bias strength (eV/A^2). | `10.0` |
| `--max-cycles INT` | Hard limit on optimization iterations. | `10000` |
| `--thresh TEXT` | Convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`). | _None_ |
| `--opt-mode {light\|heavy\|lbfgs\|rfo}` | Optimizer mode: `light` (LBFGS) or `heavy` (RFO with Hessian). | `light` |
| `--layer-opt / --no-layer-opt` | Enable LayerOpt-style microiteration (heavy mode only). | Disabled |
| `--dump/--no-dump` | Emit trajectory dumps (`optimization_trj.xyz`). | `False` |
| `--out-dir PATH` | Output directory. | `./result_opt/` |
| `--config FILE` | Base YAML configuration file. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layer information before execution. | `False` |
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

```text
out_dir/ (default:./result_opt/)
 summary.md # Quick index of key outputs
 key_opt.xyz # Shortcut to final_geometry.xyz
 key_opt.pdb # Shortcut to final_geometry.pdb (when available)
 key_opt_trj.xyz # Shortcut to optimization trajectory
 key_opt_traj.pdb # Shortcut to optimization trajectory PDB (when available)
 key_restart.yml # Shortcut to a restart snapshot (when available)
 final_geometry.xyz # Final optimized geometry (always)
 final_geometry.pdb # Converted from XYZ when the input was a PDB (B-factors annotated)
 optimization_trj.xyz # Trajectory (written when --dump or opt.dump: true)
 optimization.pdb # Converted from TRJ when input was PDB and dumping is enabled
 restart*.yml # Restart files every opt.dump_restart cycles (if enabled)
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
- UMA controls: `uma_model` (default `"uma-s-1p1"`), `uma_task_name` (default `"omol"`), `ml_hessian_mode` (`"Analytical"` or `"FiniteDifference"`), `out_hess_torch` (bool), `H_double` (bool).
- Device selection: `ml_device` (`"auto"`/`"cuda"`/`"cpu"`), `ml_cuda_idx`, `mm_device`, `mm_cuda_idx`, `mm_threads`.
- MM finite difference: `mm_fd` (bool), `mm_fd_dir` (output dir for FD info), and whether to `return_partial_hessian`.
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

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli_conventions.md).
- **Devices:** `ml_device="auto"` selects CUDA when available; `mm_device` controls MM backend device placement.
- **Hessians:** `calc.out_hess_torch=True` returns PyTorch tensors (optionally double precision via `calc.H_double`).
- **Exit codes:** `ZeroStepLength` -> exit code **2**; `OptimizationError` -> **3**; `KeyboardInterrupt` -> **130**; any other unhandled exception -> **1**.
- **Precedence:** settings are applied with **defaults < config < explicit CLI < override**.

---

## See Also

- [tsopt](tsopt.md) -- Optimize transition states (saddle points) instead of minima
- [freq](freq.md) -- Vibrational analysis to confirm a minimum was reached
- [all](all.md) -- End-to-end workflow that pre-optimizes endpoints
