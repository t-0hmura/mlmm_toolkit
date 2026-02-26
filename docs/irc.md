# `irc`

## Overview

> **Summary:** Runs EulerPC-based IRC (Intrinsic Reaction Coordinate) integration from a transition state toward reactants and products using the ML/MM calculator. By default both forward and backward branches are computed.

### At a glance
- **Use when:** You have an optimized TS and want to trace the minimum-energy path toward reactant and product basins with ML/MM.
- **Method:** EulerPC predictor-corrector integrator with full ML/MM Hessians (FAIR-Chem UMA + hessian_ff).
- **Outputs:** `finished_irc_trj.xyz`, `forward_irc_trj.xyz`, and `.pdb` companions for PDB inputs.
- **Next step:** Run [freq](freq.md) on IRC endpoints, then [opt](opt.md) to refine them to true minima.

`mlmm irc` runs IRC calculations using the EulerPC integrator with the ML/MM calculator. The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible. Inputs can be any structure readable by `pysisyphus.helpers.geom_loader` (`.pdb`, `.xyz`, `_trj.xyz`,...). If the input is `.pdb`, the generated trajectories are additionally converted to PDB.

A typical workflow is `tsopt` -> `freq` (confirm **one** imaginary mode) -> `irc`.

## Minimal example

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## Output checklist

- `result_irc/finished_irc_trj.xyz`
- `result_irc/forward_irc_trj.xyz`

## Common examples

1. Run only the forward branch.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-forward --out-dir ./result_irc_forward
```

2. Increase step size and use analytical Hessians.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.20 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
```

3. Keep both branches and raise the step limit.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 150 \
 --out-dir ./result_irc_long
```

## Workflow

1. **Input preparation** -- Any format supported by `geom_loader` is accepted. When a reference PDB is available (input is `.pdb` or `--ref-pdb` is supplied), EulerPC trajectories are converted to PDB using that topology.
2. **ML/MM calculator setup** -- Build the ML/MM calculator from `--real-parm7` and `--model-pdb`. The `--hessian-calc-mode` controls UMA Hessian evaluation.
4. **IRC integration** -- The EulerPC integrator propagates along the IRC in both directions (unless `--no-forward` disables forward). Step size and cycle count control integration length.
5. **Output & conversion** -- Trajectories are written as XYZ; PDB companions are generated when a PDB template is available and `--convert-files` is enabled.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file (`.pdb`/`.xyz`/`_trj.xyz`/...). | Required |
| `--real-parm7 PATH` | Amber topology for the full enzyme/MM region. Required unless `calc.real_parm7` is set in YAML. | _None_ |
| `--model-pdb PATH` | PDB defining the ML region. Required when `--no-detect-layer` and no `--model-indices` are given. | _None_ |
| `--model-indices TEXT` | Comma-separated ML-region atom indices (ranges allowed, e.g. `1-10,15`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer/--no-detect-layer` | Detect ML/MM layers from input PDB B-factors (`B=0/10/20`). | `True` |
| `-q, --charge INT` | Total charge; overrides `calc.charge` from YAML. | Required |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `1` |
| `--max-cycles INT` | Max number of IRC steps; overrides `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; overrides `irc.step_length`. | `0.10` |
| `--root INT` | Imaginary mode index for the initial displacement; overrides `irc.root`. | `0` |
| `--forward/--no-forward` | Run the forward IRC; overrides `irc.forward`. | `True` |
| `--out-dir PATH` | Output directory; overrides `irc.out_dir`. | `./result_irc/` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a reference PDB is available. | `True` |
| `--hessian-calc-mode CHOICE` | How UMA builds the Hessian (`Analytical` or `FiniteDifference`); overrides `calc.hessian_calc_mode`. | _Default_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## Outputs

```
out_dir/ (default: ./result_irc/)
├─ <prefix>irc_data.h5              # HDF5 dump written every irc.dump_every steps
├─ <prefix>finished_irc_trj.xyz     # Full IRC trajectory (XYZ/TRJ)
├─ <prefix>forward_irc_trj.xyz      # Forward path segment
├─ <prefix>finished_irc.pdb         # PDB conversion (only if input was .pdb)
└─ <prefix>forward_irc.pdb          # PDB conversion (only if input was .pdb)
```

## YAML configuration

Provide mappings with merge order **defaults < config < explicit CLI < override**.
Shared sections reuse [YAML Reference](yaml_reference.md) for geometry/calculator keys. For `irc`, `geom.coord_type` is forced to `cart` and `calc.return_partial_hessian` is forced to `false` after YAML/CLI merging.

### CLI-to-YAML mapping

| CLI option | YAML key |
|------------|----------|
| `--charge` | `calc.charge` |
| `--multiplicity` | `calc.spin` |
| `--step-size` | `irc.step_length` |
| `--max-cycles` | `irc.max_cycles` |
| `--root` | `irc.root` |
| `--forward` | `irc.forward` |
| `--out-dir` | `irc.out_dir` |
| `--hessian-calc-mode` | `calc.hessian_calc_mode` |

### Example YAML

```yaml
geom:
 coord_type: cart                  # forced to cart for irc (YAML value ignored)
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
 return_partial_hessian: false     # forced false for irc (full Hessian)
irc:
 step_length: 0.1                  # integration step length
 max_cycles: 125                   # maximum steps along IRC
 downhill: false                   # follow downhill direction only
 forward: true                     # propagate in forward direction
 root: 0                           # normal-mode root index
 hessian_init: calc                # Hessian initialization source
 displ: energy                     # displacement construction method
 displ_energy: 0.001               # energy-based displacement scaling
 displ_length: 0.1                 # length-based displacement fallback
 rms_grad_thresh: 0.001            # RMS gradient convergence threshold
 hard_rms_grad_thresh: null        # hard RMS gradient stop
 energy_thresh: 0.000001           # energy change threshold
 imag_below: 0.0                   # imaginary frequency cutoff
 force_inflection: true            # enforce inflection detection
 check_bonds: false                # check bonds during propagation
 out_dir: ./result_irc/            # output directory
 prefix: ""                        # filename prefix
 hessian_update: bofill            # Hessian update scheme
 hessian_recalc: null              # Hessian rebuild cadence
 max_pred_steps: 500               # predictor-corrector max steps
 loose_cycles: 3                   # loose cycles before tightening
 corr_func: mbs                    # correlation function choice
```

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli_conventions.md).
- UMA options are passed directly to the mlmm calculator. With `device: "auto"`, the calculator selects GPU/CPU automatically.
- When you have ample VRAM available, setting `--hessian-calc-mode` to `Analytical` is strongly recommended.
- `irc` is partial-first: when YAML does not explicitly set `calc.return_partial_hessian`, IRC seeds/uses a partial Hessian by default. Set `calc.return_partial_hessian: false` to force full Hessian.
- If `hessian_calc_mode: "FiniteDifference"`, `geom.freeze_atoms` can still be used to skip frozen DOF in FD Hessian construction.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used for the initial displacement.
- Standard output includes progress and timing. Exit codes: `0` on success, `130` on `KeyboardInterrupt`, `1` on unhandled exceptions.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [tsopt](tsopt.md) -- Optimize the TS before running IRC
- [freq](freq.md) -- Verify the TS candidate has one imaginary frequency; analyze IRC endpoints
- [opt](opt.md) -- Optimize IRC endpoints to true minima
- [all](all.md) -- End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml_reference.md) -- Full `irc` configuration options
- [Glossary](glossary.md) -- Definition of IRC (Intrinsic Reaction Coordinate)
