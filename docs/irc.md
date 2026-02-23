# `irc`

## Overview

> **Summary:** Trace the intrinsic reaction coordinate (IRC) from a transition state toward reactant and product using the EulerPC predictor-corrector integrator with the ML/MM calculator.

`mlmm irc` runs IRC calculations using the EulerPC integrator. The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible. Inputs can be any structure readable by `pysisyphus.helpers.geom_loader` (`.pdb`, `.xyz`, `.trj`,...). If the input is `.pdb`, the generated trajectories are additionally converted to PDB.


## Minimal example

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

```bash
# YAML-driven minimal setup (real_parm7/model_pdb provided in config)
cat > irc_min.yaml << 'YAML'
calc:
 real_parm7: real.parm7
 model_pdb: ml_region.pdb
 use_bfactor_layers: false
YAML
mlmm irc -i ts.pdb -q 0 -m 1 --config irc_min.yaml \
 --max-cycles 50 --out-dir ./result_irc_yaml
```

## Output checklist

- `result_irc/summary.md`
- `result_irc/key_irc.trj`
- `result_irc/key_irc_forward.trj`
- `result_irc/finished_irc.trj`
- `result_irc/forward_irc.trj`

## Common examples

1. Run only the forward branch.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --out-dir ./result_irc_forward
```

2. Increase step size and change the displacement root.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.20 --root 1 \
 --out-dir ./result_irc_step
```

3. Force finite-difference Hessians for stability checks.

```bash
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --hessian-calc-mode FiniteDifference \
 --max-cycles 80 --out-dir ./result_irc_fd
```

## Usage

```bash
mlmm irc -i INPUT.pdb --real-parm7 real.parm7
 [--model-pdb ml_region.pdb | --model-indices "1,2,3" | --detect-layer]
 [-q CHARGE] [-m MULT]
 [--max-cycles N] [--step-size Ds] [--root k]
 [--detect-layer/--no-detect-layer]
 [--out-dir DIR]
 [--hessian-calc-mode Analytical|FiniteDifference]
 [--show-config] [--dry-run]
```

### Examples

```bash
# Forward-only with finite-difference Hessian and custom step size
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# Use a PDB input so trajectories are also exported as PDB
mlmm irc -i ts.pdb --real-parm7 real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## Workflow

1. **Input preparation** -- Any format supported by `geom_loader` is accepted. If the input is `.pdb`, trajectories are also converted to PDB format.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file (`.pdb`/`.xyz`/`.trj`/...). | Required |
| `--real-parm7 PATH` | Amber topology for the full enzyme/MM region. Required unless `calc.real_parm7` is set in YAML. | _None_ |
| `--model-pdb PATH` | PDB defining the ML region. Required when `--no-detect-layer` and no `--model-indices` are given. | _None_ |
| `--model-indices TEXT` | Comma-separated ML-region atom indices (ranges allowed, e.g. `1-10,15`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer/--no-detect-layer` | Detect ML/MM layers from input PDB B-factors (`B=0/10/20`). | `True` |
| `-q, --charge INT` | Total charge; overrides `calc.charge` from YAML. | Strongly recommended |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `1` |
| `--max-cycles INT` | Max number of IRC steps; overrides `irc.max_cycles`. | _Default_ |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; overrides `irc.step_length`. | _Default_ |
| `--root INT` | Imaginary mode index for the initial displacement; overrides `irc.root`. | _Default_ |
| `--forward/--no-forward` | Run the forward IRC; overrides `irc.forward`. | `True` |
| `--out-dir PATH` | Output directory; overrides `irc.out_dir`. | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | How UMA builds the Hessian (`Analytical` or `FiniteDifference`); overrides `calc.hessian_calc_mode`. | _Default_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## Outputs

```text
out_dir/ (default:./result_irc/)
 summary.md # Quick index of key outputs
 key_irc.trj # Shortcut to finished_irc.trj
 key_irc_forward.trj # Shortcut to forward_irc.trj
 key_irc.pdb # Shortcut to finished_irc.pdb (when available)
 key_irc_data.h5 # Shortcut to irc_data.h5 (when available)
 <prefix>irc_data.h5 # HDF5 dump written every irc.dump_every steps
 <prefix>finished_irc.trj # Full IRC trajectory (XYZ/TRJ)
 <prefix>forward_irc.trj # Forward path segment
 <prefix>finished_irc.pdb # PDB conversion (only if input was.pdb)
 <prefix>forward_irc.pdb # PDB conversion (only if input was.pdb)
```

## YAML configuration


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

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli_conventions.md).
- UMA options are passed directly to the mlmm calculator. With `device: "auto"`, the calculator selects GPU/CPU automatically.
- If `hessian_calc_mode: "FiniteDifference"`, `geom.freeze_atoms` can be used to skip frozen DOF in FD Hessian construction.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used for the initial displacement.
- Standard output includes progress and timing. Exit codes: `0` on success, `130` on `KeyboardInterrupt`, `1` on unhandled exceptions.

---

## See Also

- [tsopt](tsopt.md) -- Optimize the TS before running IRC
- [freq](freq.md) -- Verify the TS candidate has one imaginary frequency
- [all](all.md) -- End-to-end workflow that runs IRC after tsopt
