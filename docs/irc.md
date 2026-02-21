# `irc`

## Overview

> **Summary:** Trace the intrinsic reaction coordinate (IRC) from a transition state toward reactant and product using the EulerPC predictor-corrector integrator with the ML/MM calculator.

`mlmm irc` runs IRC calculations using the EulerPC integrator. The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible. Inputs can be any structure readable by `pysisyphus.helpers.geom_loader` (`.pdb`, `.xyz`, `.trj`, ...). If the input is `.pdb`, the generated trajectories are additionally converted to PDB.

Configuration precedence is: **defaults < `--config` < explicit CLI < `--override-yaml`**. Legacy `--args-yaml` remains as an alias of `--override-yaml`.

## Minimal example

```bash
mlmm irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

## Output checklist

- `result_irc/summary.md`
- `result_irc/key_irc.trj`
- `result_irc/key_irc_forward.trj`
- `result_irc/finished_irc.trj`
- `result_irc/forward_irc.trj`
- `result_irc/backward_irc.trj`

## Common examples

1. Run only the forward branch.

```bash
mlmm irc -i ts.xyz -q -1 -m 2 --forward --no-backward \
  --out-dir ./result_irc_forward
```

2. Increase step size and change the displacement root.

```bash
mlmm irc -i ts.xyz -q 0 -m 1 --step-size 0.20 --root 1 \
  --out-dir ./result_irc_step
```

3. Force finite-difference Hessians for stability checks.

```bash
mlmm irc -i ts.pdb -q 0 -m 1 --hessian-calc-mode FiniteDifference \
  --max-cycles 80 --out-dir ./result_irc_fd
```

## Usage

```bash
mlmm irc -i INPUT.pdb -q CHARGE [-m MULT]
    [--max-cycles N] [--step-size Ds] [--root k]
    [--forward/--no-forward] [--backward/--no-backward]
    [--out-dir DIR]
    [--hessian-calc-mode Analytical|FiniteDifference]
    [--config FILE] [--override-yaml FILE|--args-yaml FILE]
    [--show-config] [--dry-run]
```

### Examples

```bash
# Forward-only with finite-difference Hessian and custom step size
mlmm irc -i ts.xyz -q -1 -m 2 --forward --no-backward \
    --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# Use a PDB input so trajectories are also exported as PDB
mlmm irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## Workflow

1. **Input preparation** -- Any format supported by `geom_loader` is accepted. If the input is `.pdb`, trajectories are also converted to PDB format.
2. **Configuration merge** -- Defaults -> `--config` -> explicit CLI -> `--override-yaml`, merging sections `geom`, `calc`, and `irc`. `--args-yaml` is a legacy alias of `--override-yaml`. It is strongly recommended to provide both `-q/--charge` and `-m/--multiplicity` explicitly.
3. **IRC integration** -- EulerPC integrates forward/backward branches according to `irc.forward/backward`, `irc.step_length`, `irc.root`, and the Hessian workflow configured through UMA.
4. **Outputs** -- Trajectories (`finished`, `forward`, `backward`) are written as `.trj` and, when the input was `.pdb`, additionally as `.pdb`.

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file (`.pdb`/`.xyz`/`.trj`/...). | Required |
| `-q, --charge INT` | Total charge; overrides `calc.charge` from YAML. | Strongly recommended |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `1` |
| `--max-cycles INT` | Max number of IRC steps; overrides `irc.max_cycles`. | _Default_ |
| `--step-size FLOAT` | Step length in mass-weighted coordinates; overrides `irc.step_length`. | _Default_ |
| `--root INT` | Imaginary mode index for the initial displacement; overrides `irc.root`. | _Default_ |
| `--forward/--no-forward` | Run the forward IRC; overrides `irc.forward`. | `True` |
| `--backward/--no-backward` | Run the backward IRC; overrides `irc.backward`. | `True` |
| `--out-dir PATH` | Output directory; overrides `irc.out_dir`. | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | How UMA builds the Hessian (`Analytical` or `FiniteDifference`); overrides `calc.hessian_calc_mode`. | _Default_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--override-yaml FILE` | Final YAML override (highest-priority YAML layer). | _None_ |
| `--args-yaml FILE` | Legacy alias of `--override-yaml`. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. | `False` |

## Outputs

```text
out_dir/ (default: ./result_irc/)
  summary.md                      # Quick index of key outputs
  key_irc.trj                     # Shortcut to finished_irc.trj
  key_irc_forward.trj             # Shortcut to forward_irc.trj
  key_irc_backward.trj            # Shortcut to backward_irc.trj
  key_irc.pdb                     # Shortcut to finished_irc.pdb (when available)
  key_irc_data.h5                 # Shortcut to irc_data.h5 (when available)
  <prefix>irc_data.h5              # HDF5 dump written every irc.dump_every steps
  <prefix>finished_irc.trj         # Full IRC trajectory (XYZ/TRJ)
  <prefix>forward_irc.trj          # Forward path segment
  <prefix>backward_irc.trj         # Backward path segment
  <prefix>finished_irc.pdb         # PDB conversion (only if input was .pdb)
  <prefix>forward_irc.pdb          # PDB conversion (only if input was .pdb)
  <prefix>backward_irc.pdb         # PDB conversion (only if input was .pdb)
```

## YAML configuration

Provide a YAML mapping with sections `geom`, `calc`, and `irc`. Merge order is defaults < config < explicit CLI < override (`--args-yaml` is the legacy alias of `--override-yaml`).

### CLI-to-YAML mapping

| CLI option | YAML key |
|------------|----------|
| `--charge` | `calc.charge` |
| `--multiplicity` | `calc.spin` |
| `--step-size` | `irc.step_length` |
| `--max-cycles` | `irc.max_cycles` |
| `--root` | `irc.root` |
| `--forward` | `irc.forward` |
| `--backward` | `irc.backward` |
| `--out-dir` | `irc.out_dir` |
| `--hessian-calc-mode` | `calc.hessian_calc_mode` |

## Notes

- For symptom-first diagnosis, start with [Common Error Recipes](recipes-common-errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Charge/multiplicity policy is documented centrally in [CLI Conventions](cli-conventions.md).
- UMA options are passed directly to the mlmm calculator. With `device: "auto"`, the calculator selects GPU/CPU automatically.
- If `hessian_calc_mode: "FiniteDifference"`, `geom.freeze_atoms` can be used to skip frozen DOF in FD Hessian construction.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used for the initial displacement.
- Standard output includes progress and timing. Exit codes: `0` on success, `130` on `KeyboardInterrupt`, `1` on unhandled exceptions.

---

## See Also

- [tsopt](tsopt.md) -- Optimize the TS before running IRC
- [freq](freq.md) -- Verify the TS candidate has one imaginary frequency
- [all](all.md) -- End-to-end workflow that runs IRC after tsopt
