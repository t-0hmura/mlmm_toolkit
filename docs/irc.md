# `irc`

## Overview

> **Summary:** Trace the intrinsic reaction coordinate (IRC) from a transition state toward reactant and product using the EulerPC predictor-corrector integrator with the ML/MM calculator.

`mlmm irc` runs IRC calculations using the EulerPC integrator. The CLI is intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible. Inputs can be any structure readable by `pysisyphus.helpers.geom_loader` (`.pdb`, `.xyz`, `.trj`, ...). If the input is `.pdb`, the generated trajectories are additionally converted to PDB.

Configuration precedence is: **built-in defaults -> YAML -> CLI**.

## Usage

```bash
mlmm irc -i INPUT.pdb -q CHARGE [-m MULT]
    [--max-cycles N] [--step-size Ds] [--root k]
    [--forward True|False] [--backward True|False]
    [--out-dir DIR]
    [--hessian-calc-mode Analytical|FiniteDifference]
    [--args-yaml FILE]
```

### Examples

```bash
# Forward-only with finite-difference Hessian and custom step size
mlmm irc -i ts.xyz -q -1 -s 2 --forward True --backward False \
    --step-size 0.2 --hessian-calc-mode FiniteDifference --out-dir ./irc_fd/

# Use a PDB input so trajectories are also exported as PDB
mlmm irc -i ts.pdb -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc/
```

## Workflow

1. **Input preparation** -- Any format supported by `geom_loader` is accepted. If the input is `.pdb`, trajectories are also converted to PDB format.
2. **Configuration merge** -- Defaults -> CLI -> YAML, merging sections `geom`, `calc`, and `irc`. It is strongly recommended to provide both `-q/--charge` and `-m/--multiplicity` explicitly.
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
| `--forward {True\|False}` | Run the forward IRC; overrides `irc.forward`. | `True` |
| `--backward {True\|False}` | Run the backward IRC; overrides `irc.backward`. | `True` |
| `--out-dir PATH` | Output directory; overrides `irc.out_dir`. | `./result_irc/` |
| `--hessian-calc-mode CHOICE` | How UMA builds the Hessian (`Analytical` or `FiniteDifference`); overrides `calc.hessian_calc_mode`. | _Default_ |
| `--args-yaml FILE` | YAML file with sections `geom`, `calc`, and `irc`. | _None_ |

## Outputs

```text
out_dir/ (default: ./result_irc/)
  <prefix>irc_data.h5              # HDF5 dump written every irc.dump_every steps
  <prefix>finished_irc.trj         # Full IRC trajectory (XYZ/TRJ)
  <prefix>forward_irc.trj          # Forward path segment
  <prefix>backward_irc.trj         # Backward path segment
  <prefix>finished_irc.pdb         # PDB conversion (only if input was .pdb)
  <prefix>forward_irc.pdb          # PDB conversion (only if input was .pdb)
  <prefix>backward_irc.pdb         # PDB conversion (only if input was .pdb)
```

## YAML configuration

Provide a YAML mapping with sections `geom`, `calc`, and `irc`. YAML values override CLI.

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

- **Strong recommendation:** Provide both `-q/--charge` and `-m/--multiplicity` explicitly to avoid running under unintended conditions.
- UMA options are passed directly to the mlmm calculator. With `device: "auto"`, the calculator selects GPU/CPU automatically.
- If `hessian_calc_mode: "FiniteDifference"`, `geom.freeze_atoms` can be used to skip frozen DOF in FD Hessian construction.
- `--step-size` is in mass-weighted coordinates; `--root` selects the imaginary-frequency index used for the initial displacement.
- Standard output includes progress and timing. Exit codes: `0` on success, `130` on `KeyboardInterrupt`, `1` on unhandled exceptions.

---

## See Also

- [tsopt](tsopt.md) -- Optimize the TS before running IRC
- [freq](freq.md) -- Verify the TS candidate has one imaginary frequency
- [all](all.md) -- End-to-end workflow that runs IRC after tsopt
