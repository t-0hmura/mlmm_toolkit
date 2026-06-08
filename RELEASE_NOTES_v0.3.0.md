# mlmm-toolkit v0.3.0

`mlmm-toolkit` is an open-source CLI for **ML/MM ONIOM** analyses of enzymatic
reactions. It replaces the QM region of conventional QM/MM with a
machine-learning interatomic potential (MLIP, default: UMA) while keeping the
surrounding protein under an analytical Amber force field (`hessian_ff`), and
chains **MM parametrisation → ML-region selection → MEP search → TS optimisation
→ IRC → frequencies → DFT single-point** in one command. A link-atom boundary
handles amino-acid residues straddling the ML/MM cut, and a microiteration
scheme makes TS optimisation and Hessian-based methods tractable on
~10 000-atom systems.

A useful initial reaction path is one command:

```bash
# Multi-structure MEP (R + P endpoints → MEP, with TS optimisation + thermo)
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' --tsopt --thermo
```

## What's new in v0.3.0

Major release: the package is refactored into a clean layered architecture
(cli / workflows / domain / backends / io / core). The CLI surface is unchanged,
so existing commands keep working.

### Breaking
- Removed the `mlmm pysis` subcommand and the `mlmm.pysis_runner` module.

### Added
- `--deterministic` on every compute subcommand for bit-reproducible GPU runs.
- `--precision fp32|fp64` on every calculator-constructing subcommand
  (fp64 also forces the Hessian to fp64).
- `--irc-pos-def` IRC convergence guard.
- TS-opt microiteration now supports internal-coordinate macro steps.
- Structured `result.json` / `summary.json` envelope (`schema_version`,
  error-class chain), mirrored to a single `summary.json` per subcommand.
- New docs: `reproducibility.md`, `output-layout.md`.

### Changed
- `--help` groups subcommands into semantic sections.
- Standalone `path-opt` / `path-search` default to `--preopt`.
- AIMNet2 rejects `--precision fp64` and `--deterministic` with clear errors.
- Cleaner CLI error rendering with a recovery hint.

See `CHANGELOG.md` for the complete list.
