# Agent Skills for `mlmm-toolkit`

This folder contains a set of skills that common AI agent interfaces
will recognize and help speed up code development by providing
concise instructions on how to use the `mlmm-toolkit` API and CLI for
common tasks. Inspired by the [nvalchemi-toolkit][nvalchemi] skill
pattern.

[nvalchemi]: https://github.com/NVIDIA/nvalchemi-toolkit

- `mlmm-overview`: what `mlmm-toolkit` is, when to use it,
  and how it differs from generic QM/MLIP path-search tools.
- `mlmm-cli`: index of the 22 subcommands plus per-subcommand
  mds (each with synopsis, key flags, examples, output, caveats).
- `mlmm-structure-io`: PDB / XYZ / GJF format references and
  the charge / multiplicity decision workflow.
- `mlmm-install-backends`: install mlmm itself, MLIP
  backends (UMA / Orb / MACE / AIMNet2), DFT (PySCF / GPU4PySCF), and
  xtb; CUDA + PyTorch pairing.
- `mlmm-workflows-output`: canonical workflows (cluster /
  multistep / scan-list / endpoint-MEP / TS-only / DFT//MLIP), output
  schema, and R/TS/P canonical paths.
- `mlmm-hpc`: PBS / SLURM preamble templates with placeholders,
  walltime guidance, monitoring, plus a flock+pbsdsh dynamic-dispatch
  recipe.
- `mlmm-env-detect`: fallback for detecting scheduler / GPU /
  CUDA / conda env when the environment is unknown.

The skills are **self-contained**: copying this `.claude/skills/`
directory into another project (or `~/.claude/skills/`) gives an agent
everything it needs to work with `mlmm-toolkit` without consulting the
main documentation.
