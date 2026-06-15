# Agent Skills for `mlmm-toolkit`

This folder contains a set of skills that common AI agent interfaces
will recognize and help speed up code development by providing
concise instructions on how to use the `mlmm-toolkit` API and CLI for
common tasks. Inspired by the [nvalchemi-toolkit][nvalchemi] skill
pattern.

[nvalchemi]: https://github.com/NVIDIA/nvalchemi-toolkit

- `mlmm-overview`: what `mlmm-toolkit` is, when to use it,
  and how it differs from generic QM/MLIP path-search tools.
- `mlmm-architecture`: 6-layer package map (cli / workflows / domain /
  backends / io / core) + 3 bundled forks (pysisyphus / thermoanalysis /
  hessian_ff); tells an agent which dir to grep for a given concern
  before touching code.
- `mlmm-cli`: index of the 22 subcommands plus per-subcommand
  mds (each with synopsis, key flags, examples, output, caveats).
- `mlmm-ts-strategy`: cross-cutting decision know-how for a
  reaction-barrier campaign — precision by GPU class (`--precision`
  fp32/fp64), the two TS-candidate routes (`path-search` MEP vs
  distance-restrained `scan`), fixing a wrong imaginary-frequency count
  (`--precision fp64` / `--coord-type dlc`), reading a barrier when the
  scan started from the Product side, staged vs concerted `--scan-lists`,
  and the same-atom-set rule for controlled mutant-vs-WT comparisons
  (B-factor layer transplant + `--detect-layer`).
- `mlmm-mcp`: how to drive `mlmm-toolkit` from any MCP client (Claude
  Desktop / Claude Code / Cursor / custom SDK) via the bundled
  `mlmm-mcp` server; lists the 22 MCP tools (including the mlmm-specific
  topology / ONIOM-layer / ONIOM-input tools) and the shared
  `SubcmdResult` return schema.
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

The skills are **self-contained**: copying this `skills/` directory
into another project as `.claude/skills/` (or merging it into
`~/.claude/skills/`) gives an agent everything it needs to work with
`mlmm-toolkit` without consulting the main documentation.
