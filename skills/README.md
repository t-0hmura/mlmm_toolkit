# Agent Skills for `mlmm-toolkit`

This folder contains skills that common AI agent interfaces recognize.
They give an agent concise instructions for using the `mlmm-toolkit` CLI:
which command to run and when, how to judge success, common pitfalls and
how to recover, and how to read the outputs. Inspired by the
[nvalchemi-toolkit][nvalchemi] skill pattern.

[nvalchemi]: https://github.com/NVIDIA/nvalchemi-toolkit

- `mlmm-overview` (start here): what `mlmm-toolkit` is and which of the three `all` modes
  fits your structures; `ts-strategy.md` for TS strategy (imaginary-mode count,
  flattening, a TS that does not come out, mutant-vs-WT comparisons);
  `outputs.md` for reading the outputs.
- `mlmm-cli`: the 22 subcommands in 17 files, each with when to use it, how to
  judge success, and pitfalls and recovery.
- `mlmm-mcp`: how to drive `mlmm-toolkit` from any MCP client (Claude
  Desktop / Claude Code / Cursor / custom SDK) via the bundled
  `mlmm-mcp` server; lists the 22 MCP tools (including the mlmm-specific
  topology / ONIOM-layer / ONIOM-input tools) and the result format
  shared by every tool.
- `mlmm-structure-io`: PDB / mmCIF / XYZ / GJF / Amber parm7/rst7
  handling and the charge / multiplicity decision workflow.
- `mlmm-model-setup`: what `extract` puts in the ML region, a hand-built
  `model.pdb` with `--model-pdb`, `--parm7`, and `-q`, the `define-layer` layers,
  trimming or enlarging the model, and the same-atom rules for R/IM/P and
  WT/mutant models; the full guide is [`docs/model-setup.md`](../docs/model-setup.md).
- `mlmm-install-backends`: install mlmm itself, MLIP backends (UMA / Orb / MACE /
  AIMNet2), DFT (PySCF / GPU4PySCF), xtb (the `--embedcharge` point-charge correction,
  not an MLIP backend), and AmberTools (tleap); CUDA + PyTorch pairing; probing an
  unknown scheduler / GPU / CUDA / conda env.
- `mlmm-hpc`: PBS / SLURM preamble templates with placeholders,
  walltime guidance, monitoring, plus a flock+pbsdsh dynamic-dispatch
  recipe.
- `colab-local-gpu-runtime`: Windows setup and operation for running the Colab
  interface on a local NVIDIA GPU through WSL2 and Docker Desktop.

## Install

Each folder here is one skill (`<name>/SKILL.md`). Copy the folders into your agent's skill directory:

- Claude Code: `.claude/skills/` in a project, or `~/.claude/skills/` for all projects.
- Codex: `.agents/skills/` in a repository, or `~/.agents/skills/` for all repositories.

For example, from the root of this repository:

```bash
mkdir -p ~/.claude/skills
cp -r skills/mlmm-* skills/colab-local-gpu-runtime ~/.claude/skills/
```

For exact flags and defaults, check the installed CLI (`mlmm <subcommand> --help-advanced`) and the [command reference](../docs/reference/commands/index.md).
