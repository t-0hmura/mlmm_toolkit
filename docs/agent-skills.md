# Skills for AI Agents

`mlmm-toolkit` ships instructions for AI agents (Claude Code, Codex, Cursor, and others) in the `skills/` directory. With the skills loaded, an agent can run and analyze calculations from plain-language instructions.

## First example

To install them, tell your AI agent:

> Import `https://github.com/t-0hmura/mlmm_toolkit/tree/main/skills` as skills, and install mlmm-toolkit by following `mlmm-install`.

If you cloned the GitHub repository, you can give the local `skills/` path instead. Then you can ask, for example:

> Read *the paper*, build a model from the PDB structure *PDB ID*, and study the mechanism of *the reaction step* with the mlmm-toolkit skills.

## Skills

Each folder in `skills/` is one skill (`<name>/SKILL.md`).

| Skill | What it covers |
| --- | --- |
| `mlmm-overview` (start here) | What mlmm-toolkit is, which of the three `all` modes fits your structures, TS strategy, and how to read the outputs |
| `mlmm-cli` | When to use each subcommand, how to judge success, and pitfalls and recovery |
| `mlmm-model-setup` | Structure files and Amber parm7/rst7, the charge and multiplicity, the ML region and layers, and keeping the same atoms across states and variants |
| `mlmm-install` | Installing mlmm-toolkit, the MLIP backends, DFT, xTB, and AmberTools, and probing an unknown environment |
| `mlmm-hpc` | PBS / SLURM job scripts, walltime, monitoring, and dynamic dispatch |
| `colab-local-gpu-runtime` | Running the Colab interface on a local NVIDIA GPU through WSL2 and Docker Desktop on Windows |

## Install by hand

Copy the skill folders into your agent's skill directory:

* Claude Code: `.claude/skills/` in a project, or `~/.claude/skills/` for all projects.
* Codex: `.agents/skills/` in a repository, or `~/.agents/skills/` for all repositories.

For example, from the root of the repository:

```bash
mkdir -p ~/.claude/skills
cp -r skills/mlmm-* skills/colab-local-gpu-runtime ~/.claude/skills/
```

The skills leave exact flags and defaults to the installed CLI (`mlmm <subcommand> --help-advanced`) and the [command reference](reference/commands/index.md).

## See also

* [Getting Started](getting-started.md) — the shortest run, and which page to read next
* [Installation](installation.md) — setting up the environment
* [`skills/README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/skills/README.md) — the skill list on GitHub
