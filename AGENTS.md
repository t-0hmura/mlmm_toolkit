# Working on mlmm-toolkit

This file is for coding agents that change the code, docs, or skills in this repository.
To run mlmm-toolkit for a user, load the skills in `skills/` instead (see `skills/README.md`).
The full rules and their reasons are in `CONTRIBUTING.md`; § numbers below refer to it.

## Before you finish

1. `pytest -q`. Fix a failing test at its cause; never delete or skip it.
2. `python .github/scripts/check_engineering_markers.py`.
3. `python .github/scripts/run_docs_quality.py`. It stops at the first failing step; to see every failure, run the scripts it lists one by one.
4. `python .github/scripts/check_help_registry.py` after CLI changes.
5. Smoke: copy `tests/smoke/` outside the repository and run `bash run.sh` through your scheduler.

After changing CLI options, regenerate `docs/reference/` with `python .github/scripts/generate_reference.py`.
Write boolean options as `--flag/--no-flag`.

## Do not touch

Read §4 before changing any of these:

- The nine chemistry rules marked `# CHEMISTRY-RULE:N`, and `# DOMAIN_PURE` modules (§4.1).
- The `del` chains that release GPU memory (§4.2).
- Divergent files in the bundled forks, including `hessian_ff/`, which has no upstream package; read its `README.md` first (§4.3).
- Package discovery and runtime dependencies (§4.4).
- Absolute paths in `_LAZY_SUBCOMMANDS` (§4.5).
- Chemistry default choices, including the ONIOM region shell radii (§4.6).
- Output that downstream parsers read, including `summary.json` (§4.7, §1.5).

## Editing skills

Each skill is `skills/<name>/SKILL.md`.
After editing skills, run `python .github/scripts/run_docs_quality.py` (it runs the skill checkers) and `pytest tests/test_skill_command_checker.py tests/test_skill_drift_checker.py tests/test_dynamic_dispatch_skill.py tests/test_cli_completion.py -q`.
The checkers require:

- Frontmatter with `name` and `description`; `name` is lowercase hyphen-case, at most 64 characters, and equals the folder name.
- A `description` of at most 1024 characters, without `<` or `>`.
- Real flags in shell code blocks (no language, `bash`, `console`, or `sh`) and in inline `mlmm …` (`check_skill_commands.py`).
- In a concrete `mlmm` example with `-i`: the required topology option (`--parm7`); a charge option (`-q` or `-l`) unless every input is `.gjf` or `--config` is given; and `--ref-pdb` for XYZ input (for `path-search`, one per XYZ position). Write placeholders as `<...>`.
- Real flags in backticks in prose in the folders listed in `MLMM_CLI_DIRS` (`check_skill_drift.py`). Add a new skill folder that documents the CLI there.
- The sentences that `check_docs_contract.py` pins in `docs/` and `skills/` (some in a fixed order), and none of the phrases it bans. `tests/test_dynamic_dispatch_skill.py` also reads the `## dispatcher.sh` section of `skills/mlmm-hpc/dynamic-dispatch.md`. When you move or reword pinned text, update the checker or test in the same change.

Also write only `name` and `description` in the frontmatter, state TRIGGER and SKIP in the `description`, and keep each `SKILL.md` within 500 lines.

Do not add flag tables to skills; point to `--help-advanced` and `docs/reference/commands/`. Describe current behaviour only.
When you add, remove, or rename a skill, update `skills/README.md` and the `[Unreleased]` section of `CHANGELOG.md` in the same change.

## More

[CONTRIBUTING.md](CONTRIBUTING.md) · [Architecture](docs/architecture.md) · [Skills index](skills/README.md)
