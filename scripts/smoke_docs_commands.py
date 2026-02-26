#!/usr/bin/env python3
"""Lightweight smoke tests for commands embedded in docs markdown files."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from click.testing import CliRunner

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPO_ROOT / "docs"
TOOL_NAME = "mlmm"
CLI_MODULE = "mlmm_toolkit"
DOCS_SMOKE_COMMAND_TIMEOUT_SEC = float(os.environ.get("DOCS_SMOKE_COMMAND_TIMEOUT_SEC", "120"))

sys.path.insert(0, str(REPO_ROOT))

from mlmm_toolkit.cli import cli as root_cli  # noqa: E402


_CODE_LANGS = {"", "bash", "sh", "shell", "console"}
_ALL_ONLY_PATH_EXTS = {".pdb", ".xyz", ".gjf", ".yaml", ".yml", ".json"}


def _extract_commands_from_block(lines: list[str]) -> list[str]:
    commands: list[str] = []
    current = ""
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if " #" in line:
            line = line.split(" #", 1)[0].rstrip()
            if not line:
                continue
        if line.startswith("$"):
            line = line[1:].strip()
        current = f"{current} {line}".strip() if current else line
        if current.endswith("\\"):
            current = current[:-1].rstrip()
            continue
        commands.append(current)
        current = ""
    if current:
        commands.append(current)
    filtered: list[str] = []
    for cmd in commands:
        if not cmd.startswith(TOOL_NAME):
            continue
        if any(mark in cmd for mark in ("<", ">", "[", "]")):
            continue
        filtered.append(cmd)
    return filtered


def _extract_docs_commands() -> list[str]:
    commands: list[str] = []
    for path in sorted(DOCS_ROOT.rglob("*.md")):
        lines = path.read_text(encoding="utf-8").splitlines()
        in_fence = False
        fence_lang = ""
        block: list[str] = []
        for line in lines:
            stripped = line.strip()
            if stripped.startswith("```"):
                marker = stripped[3:].strip().lower()
                if not in_fence:
                    in_fence = True
                    fence_lang = marker
                    block = []
                else:
                    if fence_lang in _CODE_LANGS:
                        commands.extend(_extract_commands_from_block(block))
                    in_fence = False
                    fence_lang = ""
                    block = []
                continue
            if in_fence:
                block.append(line)
    return commands


def _subcommand_from_tokens(tokens: list[str]) -> str:
    if len(tokens) < 2 or tokens[1].startswith("-"):
        return "all"
    return tokens[1]


def _prepare_fixture_files(tmp: Path) -> dict[str, Path]:
    pdb_text = (
        "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    r_pdb = tmp / "R.pdb"
    p_pdb = tmp / "P.pdb"
    xyz = tmp / "input.xyz"
    gjf = tmp / "input.gjf"
    cfg = tmp / "config.yaml"
    out_dir = tmp / "result_all"

    r_pdb.write_text(pdb_text, encoding="utf-8")
    p_pdb.write_text(pdb_text, encoding="utf-8")
    xyz.write_text("1\n\nC 0.0 0.0 0.0\n", encoding="utf-8")
    gjf.write_text("%chk=test\n#p hf/3-21g\n\nTitle\n\n0 1\nC 0.0 0.0 0.0\n\n", encoding="utf-8")
    cfg.write_text("extract:\n  radius: 2.6\n", encoding="utf-8")
    out_dir.mkdir(parents=True, exist_ok=True)

    return {
        "r_pdb": r_pdb,
        "p_pdb": p_pdb,
        "xyz": xyz,
        "gjf": gjf,
        "config": cfg,
        "out_dir": out_dir,
    }


def _sanitize_all_args(args: list[str], fixture: dict[str, Path]) -> list[str]:
    out: list[str] = []
    saw_input = False
    saw_dry_run = False
    saw_center = False
    i = 0
    while i < len(args):
        tok = args[i]
        if tok in {"--version", "-h", "--help"}:
            i += 1
            continue
        if tok in {"-c", "--center"}:
            saw_center = True
            out.append(tok)
            if i + 1 < len(args):
                out.append(args[i + 1])
            i += 2
            continue
        if tok in {"-i", "--input"}:
            saw_input = True
            out.extend([tok, str(fixture["r_pdb"]), str(fixture["p_pdb"])])
            i += 1
            while i < len(args) and not args[i].startswith("-"):
                i += 1
            continue
        if tok == "--config":
            out.extend([tok, str(fixture["config"])])
            i += 2
            continue
        if tok == "--out-dir":
            out.extend([tok, str(fixture["out_dir"])])
            i += 2
            continue
        if tok == "--dry-run":
            saw_dry_run = True
            out.append(tok)
            i += 1
            continue
        if tok == "--no-dry-run":
            saw_dry_run = True
            out.append("--dry-run")
            i += 1
            continue
        if (not tok.startswith("-")) and Path(tok).suffix.lower() in _ALL_ONLY_PATH_EXTS:
            ext = Path(tok).suffix.lower()
            if ext == ".pdb":
                out.append(str(fixture["r_pdb"]))
            elif ext == ".xyz":
                out.append(str(fixture["xyz"]))
            elif ext == ".gjf":
                out.append(str(fixture["gjf"]))
            else:
                out.append(str(fixture["config"]))
            i += 1
            continue

        out.append(tok)
        i += 1

    if not saw_input:
        out.extend(["-i", str(fixture["r_pdb"]), str(fixture["p_pdb"])])
    if not saw_center:
        out.extend(["-c", "LIG"])
    if not saw_dry_run:
        out.append("--dry-run")
    if "--out-dir" not in out:
        out.extend(["--out-dir", str(fixture["out_dir"])])
    return out


def _run_help_smoke(commands: list[str]) -> None:
    runner = CliRunner()
    subcommands = sorted({_subcommand_from_tokens(shlex.split(cmd)) for cmd in commands})
    for subcmd in subcommands:
        result = runner.invoke(root_cli, [subcmd, "--help"], catch_exceptions=False)
        if result.exit_code != 0:
            raise RuntimeError(
                f"[help-smoke] failed for '{TOOL_NAME} {subcmd} --help':\n{result.output}"
            )
    print(f"[help-smoke] validated {len(subcommands)} subcommands from docs.")


def _run_all_dry_run_smoke(commands: list[str]) -> None:
    try:
        probe = subprocess.run(
            [sys.executable, "-m", CLI_MODULE, "all", "--help"],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            "[dry-run-smoke] timeout while probing availability for "
            f"'{TOOL_NAME} all --help' ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
        ) from exc
    probe_output = f"{probe.stdout}\n{probe.stderr}"
    if "Command 'all' is unavailable" in probe_output or "Missing dependency:" in probe_output:
        print("[dry-run-smoke] skipped: 'all' command is unavailable in this environment.")
        return

    all_cmds: set[str] = set()
    for cmd in commands:
        tokens = shlex.split(cmd)
        if not tokens or tokens[0] != TOOL_NAME:
            continue
        if len(tokens) >= 2 and tokens[1] == "all":
            all_cmds.add(cmd)
            continue
        if len(tokens) >= 2 and tokens[1].startswith("-"):
            if any(tok in {"-i", "--input"} for tok in tokens[1:]):
                all_cmds.add(cmd)
    all_cmds = sorted(all_cmds)
    if not all_cmds:
        raise RuntimeError("No 'all' command examples found in docs.")

    with tempfile.TemporaryDirectory(prefix=f"{TOOL_NAME}_docs_smoke_") as tmpdir:
        fixture = _prepare_fixture_files(Path(tmpdir))
        for raw in all_cmds:
            tokens = shlex.split(raw)
            if not tokens or tokens[0] != TOOL_NAME:
                continue
            args = tokens[1:]
            if not args or args[0].startswith("-"):
                args = ["all", *args]
            if args[0] != "all":
                continue
            dry_args = _sanitize_all_args(args, fixture)
            try:
                completed = subprocess.run(
                    [sys.executable, "-m", CLI_MODULE, *dry_args],
                    cwd=REPO_ROOT,
                    text=True,
                    capture_output=True,
                    timeout=DOCS_SMOKE_COMMAND_TIMEOUT_SEC,
                )
            except subprocess.TimeoutExpired as exc:
                raise RuntimeError(
                    f"[dry-run-smoke] timeout for docs command ({DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s):\n"
                    f"  {raw}\n"
                    f"sanitized args:\n"
                    f"  {TOOL_NAME} {' '.join(dry_args)}"
                ) from exc
            if completed.returncode != 0:
                raise RuntimeError(
                    f"[dry-run-smoke] failed for docs command:\n  {raw}\n"
                    f"sanitized args:\n  {TOOL_NAME} {' '.join(dry_args)}\n\n"
                    f"stdout:\n{completed.stdout}\n\nstderr:\n{completed.stderr}"
                )
    print(
        f"[dry-run-smoke] validated {len(all_cmds)} docs examples "
        f"(timeout={DOCS_SMOKE_COMMAND_TIMEOUT_SEC:g}s)."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()

    commands = _extract_docs_commands()
    if not commands:
        raise RuntimeError("No commands were extracted from docs markdown code fences.")

    _run_help_smoke(commands)
    _run_all_dry_run_smoke(commands)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
