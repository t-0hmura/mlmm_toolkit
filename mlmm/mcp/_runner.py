"""Shared subprocess-runner + summary-parser used by every MCP tool body.

Design intent (see docs/mcp_server.md):
- Each MCP tool body assembles a CLI argv list (`["mlmm", "opt", ...]`)
  and calls `run_subcmd(argv, out_dir)`.
- The runner spawns the CLI, captures stdout / stderr, parses the produced
  `summary.json` if any, and returns a structured `SubcmdResult`.
- Errors (non-zero exit code, missing summary.json, parse error) are surfaced
  as structured fields so the calling agent sees a stable schema, not a
  Python exception traceback.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Sequence, TypedDict


# Schema version for the MCP tool return envelope. Bump when the field
# set / value types in `SubcmdResultDict` change. 1.0 matched the
# baseline (status / exit_code / out_dir / summary / stderr_tail /
# stdout_tail / hint / argv / schema_version); 1.1 adds ``run_id`` and
# current-run summary validation.
MCP_SUBCMD_RESULT_SCHEMA_VERSION = "1.1"

# Allowed values for the `status` field. Documented in docs/mcp_server.md.
MCP_SUBCMD_RESULT_STATUSES = (
    "ok",
    "failed",
    "summary_missing",
    "summary_parse_error",
    "summary_run_mismatch",
)

_SUMMARY_ONLY_COMMANDS = frozenset({"all", "path-search"})


class SubcmdResultDict(TypedDict, total=False):
    """TypedDict shape returned by every MCP tool body.

    Mirrors :class:`SubcmdResult.to_dict`. Total=False so consumers can
    omit optional fields (out_dir / summary / hint) on early-failure paths.
    """

    schema_version: str
    status: str
    exit_code: int
    out_dir: Optional[str]
    summary: dict[str, Any]
    stderr_tail: str
    stdout_tail: str
    hint: Optional[str]
    argv: list[str]
    run_id: str


@dataclass
class SubcmdResult:
    """Structured result of a single mlmm subcmd invocation.

    `status` is one of :data:`MCP_SUBCMD_RESULT_STATUSES`. The envelope
    carries a versioned ``schema_version`` so MCP clients can pin the contract
    and migrate when the structure changes.
    """

    status: str  # see MCP_SUBCMD_RESULT_STATUSES
    exit_code: int
    out_dir: Optional[str] = None
    summary: dict[str, Any] = field(default_factory=dict)
    stderr_tail: str = ""
    stdout_tail: str = ""
    hint: Optional[str] = None
    argv: list[str] = field(default_factory=list)
    run_id: str = ""

    def to_dict(self) -> SubcmdResultDict:
        """Serialise to a plain dict the MCP framework can ship over JSON-RPC."""
        return {
            "schema_version": MCP_SUBCMD_RESULT_SCHEMA_VERSION,
            "status": self.status,
            "exit_code": self.exit_code,
            "out_dir": self.out_dir,
            "summary": self.summary,
            "stderr_tail": self.stderr_tail,
            "stdout_tail": self.stdout_tail,
            "hint": self.hint,
            "argv": self.argv,
            "run_id": self.run_id,
        }


def _tail(text: str, max_lines: int = 60) -> str:
    """Return at most the last `max_lines` lines of `text`."""
    if not text:
        return ""
    lines = text.rstrip().splitlines()
    if len(lines) <= max_lines:
        return text.rstrip()
    return "...\n" + "\n".join(lines[-max_lines:])


def _extract_hint(stderr: str) -> Optional[str]:
    """Extract the most recent `; recover: <hint>` suffix from stderr, if any.

    Subcommands emit recovery hints in this form so MCP clients can surface
    them to the agent without re-parsing the full error.
    """
    hint = None
    for line in stderr.splitlines():
        marker = "; recover:"
        if marker in line:
            hint = line.split(marker, 1)[1].strip()
    return hint


def _bind_server_argv(argv: Sequence[str]) -> list[str]:
    """Bind the product CLI alias to this interpreter and imported package."""

    requested = list(argv)
    if requested and Path(requested[0]).name == "mlmm":
        return [sys.executable, "-m", "mlmm", *requested[1:]]
    return requested


def _expected_primary_filename(argv: Sequence[str]) -> Optional[str]:
    """Return the leaf-pair primary, excluding aggregate one-file commands."""

    requested = list(argv)
    if len(requested) < 2 or Path(requested[0]).name != "mlmm":
        return None
    command = requested[1]
    return None if command in _SUMMARY_ONLY_COMMANDS else "result.json"


def _child_env(
    *,
    run_id: str,
    env_overrides: Optional[dict[str, str]],
) -> dict[str, str]:
    """Build a child environment with the imported source root first."""

    from mlmm.core.result_commit import MLMM_RUN_ID_ENV

    env = os.environ.copy()
    if env_overrides:
        env.update({str(key): str(value) for key, value in env_overrides.items()})
    source_root = str(Path(__file__).resolve().parents[2])
    inherited = env.get("PYTHONPATH", "")
    inherited_parts = [
        part
        for part in inherited.split(os.pathsep)
        if part and part != source_root
    ]
    env["PYTHONPATH"] = os.pathsep.join([source_root, *inherited_parts])
    env[MLMM_RUN_ID_ENV] = run_id
    return env


def _read_current_summary(
    summary_path: Path,
    *,
    expected_run_id: str,
) -> tuple[str, dict[str, Any], Optional[str]]:
    """Read only a well-formed summary produced by this invocation."""

    if not summary_path.exists():
        return "summary_missing", {}, None
    try:
        loaded = json.loads(summary_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return (
            "summary_parse_error",
            {},
            f"{summary_path.name} present but not valid JSON: {exc}",
        )
    if not isinstance(loaded, dict) or loaded.get("run_id") != expected_run_id:
        actual = loaded.get("run_id") if isinstance(loaded, dict) else None
        return (
            "summary_run_mismatch",
            {},
            f"{summary_path.name} run_id {actual!r} does not match current run",
        )
    return "ok", loaded, None


def run_subcmd(
    argv: Sequence[str],
    *,
    out_dir: Optional[Path] = None,
    timeout: Optional[float] = None,
    env_overrides: Optional[dict[str, str]] = None,
    summary_filename: str = "summary.json",
) -> SubcmdResult:
    """Spawn a mlmm subcmd and collect a structured result.

    Parameters
    ----------
    argv
        Argv list (e.g. ``["mlmm", "opt", "-i", "r.pdb", "-q", "-1"]``).
        The known ``mlmm`` alias is bound to this interpreter/module.
    out_dir
        Expected output directory (must match the `--out-dir` argv entry).
        Used to locate `summary.json`. If None, the runner does not parse a
        summary and only reports exit code + stderr/stdout.
    timeout
        Subprocess timeout in seconds. None = no timeout.
    env_overrides
        Optional environment variables to set for the subprocess.
    summary_filename
        Override for the summary file (default ``summary.json``).
    """
    run_id = str(uuid.uuid4())
    executed_argv = _bind_server_argv(argv)
    expected_primary_filename = _expected_primary_filename(argv)
    env = _child_env(run_id=run_id, env_overrides=env_overrides)

    try:
        proc = subprocess.run(
            executed_argv,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    except FileNotFoundError as exc:
        return SubcmdResult(
            status="failed",
            exit_code=127,
            stderr_tail=str(exc),
            hint=(
                "The requested subprocess executable is unavailable in the "
                "environment that hosts the MCP server."
            ),
            argv=executed_argv,
            run_id=run_id,
        )
    except subprocess.TimeoutExpired as exc:
        return SubcmdResult(
            status="failed",
            exit_code=124,
            stderr_tail=f"TIMEOUT after {timeout}s: {exc}",
            hint=(
                "Increase the `timeout` parameter or rerun with a smaller "
                "system / fewer cycles."
            ),
            argv=executed_argv,
            run_id=run_id,
        )

    exit_code = proc.returncode
    stderr_tail = _tail(proc.stderr)
    stdout_tail = _tail(proc.stdout)
    hint = _extract_hint(proc.stderr)

    summary: dict[str, Any] = {}
    summary_status: str = "summary_missing"
    if out_dir is not None:
        summary_path = Path(out_dir) / summary_filename
        summary_status, summary, summary_hint = _read_current_summary(
            summary_path,
            expected_run_id=run_id,
        )
        hint = hint or summary_hint
        if summary_status == "ok" and expected_primary_filename:
            primary_path = Path(out_dir) / expected_primary_filename
            primary_status, primary, primary_hint = _read_current_summary(
                primary_path,
                expected_run_id=run_id,
            )
            try:
                pair_bytes_match = (
                    primary_path.read_bytes() == summary_path.read_bytes()
                )
            except OSError:
                pair_bytes_match = False
            if primary_status != "ok" or primary != summary or not pair_bytes_match:
                summary_status = "summary_run_mismatch"
                summary = {}
                hint = hint or primary_hint or (
                    f"{expected_primary_filename} does not match the current "
                    f"{summary_filename} generation"
                )

    if exit_code != 0:
        status = "failed"
    elif out_dir is not None and summary_status != "ok":
        status = summary_status
    else:
        status = "ok"

    return SubcmdResult(
        status=status,
        exit_code=exit_code,
        out_dir=str(out_dir) if out_dir else None,
        summary=summary,
        stderr_tail=stderr_tail,
        stdout_tail=stdout_tail,
        hint=hint,
        argv=executed_argv,
        run_id=run_id,
    )
