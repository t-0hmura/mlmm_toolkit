"""Unit tests for the MCP runner envelope (SubcmdResult / TypedDict)."""

from __future__ import annotations

from mlmm.mcp._runner import (
    SubcmdResult,
    SubcmdResultDict,
    MCP_SUBCMD_RESULT_SCHEMA_VERSION,
    MCP_SUBCMD_RESULT_STATUSES,
)


def test_subcmd_result_to_dict_carries_schema_version() -> None:
    r = SubcmdResult(status="ok", exit_code=0, argv=["mlmm", "opt"])
    d = r.to_dict()
    assert d["schema_version"] == MCP_SUBCMD_RESULT_SCHEMA_VERSION
    assert d["status"] == "ok"
    assert d["exit_code"] == 0
    assert d["argv"] == ["mlmm", "opt"]


def test_subcmd_result_status_enum() -> None:
    assert MCP_SUBCMD_RESULT_STATUSES == (
        "ok",
        "failed",
        "summary_missing",
        "summary_parse_error",
    )


def test_subcmd_result_dict_keys_match_to_dict() -> None:
    # SubcmdResultDict is a TypedDict; verify the runtime keys are a subset.
    typed_keys = set(SubcmdResultDict.__annotations__.keys())
    r = SubcmdResult(status="ok", exit_code=0)
    runtime_keys = set(r.to_dict().keys())
    assert runtime_keys <= typed_keys, (runtime_keys - typed_keys)
