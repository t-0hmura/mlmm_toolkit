"""MCP (Model Context Protocol) server for mlmm_toolkit.

Exposes every `mlmm <subcmd>` CLI as an MCP tool, callable by LLM agents
(Claude Desktop, Cursor, Codeium, custom MCP clients) via JSON-RPC over
stdio.

Each tool wraps the corresponding CLI invocation via `subprocess.run`,
parses the produced `summary.json`, and returns a structured dict.

Entry point: console_script `mlmm-mcp`.
"""
from __future__ import annotations

__all__ = ["serve"]


def serve() -> None:
    """Console-script entry point.

    Lazily imports `mlmm.mcp.server` so that `mlmm --help` does not pull
    `mcp` into the import graph.
    """
    from mlmm.mcp.server import serve as _serve

    _serve()
