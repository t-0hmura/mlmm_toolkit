"""FastMCP server entry point for mlmm_toolkit.

Run as console_script:

    mlmm-mcp        # equivalent to: python -m mlmm.mcp.server

The server speaks JSON-RPC over stdio. Register in an MCP client (Claude
Desktop, Cursor, etc.) via:

    {
      "mcpServers": {
        "mlmm": {
          "command": "mlmm-mcp",
          "args": []
        }
      }
    }

See `docs/mcp_server.md` for the full tool list, parameter schemas, and
return-value structure.
"""
from __future__ import annotations

from mcp.server.fastmcp import FastMCP

from mlmm.mcp._tools import register_all

mcp = FastMCP("mlmm")
register_all(mcp)


def serve() -> None:
    """Run the MCP server on stdio (blocking call)."""
    mcp.run()


if __name__ == "__main__":
    serve()
