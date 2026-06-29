---
name: mlmm-mcp
description: How to drive `mlmm-toolkit` from an MCP-speaking agent (Claude Desktop / Claude Code / Cursor / Codeium / custom Python or TypeScript MCP SDK clients) via the bundled `mlmm-mcp` server. Lists the 22 MCP tools (one per CLI subcommand, including mlmm-specific topology / ONIOM-layer / ONIOM-input tools) and the shared `SubcmdResult` schema. TRIGGER on questions about MCP setup, agent integration, tool names, return-value shape, or "how does an agent invoke mlmm". SKIP for direct CLI usage — that's `mlmm-cli`.
---

# mlmm-toolkit MCP server

`mlmm-mcp` exposes every CLI subcommand as an MCP tool over stdio JSON-RPC. Any MCP client can drive the full ONIOM reaction-path pipeline without spawning shell processes manually.

## Install + register

```bash
pip install "mlmm[mcp]"   # adds the `mcp[cli]` extra
```

This registers the `mlmm-mcp` console script. Drop the snippet below into your client's MCP config file:

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) / `%APPDATA%\Claude\claude_desktop_config.json` (Windows)
- Cursor — `~/.cursor/mcp.json`
- Other MCP clients — consult the client's own MCP-server docs

```json
{
  "mcpServers": {
    "mlmm": {
      "command": "mlmm-mcp",
      "args": []
    }
  }
}
```

See [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json) for a full sample with explicit `PATH` / `AMBERHOME` / `CUDA_VISIBLE_DEVICES` env-var overrides — AmberTools (antechamber / parmchk2 / tleap) must be on `PATH` for the `prepare_amber_topology` tool to work.

## 22 tools (1 per CLI subcommand)

### Topology / layer prep (mlmm-specific)

| MCP tool | wraps | purpose |
|---|---|---|
| `prepare_amber_topology` | `mlmm mm-parm` | Generate AMBER parm7 / rst7 via AmberTools |
| `define_layer` | `mlmm define-layer` | Assign ML / MM-movable / MM-frozen B-factor layers |
| `extract_pocket` | `mlmm extract` | Cut a sphere around a ligand to make an active-site model |

### Stage runners (ONIOM-aware)

| MCP tool | wraps | purpose |
|---|---|---|
| `optimize_geometry` | `mlmm opt` | ONIOM geometry opt (microiter macro / micro) |
| `find_transition_state` | `mlmm tsopt` | ONIOM TS search (RS-I-RFO / Dimer / TRIM / RS-P-RFO) |
| `run_irc` | `mlmm irc` | ONIOM IRC integration from a TS geometry |
| `compute_frequencies` | `mlmm freq` | ONIOM vibrational analysis + thermochemistry |
| `run_single_point_oniom` | `mlmm sp` | ONIOM single-point energy + forces (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | wraps | purpose |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `mlmm scan{,2d,3d}` | ONIOM restraint scans |
| `optimize_path` | `mlmm path-opt` | Two-endpoint ONIOM MEP optimisation |
| `search_paths` | `mlmm path-search` | Recursive ONIOM pathway search |
| `run_full_pipeline` | `mlmm all` | End-to-end (extract → MEP → TS → IRC → freq → DFT) |
| `run_single_point_dft` | `mlmm dft` | ONIOM-embedded single-point DFT via gpu4pyscf |

### ONIOM I/O (Gaussian / ORCA)

| MCP tool | wraps | purpose |
|---|---|---|
| `export_oniom_input` | `mlmm oniom-export` | Write a Gaussian g16 / ORCA ONIOM input deck |
| `import_oniom_input` | `mlmm oniom-import` | Read a Gaussian / ORCA ONIOM input back to XYZ / layered PDB |

### Structure / I/O helpers

| MCP tool | wraps | purpose |
|---|---|---|
| `add_element_info` | `mlmm add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `mlmm fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `mlmm trj2fig` | Energy profile PNG / HTML / SVG / PDF |
| `plot_energy_diagram` | `mlmm energy-diagram` | Categorical energy diagram |
| `detect_bond_changes` | `mlmm bond-summary` | Bond-change diff between two PDBs |

## `SubcmdResult` return schema

Every tool returns the same structured dict so the calling agent can dispatch on `status` without parsing stderr:

```python
{
    "schema_version": "1.0",         # pin to MCP_SUBCMD_RESULT_SCHEMA_VERSION
    "status": "ok" | "failed" | "summary_missing" | "summary_parse_error",
    "exit_code": int,                # subprocess exit code
    "out_dir": str | None,           # working directory the CLI wrote to
    "summary": dict,                 # parsed summary.json (CLI output schema)
    "stderr_tail": str,              # last ~60 lines of stderr
    "stdout_tail": str,              # last ~60 lines of stdout
    "hint": str | None,              # parsed `; recover: <hint>` suffix
    "argv": list[str],               # full argv that was executed (reproducible)
}
```

A failed subcommand additionally surfaces a structured exception envelope inside `summary`:

- `error_class_chain` — list of class names walking the MRO (e.g. `["CudaOutOfMemoryError", "RuntimeError", "Exception"]`)
- `error_module` — module path of the originating exception class
- `error_label` — the high-level CLI stage label (e.g. `"optimization"`, `"TS optimization"`, `"IRC"`, `"single-point"`, or the fallback `"UnhandledError"`). This is a coarse stage label; conditions such as OOM must be read from `error_class_chain` (e.g. a CUDA out-of-memory class), not from `error_label`

so an MCP client can pattern-match on the hierarchy instead of substring-matching `stderr_tail`.

## Opt-in IRC convergence guard

`run_irc` accepts `irc_pos_def: bool` — IRC convergence then additionally requires a positive-definite mass-weighted Hessian, blocking the IRC "shoulder" false-convergence where the rms-only criterion calls success before reaching the local minimum. Defaults to `None` (rms-only, legacy).

`find_transition_state` accepts `opt_mode="trim"` (Helgaker 1991) / `opt_mode="rsprfo"` (Banerjee 1985) as alternative TS optimisers; the server passes `--no-microiter` automatically since those modes are not microiter-capable.

## Sandbox / safety notes

- The MCP server inherits the calling environment's PATH, conda env, CUDA setup, and AmberTools path. Long-running tools (opt / tsopt / irc / scan) launch the `mlmm` CLI as a subprocess — set `timeout_seconds` on each call.
- Output files land under the `out_dir` kwarg (defaults to a unique `tempfile.mkdtemp("mlmm_mcp_<subcmd>_…")`).
- The server does not modify `~/.bashrc` / login env, install software, or write outside `out_dir`. Required inputs (PDBs, parm7, ML weights) must already exist on disk.

## See also
- Full MCP server doc: [`docs/mcp_server.md`](../../docs/mcp_server.md)
- Sample MCP-client config: [`examples/mcp_client_config.json`](../../examples/mcp_client_config.json)
- Server / tool source: `mlmm/mcp/{server,_tools,_runner}.py`
