# mlmm MCP server

`mlmm-mcp` is an [MCP](https://modelcontextprotocol.io/) server that lets any
MCP-speaking agent drive every `mlmm` CLI subcommand via JSON-RPC over stdio
— including Claude Desktop / Claude Code / Cursor / Codeium, and any custom
agent built on the official Python or TypeScript MCP SDKs.

## Install

```bash
pip install "mlmm-toolkit[mcp]"
```

This adds the `mcp[cli]` dependency and registers the `mlmm-mcp` console
script.

## Tools

22 tools, one per CLI subcommand. Each tool returns a structured dict (`SubcmdResultDict` in `mlmm.mcp._runner`) with:

- `schema_version`: envelope version. Live value: `mlmm.mcp._runner.MCP_SUBCMD_RESULT_SCHEMA_VERSION`. Bumps signal a field-set / value-type change; pin against the constant rather than the literal in this doc.
- `status`: `ok` | `failed` | `summary_missing` | `summary_parse_error`
- `exit_code`: subprocess exit code
- `out_dir`: working directory the CLI wrote to
- `summary`: parsed `summary.json` (CLI output schema; see [JSON Output Reference](json-output.md) for the per-stage shape)
- `stderr_tail` / `stdout_tail`: last ~60 lines of process output
- `hint`: parsed `; recover: <hint>` suffix from CLI error messages, if any
- `argv`: the full argv that was executed (for reproducibility)

For typed-Python consumers, `mlmm.mcp._runner` also exposes `SubcmdResultDict` (a `TypedDict` mirroring the runtime payload) and `MCP_SUBCMD_RESULT_STATUSES` (the enum tuple of allowed `status` strings).

### Structured error envelope

When a subcommand fails, the parsed `summary` (or sibling `result.json`) carries an extended error envelope so agents can pattern-match the exception class hierarchy without parsing text:

- `error`: `str(exc)` of the original exception
- `error_type`: exception class name (e.g. `"OptimizationError"`)
- `error_class_chain`: MRO class names (e.g. `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`)
- `error_module`: defining module of the exception class
- `error_label`: the high-level CLI stage label (e.g. `"opt"`, `"tsopt-stage"`)

### Topology / layer prep (mlmm-specific)

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `prepare_amber_topology` | `mlmm mm-parm` | Generate AMBER parm7/rst7 via AmberTools |
| `define_layer` | `mlmm define-layer` | Assign ML / MM-movable / MM-frozen B-factor layers |
| `extract_pocket` | `mlmm extract` | Cut a sphere around a ligand to make an active-site model |

### Stage runners (ONIOM-aware)

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `optimize_geometry` | `mlmm opt` | ONIOM geometry opt (microiter macro / micro) |
| `find_transition_state` | `mlmm tsopt` | ONIOM TS search (RS-I-RFO / Dimer / TRIM / RS-P-RFO) |
| `run_irc` | `mlmm irc` | ONIOM IRC integration from a TS geometry |
| `compute_frequencies` | `mlmm freq` | ONIOM vibrational analysis + thermochemistry |
| `run_single_point_oniom` | `mlmm sp` | ONIOM single-point energy + forces (+optional Hessian) |

### Scans / paths / pipeline

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `mlmm scan{,2d,3d}` | ONIOM restraint scans |
| `optimize_path` | `mlmm path-opt` | Two-endpoint ONIOM MEP optimization |
| `search_paths` | `mlmm path-search` | Recursive ONIOM pathway search |
| `run_full_pipeline` | `mlmm all` | End-to-end: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `mlmm dft` | ONIOM-embedded single-point DFT via gpu4pyscf |

### ONIOM I/O (Gaussian / ORCA)

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `export_oniom_input` | `mlmm oniom-export` | Write a Gaussian g16 / ORCA ONIOM input deck |
| `import_oniom_input` | `mlmm oniom-import` | Read a Gaussian / ORCA ONIOM input back to XYZ / layered PDB |

### Structure / I/O helpers

| MCP tool | CLI subcmd | Purpose |
|---|---|---|
| `add_element_info` | `mlmm add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `mlmm fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `mlmm trj2fig` | Energy profile PNG / HTML / SVG / PDF |
| `plot_energy_diagram` | `mlmm energy-diagram` | Categorical energy diagram |
| `detect_bond_changes` | `mlmm bond-summary` | Bond-change diff between two PDBs |

## Opt-in IRC convergence guard

`run_irc` accepts `irc_pos_def: bool` — IRC convergence then additionally
requires a positive-definite mass-weighted Hessian, blocking the IRC
"shoulder" false-convergence where the rms-only criterion calls success
before reaching the local minimum. Defaults to `None` (rms-only, legacy).

`find_transition_state` accepts `opt_mode="trim"` (Helgaker 1991) /
`opt_mode="rsprfo"` (Banerjee 1985) as alternative TS optimizers; the
server passes `--no-microiter` automatically since those modes are not
microiter-capable.

## Client configuration

Every MCP client takes the same `mcpServers` schema. Drop the snippet below
into your client's MCP config file (Claude Desktop's
`~/Library/Application Support/Claude/claude_desktop_config.json` on macOS,
Cursor's `~/.cursor/mcp.json`, etc.):

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

See [`examples/mcp_client_config.json`](../examples/mcp_client_config.json)
for a full example with explicit env-var overrides (PATH / AMBERHOME /
CUDA_VISIBLE_DEVICES).

## Sandbox / safety notes

- The MCP server inherits the calling environment's PATH, conda env, CUDA
  setup, and AmberTools path. Each tool spawns the `mlmm` CLI as a subprocess,
  so long-running tools (opt / tsopt / irc / scan) run out-of-process — set
  `timeout_seconds` on each call to bound them.
- Output files land under the `out_dir` kwarg (defaults to a unique
  `tempfile.mkdtemp("mlmm_mcp_<subcmd>_…")`).
- The server does not modify `~/.bashrc` / login env, install software,
  or write outside `out_dir`. Required inputs (PDBs, parm7, ML weights)
  must already exist on disk.
