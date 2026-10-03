# mlmm MCP server

This page explains how to call the 22 mlmm-toolkit tools from an AI agent over
MCP (Model Context Protocol): installation, the list of tools, and client
configuration. The server, `mlmm-mcp`, speaks JSON-RPC over stdio, so any
[MCP](https://modelcontextprotocol.io/) client can use it, including Claude
Desktop, Claude Code, Cursor, Codeium, and agents built on the official Python
or TypeScript MCP SDKs.

## Install

```bash
pip install "mlmm-toolkit[mcp]"
```

This adds the `mcp[cli]` dependency and registers the `mlmm-mcp` console
script.

## Tools

22 tools, one per CLI subcommand. Each tool returns a structured dict with:

- `schema_version`: version of the result format (`"2.0"`); read it from each response to know which fields to expect.
- `execution_status`: `completed` | `failed`
- `scientific_status`: `success` | `partial` | `failed`
- `summary_status`: only `ok` comes with a non-empty `summary`
  - `ok`: the `summary` of this call was read
  - `not_required`: a tool without `out_dir`, which writes no summary
  - `summary_missing`: no `summary.json` in `out_dir`
  - `summary_parse_error`: `summary.json`, or the JSON printed for `detect_bond_changes`, cannot be read as a JSON object
  - `summary_run_mismatch`: the file belongs to another run, or the `result.json` next to it is not identical
- `exit_code`: exit code of the CLI process
- `out_dir`: output directory of the stage runners and the scan / path / pipeline tools, or null for the other tools
- `summary`: parsed `summary.json`; for `detect_bond_changes`, the JSON that `mlmm bond-summary --json` prints; an empty object for the other tools without `out_dir`
- `stderr_tail` / `stdout_tail`: last ~60 lines of process output
- `hint`: the recovery hint from a `; recover: <hint>` suffix of a CLI error message, if any
- `argv`: the full command line that was run (for reproducibility)
- `run_id`: UUID of this call

The tables below list each tool's required arguments. The
input paths among them (`input_pdb`, `reactant_pdb`, and so on) are passed to
the command's `-i`; `parm7` is the full-system topology (`--parm7`). The
optional arguments set the command's CLI options, for example `charge` (`-q`),
`ligand_charge` (`-l`), and `max_cycles` (`--max-cycles`). Every tool also
takes `extra_args` (further CLI flags as a list of strings) and
`timeout_seconds`, and the stage runners and the scan / path / pipeline tools
take `out_dir`. The client receives every argument with its type in the tool's
input schema.

### Structured error envelope

When a stage runner or a scan / path tool fails, its `summary` carries the error fields below, so agents can match the error class without parsing text.

- `error`: the error message
- `error_type`: exception class name
- `error_class_chain`: the class and its parent classes, most specific first (e.g. `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`)
- `error_module`: module that defines the exception class
- `error_label`: the high-level CLI stage label

### Topology / layer prep

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `prepare_amber_topology` | `input_pdb`, `output_prefix` | `mlmm mm-parm` | AMBER parm7/rst7 for the whole system, built with AmberTools |
| `define_layer` | `input_pdb`, `output_pdb`, and `model_pdb` or `model_indices` | `mlmm define-layer` | Write the ML / Movable-MM / Frozen layers as B-factors |
| `extract_pocket` | `complex_pdb`, `ligand_id`, `radius_angstrom`, `output_pdb` | `mlmm extract` | Active-site model: residues within `radius_angstrom` of the centers given in `ligand_id` (`-c`) |

### Stage runners

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `optimize_geometry` | `input_pdb`, `parm7` | `mlmm opt` | ONIOM geometry optimization (L-BFGS by default, or RFO with microiteration) |
| `find_transition_state` | `input_pdb`, `parm7` | `mlmm tsopt` | ONIOM TS search (RS-P-RFO / Dimer / RS-I-RFO / TRIM) |
| `run_irc` | `input_pdb`, `parm7` | `mlmm irc` | ONIOM IRC integration from a TS geometry |
| `compute_frequencies` | `input_pdb`, `parm7` | `mlmm freq` | ONIOM vibrational analysis + thermochemistry |
| `run_single_point_oniom` | `input_pdb`, `parm7` | `mlmm sp` | ONIOM single-point energy + forces (+optional Hessian with `do_hess`) |

### Scans / paths / pipeline

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `scan_1d` / `scan_2d` / `scan_3d` | `input_pdb`, `parm7`, `scan_lists` | `mlmm scan` / `mlmm scan2d` / `mlmm scan3d` | ONIOM scans with harmonic restraints |
| `optimize_path` | `reactant_pdb`, `product_pdb`, `parm7` | `mlmm path-opt` | Two-endpoint ONIOM MEP optimization |
| `search_paths` | `input_pdb`, `product_pdb`, `parm7` | `mlmm path-search` | Recursive ONIOM reaction-pathway search |
| `run_full_pipeline` | `reactant_complex_pdb` | `mlmm all` | End-to-end: extract → MEP → TS → IRC → freq → DFT |
| `run_single_point_dft` | `input_pdb`, `parm7` | `mlmm dft` | Single-point DFT of the ML region, combined with the MM energy (GPU4PySCF or PySCF) |

### ONIOM I/O (Gaussian / ORCA)

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `export_oniom_input` | `input_layered_pdb`, `parm7`, `charge`, `multiplicity`, `output_file` | `mlmm oniom-export` | Write a Gaussian g16 or ORCA ONIOM input (`format_engine`: `g16` by default, or `orca`) |
| `import_oniom_input` | `input_file`, `output_prefix` | `mlmm oniom-import` | Read a Gaussian / ORCA ONIOM input back to XYZ and a layered PDB |

### Structure / I/O helpers

| MCP tool | Required arguments | CLI subcmd | Purpose |
|---|---|---|---|
| `add_element_info` | `input_pdb`, `output_pdb` | `mlmm add-elem-info` | Repair PDB element columns |
| `fix_altloc` | `input_pdb`, `output_pdb` | `mlmm fix-altloc` | Resolve PDB alternate locations |
| `plot_trajectory` | `input_trj_xyz`, `output_png` | `mlmm trj2fig` | Energy profile figure (PNG; also JPEG/SVG/PDF/HTML/CSV) |
| `plot_energy_diagram` | `energies`, `output_png` | `mlmm energy-diagram` | State energy diagram from given energies |
| `detect_bond_changes` | `reactant_pdb`, `product_pdb` | `mlmm bond-summary` | Bond changes between two structures (XYZ / PDB / GJF) |

### Charge and ordered inputs

`charge` is the charge of the ML region (`-q`), not of the whole system. To
derive it from residue names instead, omit `charge` and give a per-resname
`ligand_charge` mapping such as `"SAM:1,GPP:-3"`; the sum then covers the ML
region taken from the B-factor layers of the input or from `--model-pdb`, which
`run_single_point_oniom` takes as `model_pdb` and the other tools through
`extra_args` ({ref}`ML/MM options <mlmm-options>`).

`search_paths` takes any ordered intermediates between `input_pdb` (reactant)
and `product_pdb` as `intermediate_pdbs`. `scan_1d`, `scan_2d`, and `scan_3d`
pass `scan_lists` as one value of `--scan-lists`; the [`scan`](scan.md) page
gives its format. With only `reactant_complex_pdb`, `run_full_pipeline` needs
`do_tsopt=True` or a scan given through `extra_args`
(`["--scan-lists", "…"]`), as `mlmm all` does with one input.

## IRC and TS optimization settings

The IRC and TS arguments are the CLI options of the same name; the command
pages give their meaning and defaults.

- `run_irc` (`step_size`, `irc_pos_def`): [`irc`](irc.md); `--irc-pos-def` is in the [generated reference](reference/commands/irc.md)
- `find_transition_state` (`opt_mode`, `microiter`, `flatten`): [`tsopt`](tsopt.md) `--opt-mode`, default `hess` (RS-P-RFO); `microiter=False` turns off microiteration; see {ref}`--opt-mode by command <opt-mode-semantics>`
- `run_full_pipeline` (`refine_path`, `do_tsopt`, `do_thermo`, `do_dft`, `thresh_post`): [`all`](all.md) `--refine-path`, `--tsopt`, `--thermo`, and `--dft`; `--thresh-post` is in the [generated reference](reference/commands/all.md)

## Client configuration

Client configuration schemas differ. The snippet below applies to clients that
accept a top-level `mcpServers` object; check the client's own MCP
documentation before choosing the config file and schema.

- Claude Desktop — `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) / `%APPDATA%\Claude\claude_desktop_config.json` (Windows)
- Cursor — `~/.cursor/mcp.json`
- Claude Code — no file to edit: run `claude mcp add mlmm -- mlmm-mcp`

Once the client has started the server, its tool list shows the 22 tools.

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
for a full example that sets environment variables (PATH / AMBERHOME /
CUDA_VISIBLE_DEVICES).

VS Code instead uses a top-level `servers` object in `.vscode/mcp.json`
([VS Code MCP configuration reference](https://code.visualstudio.com/docs/agents/reference/mcp-configuration)):

```json
{
  "servers": {
    "mlmm": {
      "command": "mlmm-mcp",
      "args": []
    }
  }
}
```

### Custom Python MCP client

```python
import asyncio

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

async def main():
    server_params = StdioServerParameters(command="mlmm-mcp")
    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            result = await session.call_tool(
                "optimize_geometry",
                arguments={
                    "input_pdb": "r_complex_layered.pdb",
                    "parm7": "real.parm7",
                    "charge": 0,
                    "max_cycles": 50,
                },
            )
            print(result.content)

asyncio.run(main())
```

## Sandbox / safety notes

- Each tool runs `mlmm` with the server's own Python interpreter
  (`python -m mlmm`) in a subprocess that keeps the caller's working directory,
  so relative input paths keep their meaning and another `mlmm` earlier on
  `PATH` is not used.
- The server inherits the caller's PATH, conda environment, CUDA setup, and
  AmberTools path; `prepare_amber_topology` needs AmberTools (antechamber,
  parmchk2, tleap) on PATH. Long tools (opt / tsopt / irc / scan) run the CLI
  in a subprocess, so set `timeout_seconds` on each call to stop runaway
  calculations (default: no timeout).
- The stage runners and the scan / path / pipeline tools write under `out_dir`.
  When it is not given, each call gets its own temporary directory
  (e.g. `mlmm_mcp_opt_…`), so parallel calls do not collide.
- The other tools have no `out_dir` and write to their explicit output path.
  `extra_args` passes additional CLI flags but cannot override the typed output
  paths, `--out-dir`, `--out-json/--no-out-json`, or `--json/--no-json` for
  `detect_bond_changes`, also in the attached forms such as `--out-dir=…`;
  such a call is rejected before the CLI starts. The
  returned `argv` shows every path the command was given.
- The server is not a filesystem sandbox: external programs and model caches
  can write outside `out_dir`. Input structures, parm7 files, and MLIP weights
  must already be on disk; use client-side path restrictions or OS / container
  isolation when file access must be confined.

## Notes

- For `run_full_pipeline`, the tools without `out_dir`, and a run that stops before it writes `summary.json` (`summary_missing`), read `stderr_tail` and `hint` instead of the structured error fields.
- The `ligand_charge` derivation does not work with `--model-indices`, or with `--no-detect-layer` and no `--model-pdb`; give `charge` in those cases.
- If both `charge` and `ligand_charge` are given, the explicit `charge` wins.

## See also

* [JSON Output Reference](json-output.md) — the status fields and `summary.json` that the tools return
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* [Command Reference](reference/commands/index.md) — the CLI options behind each tool, for `extra_args`
