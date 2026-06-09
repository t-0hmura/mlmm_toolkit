# Output Directory Layout

Each `mlmm` subcommand writes to its output directory under the filename conventions below, which agents and downstream scripts can rely on.

## Filename conventions

| Filename | Written by | Purpose |
|---|---|---|
| `summary.json` | `all`, `path-search`, and every per-stage subcommand that runs `write_result_json` | Authoritative JSON envelope (see [JSON Output Reference](json-output.md)). Read this first. Pure utility subcommands that do not call `write_result_json` (e.g. `fix-altloc`, `add-elem-info`, `bond-summary`) do not emit it. |
| `result.json` | per-stage subcommands that call `write_result_json` (`opt`, `tsopt`, `freq`, `irc`, `sp`, `scan` / `scan2d` / `scan3d`, `path-opt`, `dft`, `extract`) | Alternate filename — identical payload to `summary.json`. Read `summary.json` for the single-filename convention; `result.json` carries the same content and can be deleted if you only consume `summary.json`. |
| `summary.log` | `path-search`, `all` | Human-readable run log (one row per segment / stage). |
| `final_geometry.xyz` | `opt`, `tsopt` | Optimized geometry (XYZ, full precision). |
| `mep.pdb` / `mep_trj.xyz` | `path-search`, `all` | Reaction path frames (PDB / XYZ); standalone `path-opt` writes `final_geometries_trj.xyz` / `final_geometries.pdb` instead. |
| `mep_plot.png` | `path-search`, `all` | Raw MEP energy profile (PNG). `all` copies it to the root from the engine output. |
| `forward_irc_trj.xyz` / `backward_irc_trj.xyz` (and `finished_irc_trj.xyz`) | `irc` | IRC trajectories (XYZ); companion `*_irc.pdb` files carry the same frames in PDB form. |
| `frequencies_cm-1.txt` | `freq` | Vibrational frequency listing (cm⁻¹). |
| `*.gjf` | various (when `--convert-files`) | Gaussian-format companion structure. |

## Default `--out-dir`

| Subcommand | Default `--out-dir` |
|---|---|
| `all` | `./result_all/` |
| `opt` | `./result_opt/` |
| `tsopt` | `./result_tsopt/` |
| `freq` | `./result_freq/` |
| `irc` | `./result_irc/` |
| `dft` | `./result_dft/` |
| `scan` / `scan2d` / `scan3d` | `./result_scan*/` |
| `path-opt` / `path-search` | `./result_path_*/` |
| `sp` | `./result_sp/` |
| `extract` | `./` (writes `pocket.pdb`, or `pocket_<input>.pdb` for multiple inputs, in the working directory) |
| `mm-parm` | `./` (writes `<prefix>.parm7` / `<prefix>.rst7`) |
| `define-layer` | `./` (writes `<input>_layered.pdb`) |

Override with `--out-dir <path>` (or `-o`); explicit paths take precedence over both per-stage defaults and YAML.

## Standalone vs `all`

A subcommand run on its own writes a **flat** result directory. The same writer, when orchestrated by `all`, nests into a structured tree:

- **Standalone subcommand** → flat `result_<subcmd>/` with the files above. There is no `segments/` and no `_work/` — those appear only when `all` coordinates several writers in one run.
- **Inside `all`, leaf writers nest unchanged.** A per-segment leaf output at `segments/seg_NN/<subcmd>/` is structurally identical to the standalone `result_<subcmd>/`; `all` just points the writer at a different directory.
- **`path-search` / `path-opt` are the engine exception.** Run standalone, `path-search` is itself a deliverable (`result_path_search/` with its own `summary.log`, `mep.pdb`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`). Inside `all`, its raw output is engine scratch under `_work/path_search/`; the merged products (`mep.pdb`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`) are moved to the pipeline root and `summary.{json,log}` copied there. This asymmetry is intentional.

The `all` tree therefore has three zones:

```text
result_all/
├─ summary.log · summary.json                 # copied to the root
├─ mep.pdb · mep_trj.xyz · mep_plot.png · energy_diagram_MEP.png   # MEP products moved from the engine
├─ energy_diagram_*_all.png · irc_plot_all.png
├─ ml_region.pdb                              # ML-region definition (reusable as --model-pdb)
├─ mm_parm/                                   # MM topology <input>.parm7 / .rst7 (reusable as --parm)
├─ layered/                                   # layered full-system PDBs (B-factor annotated; reusable inputs)
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.pdb · ts.pdb · product.pdb         # canonical R/TS/P
│     └─ ts/ · irc/ · freq/ · dft/ · structures/    # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to remove)
   ├─ pockets/ · scan/
   └─ path_search/                             # raw MEP-engine output (path_opt/ with --no-refine-path)
```

In TSOPT-only mode there is no MEP stage, so `_work/path_search/` is absent and the deliverables live under `segments/seg_01/`. See [all](all.md) for the full per-mode breakdown.

## Agent recipe

```python
# Read whichever subcommand's output, single filename across the board.
import json
from pathlib import Path

summary = json.loads((Path(out_dir) / "summary.json").read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

`summary.json` is guaranteed to exist on every subcommand that called `write_result_json` (including failure paths — failure envelopes also carry the schema version + error class chain).