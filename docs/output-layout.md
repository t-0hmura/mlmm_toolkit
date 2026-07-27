# Output Directory Layout

Each `mlmm` subcommand writes to its output directory following the filename conventions below, which agents and downstream scripts can rely on.

## Filename conventions

| Filename | Written by | Purpose |
|---|---|---|
| `summary.json` | `all` and `path-search` after their summary writer is reached | Authoritative aggregate JSON envelope (see [JSON Output Reference](json-output.md)). Early CLI/input validation may fail before it exists. |
| `summary.json` | successful per-stage/report runs with `--out-json` (default `--no-out-json`); caught runtime errors may write a best-effort envelope without the flag | Compatibility mirror of leaf `result.json`. A successful writer return guarantees identical bytes. Pure utilities such as `fix-altloc`, `add-elem-info`, and `bond-summary` never emit it. |
| `result.json` | same conditions as the per-stage `summary.json` (`opt`, `tsopt`, `freq`, `irc`, `sp`, scan variants, `path-opt`, `dft`, `extract`, `trj2fig`, `energy-diagram`) | Authoritative leaf/report envelope, published after its compatibility mirror. Consume this file when distinguishing interrupted generations. |
| `summary.log` | `path-search`, `all` | Human-readable run log (one row per segment / stage). |
| `final_geometry.xyz` | `opt`, `tsopt` | Optimized geometry (XYZ, full precision). |
| `mep.pdb` / `mep.cif` / `mep_trj.xyz` | `path-search`, `all` | Reaction path frames; `mep.cif` restores original IDs for bridged input. Standalone `path-opt` writes `final_geometries_trj.xyz` / `final_geometries.pdb` instead. |
| `mep_plot.png` | `path-search`, `all` | Raw MEP energy profile (PNG). `all` copies it to the root from the engine output. |
| `forward_irc_trj.xyz` / `backward_irc_trj.xyz` (and `finished_irc_trj.xyz`) | `irc` | IRC trajectories (XYZ); companion `*_irc.pdb` files carry the same frames in PDB form. |
| `frequencies_cm-1.txt` | `freq` | Vibrational frequency listing (cm⁻¹). |
| `*.gjf` | various (when `--convert-files`) | Gaussian-format companion structure. |
| `ml_region_without_linkH.{xyz,pdb}` / `ml_region_with_linkH.{xyz,pdb}` | `all`, `dft` | Directly inspectable ML model before/after parm7-derived link-H insertion. PDB companions are written for PDB input. |

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
- **`path-search` / `path-opt` are the engine exception.** Run standalone, `path-search` is itself a deliverable (`result_path_search/` with its own `summary.log`, `mep.pdb`, optional `mep.cif`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`). Inside `all`, its raw output is engine scratch under `_work/path_opt/` (`_work/path_search/` only with `--refine-path`); the merged products (`mep.pdb`, optional `mep.cif`, `mep_trj.xyz`, `mep_plot.png`, `energy_diagram_MEP.png`) are moved to the pipeline root and `summary.{json,log}` copied there. This asymmetry is intentional.

The `all` tree therefore has three zones:

```text
result_all/
├─ summary.log · summary.json                 # copied to the root
├─ mep.pdb · mep.cif · mep_trj.xyz · mep_plot.png · energy_diagram_MEP.png
├─ energy_diagram_*_all.png · irc_plot_all.png
├─ ml_region.pdb                              # ML-region definition (reusable as --model-pdb)
├─ ml_region_without_linkH.{xyz,pdb} · ml_region_with_linkH.{xyz,pdb}
├─ mm_parm/                                   # MM topology <input>.parm7 / .rst7 (reusable as --parm)
├─ layered/                                   # layered full-system PDBs (B-factor annotated; reusable inputs)
├─ segments/
│  └─ seg_NN/                                  # 2-digit per-reactive-segment deliverables
│     ├─ reactant.{pdb,cif} · ts.{pdb,cif} · product.{pdb,cif} # CIF for bridged input
│     └─ ts/ · irc/ · freq/ · dft/ · structures/    # per-stage working files (--tsopt / --thermo / --dft)
└─ _work/                                      # pipeline scratch (safe to remove)
   ├─ pockets/ · scan/
   └─ path_opt/                                # raw MEP-engine output (path_search/ with --refine-path)
```

In TSOPT-only mode there is no MEP stage, so `_work/path_opt/` is absent and the deliverables live under `segments/seg_01/`. See [all](all.md) for the full per-mode breakdown.

## Agent recipe

```python
# Select the authoritative name for the command that produced out_dir.
import json
from pathlib import Path

subcommand = "opt"  # replace with the command you ran
primary = "summary.json" if subcommand in {"all", "path-search"} else "result.json"
summary = json.loads((Path(out_dir) / primary).read_text())

if summary["status"] == "error":
    chain = summary.get("error_class_chain", [])
    if "OptimizationError" in chain:
        # retry with looser convergence threshold
        ...
    else:
        raise RuntimeError(summary["error"])
```

`all` / `path-search` write aggregate `summary.json` after reaching their summary writer. Per-stage/report commands write `result.json` plus the mirror on a successful `--out-json` run; caught runtime exceptions may write a best-effort error envelope even without the flag. Do not assume a per-stage JSON file exists after usage validation or before its output directory is resolved.
