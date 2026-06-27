# `mlmm all` — base orientation

## Purpose

`all` is the meta-command that chains the entire workflow:
extract → path-opt → tsopt → irc → freq → (optional) dft. The MEP stage
is single-pass `path-opt` by default; recursive `path-search` runs only
with `--refine-path`. It
resolves three input modes via flag context (see the three companion
mds: `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md`).

Use `all` when you want a **single qsub-able invocation** that
produces R / TS / P / IM coordinates plus barrier numbers for one or
more elementary steps.

## Synopsis

```bash
mlmm all [--parm enzyme.parm7] -i <input(s)> [-c <substrate>] [-l 'RES:Q,...'] \
    [--scan-lists '...'] [--tsopt True] [--thermo True] [--dft True] \
    [-b uma|orb|mace|aimnet2] [-o result_all/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology — optional; when omitted, `mm_parm` generates a parm7 from the input PDB |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | Force total system charge (highest priority over derived charges) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags (cross-mode)

| Flag | Type | Default | Description |
|---|---|---|---|
| `-i, --input` | path(s) | required | One or more reaction-ordered structures, or a TS-candidate alone |
| `-c, --center` | str | (uses input as-is) | Substrate selector: `'RES1,RES2,...'`, PDB path, or `'A:44,B:SAM'` |
| `-l, --ligand-charge` | str | none | Per-residue charges, e.g. `'SAM:1,GPP:-3'` |
| `-q, --charge` | int | derived from `-l` | Total cluster charge override |
| `-m, --multiplicity` | int | 1 | Spin multiplicity (2S+1) |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) when `-c` triggers extraction |
| `--scan-lists` | repeated | none | Staged distance scans (mode 2 — `all-scan-list.md`) |
| `--tsopt BOOL` | BOOL | `False` | Run TS optimization after the MEP stage (path-opt/path-search) |
| `--thermo BOOL` | BOOL | `False` | Run freq + thermochemistry |
| `--dft BOOL` | BOOL | `False` | Run DFT single point on R / TS / P |
| `--dft-func-basis` | str | `wb97m-v/def2-tzvpd` | DFT functional/basis (when `--dft True`) |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_all/` | Top-level output directory |
| `--config` | path | none | YAML config applied before CLI flags |
| `--show-config` / `--dry-run` | flag | off | Print resolved config and exit |
| `--help-advanced` | flag | — | Reveal advanced flags (extraction / `mm_parm` / freeze-atom / stage-specific overrides) |

Run `mlmm all --help-advanced` for the full list (it changes
between versions).

## Mode selection cheatsheet

```
Single -i input.{xyz,pdb,gjf} (no --scan-lists, no extra inputs)
    └── all-ts-only.md     (treat input as TS candidate; tsopt+irc+freq)

Single -i input.pdb + --scan-lists '...'
    └── all-scan-list.md   (single reactant + staged distance scans)

Multiple -i 1.R.pdb [2.IM.pdb ...] N.P.pdb (reaction-ordered)
    └── all-endpoint-mep.md (multi-endpoint MEP)
```

## Output tree (typical)

```
result_all/
├── summary.json                    # machine-readable per-stage results
├── summary.log                     # human-readable text + dir tree
├── mep.pdb / mep_trj.xyz           # full path (copied to the root)
├── mep_plot.png / energy_diagram_MEP.png
├── ml_region.pdb / mm_parm/ / layered/   # reusable ONIOM setup (--model-pdb / --parm inputs)
├── segments/
│   └── seg_NN/                     # canonical R/TS/P/IM + per-stage output
│       ├── reactant.pdb / .xyz
│       ├── ts.pdb / .xyz
│       ├── product.pdb / .xyz
│       ├── tsopt/                  # TS optimization output (--tsopt)
│       ├── irc/                    # forward/backward IRC trajectories
│       ├── freq/                   # frequencies + thermo (--thermo)
│       └── dft/                    # single-point DFT (--dft)
└── _work/                          # pipeline scratch (safe to delete)
    ├── pockets/ / scan/ / add_elem_info/
    └── path_opt/                   # raw MEP-engine output (path_search/ with --refine-path)
        ├── seg_NN_mep/ / hei_seg_NN.* / mep_trj.xyz / mep.pdb
        └── energy_diagram_*.png
```

`segments/seg_NN/` is the primary place to look for R/TS/P/IM
coordinates after a successful run; per-stage working files live in its
`tsopt/`, `irc/`, `freq/`, `dft/` subdirectories. See
`mlmm-workflows-output/SKILL.md` for canonical path conventions and the
bond-change interpretation.

## Output keys (summary.json — top level)

```python
import json
d = json.load(open("result_all/summary.json"))
print(d["status"])                    # "success" / "partial" / "failed"
print(d["mlmm_toolkit_version"])
print(d["charge"], d["spin"])
print(d["rate_limiting_step"])        # which segment is rate-limiting
print(len(d["segments"]))             # number of elementary steps
for seg in d["segments"]:
    print(seg["index"], seg["barrier_kcal"], seg["delta_kcal"])
```

Per-segment fields include `barrier_kcal`, `delta_kcal`, `bond_changes`,
`structures` (paths to reactant/ts/product files), `tsopt` /
`irc` / `freq` / `dft` sub-objects with their own `status` and energies.

## Resume / restart

`mlmm all` writes restart files when `--dump` is on (see
`--help-advanced`). To rerun only a failed segment:

```bash
mlmm tsopt -i _work/path_opt/hei_seg_03.xyz -o segments/seg_03/tsopt -b uma
mlmm irc   -i segments/seg_03/tsopt/final_geometry.xyz -o segments/seg_03/irc -b uma
mlmm freq  -i segments/seg_03/irc/finished_irc_trj.xyz -o segments/seg_03/freq -b uma
```

The directory layout matches what `all` produces, so downstream
analysis scripts keep working.

## Caveats

- `--scan-lists` is a Python literal-eval expression. Most
  shell-quoting trouble traces back to single vs double quotes.
- If `summary.json` shows `"status": "failed"` for any segment, look
  at the corresponding `summary.log` block; per-stage errors are also
  duplicated into `segments/seg_NN/<stage>/result.json`.
- The `segments/seg_NN/` deliverable directory is **only populated on success**
  for that segment. Failed segments leave `_work/path_opt/seg_NN_mep/`
  scratch (`_work/path_search/...` under `--refine-path`) but not the
  `segments/seg_NN/` copy.

## See also

- `all-endpoint-mep.md`, `all-scan-list.md`, `all-ts-only.md` — three
  invocation modes.
- `extract.md`, `path-search.md`, `tsopt.md`, `irc.md`, `freq.md`,
  `dft.md` — the underlying subcommands.
- `mlmm-workflows-output/SKILL.md` — output schema and
  R/TS/P coordinate conventions.
- Defaults: `import mlmm.core.defaults` (`OUT_DIR_ALL`, plus the
  per-stage `*_KW` dicts).