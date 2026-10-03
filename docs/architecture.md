# Architecture: mlmm-toolkit

## 1. Overview

This page is for people who change the mlmm-toolkit code: it describes the
package layers, where each file lives, and the constraints to respect before
you patch. After a patch, run the checks in
[Required validation](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md#11-required-validation)
of CONTRIBUTING. To run calculations instead, start from
[Getting Started](getting-started.md).

`mlmm-toolkit` is a Python CLI that performs **ML/MM (ONIOM) enzymatic reaction-path analysis** on a complete protein environment: a small reaction core is computed with a machine-learning interatomic potential (MLIP) and the surrounding protein with a molecular-mechanics (MM) force field, combined by subtractive ONIOM.

The `all` workflow chains `extract`, `mm-parm`, MEP search, TS optimization,
IRC, frequency analysis, and single-point DFT; the TS, thermochemistry, and DFT
stages are optional.

Three bundled forks, `pysisyphus/`, `thermoanalysis/`, and `hessian_ff/`, live at the repo top (§5.3, §6).

---

## 2. Layered structure (6 physical directories)

### 2.1 Layer table

| layer | dir | responsibility | may depend on |
|---|---|---|---|
| **L1 Interface** | `mlmm/cli/` | Click root group, decorator factories, `--help-advanced`, bool flag normalization, subcommand resolver, AmberTools preflight | `workflows/`, `core/` |
| **L2 Application** | `mlmm/workflows/` | per-subcommand orchestration plus shared workflow helpers (`_all_helpers.py`, `_opt_freq_common.py`, `_run_session.py`, …) | `domain/`, `backends/`, `io/`, `core/` |
| **L3 Domain** | `mlmm/domain/` | chemistry-aware helper logic (bond change detection, bond summary, element-info propagation) | `core/` |
| **L4a Infra (MLIP + ONIOM)** | `mlmm/backends/` | MLIP backend dispatch, inline backend integrations, and the ML/MM ONIOM calculator core | `core/` |
| **L4b Infra (I/O)** | `mlmm/io/` | output layout, summary, trajectory, PDB fix, energy diagram, Hessian cache, analytical-Hessian glue | `core/` |
| **L5 Foundation** | `mlmm/core/` | shared defaults, PDB/XYZ/plot helpers, result commit/output support, and residue tables | `backends/`, `domain/`, `io/` (a few upward imports, see below) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/`, `<repo>/hessian_ff/` | repo-internal forks (optimizer / thermochemistry / analytical MM Hessian) | (sibling, layer-external) |

**Dependency direction (design goal)**: `L1 → L2 → {L3, L4} → L5`. Shared charge/spin preparation and layer helpers live in `workflows/charge_prep.py` and `workflows/_opt_freq_common.py`. Bundled forks sit outside the layer graph and may be imported from any layer (`from pysisyphus.X import Y`).

CI checks only part of this direction:

- `.github/scripts/check_import_graph.py` forbids import cycles among `mlmm` modules, imports of `workflows` from `core` or `domain`, and imports of `mlmm` from the bundled forks.
- `.github/scripts/check_engineering_markers.py` checks the `# CHEMISTRY-RULE` and `# DOMAIN_PURE` markers (§5.1) and that MLIP runtimes are imported only under `backends/`.

### 2.2 ASCII map of the package tree

```
mlmm_toolkit/ [GH: t-0hmura/mlmm_toolkit]
├── pyproject.toml packages.find = ["mlmm*",...] (package-discovery glob)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx documentation site)
├── mlmm/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_IMPORTS + __getattr__
│ ├── __main__.py `from mlmm.cli.app import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── common_options.py @add_precision_option / @add_backend_model_option / @add_ml_charge_spin_options et al.
│ │ ├── decorators.py make_is_param_explicit, bool/YAML helpers, render_cli_exception
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalization
│ │ ├── default_group.py subcommand resolver, lazy module import
│ │ └── preflight.py AmberTools / conda env / GPU preflight
│ │
│ ├── workflows/ # === L2 Application ===
│ │ ├── all.py full pipeline orchestrator (extract → … → DFT)
│ │ ├── path_search.py / path_opt.py MEP search / COS wrapper
│ │ ├── tsopt.py / freq.py / irc.py / dft.py per-stage runners
│ │ ├── opt.py / scan.py / scan2d.py /
│ │ │ scan3d.py / scan_common.py ONIOM geometry opt / scans
│ │ ├── extract.py active-site extraction CLI
│ │ ├── define_layer.py ML / Movable-MM / Frozen B-factor assignment
│ │ ├── mm_parm.py AmberTools-driven parm7 / rst7 generation
│ │ ├── oniom_export.py ONIOM input writer (Gaussian / ORCA)
│ │ ├── oniom_import.py ONIOM input reader (sanity / atom-name diff)
│ │ ├── align_freeze.py Kabsch + frozen-subset rmsd
│ │ └── _all_helpers.py / _opt_freq_common.py / _run_session.py /
│ │     restraints.py shared workflow helpers
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ └── add_elem_info.py PDB element column normalizer
│ │
│ ├── backends/ # === L4a Infra (MLIP + ONIOM) ===
│ │ ├── __init__.py --precision routing (apply_precision_to_calc_cfg)
│ │ ├── mlmm_calc.py ML/MM ONIOM calculator core (4 MLIP backends UMA / ORB / MACE / AIMNet2
│ │ inline; CHEMISTRY-RULE:1 / 2 / 8 host)
│ │ ├── custom.py user ASE calculator loaded from --calc-file (custom backend)
│ │ ├── pyscf_dft.py optional PySCF/GPU4PySCF high-level adapter
│ │ └── _determinism.py strict-determinism setup (--deterministic)
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / HTML / SVG / PDF
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── pdb_indexing.py parm7 atom indexing (CHEMISTRY-RULE:9)
│ │ ├── hessian_cache.py in-memory Hessian cache
│ │ └── hessian_calc.py numerical-Hessian build + frequency / vibrational I/O helpers
│ │
│ ├── core/ # === L5 Foundation ===
│ │ ├── defaults.py shared workflow/calculator defaults
│ │ ├── dft_settings.py DFT settings (CHEMISTRY-RULE:4)
│ │ ├── utils.py PDB / XYZ / plot helpers
│ │ ├── logging.py -v/--verbose LEVEL (0–3) logging wiring
│ │ ├── calc_eval.py per-stage calc evaluation
│ │ ├── output.py / result_commit.py output/result commit helpers
│ │ ├── pes_composition.py energy-component composition
│ │ └── residue_data.py residue tables
│ │
│ └── mcp/ # non-layer subpackage: MCP server exposing every CLI subcommand
│   ├── server.py / _runner.py
│   └── _tools.py
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (CI, release, engineering, and docs checks)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ repo-internal optimizer, TS, IRC, COS, and calculator fork
 thermoanalysis/ repo-internal fork
 hessian_ff/ repo-internal native Hessian/MM support, NO upstream PyPI, mandatory bundling
```

### 2.3 Per-layer responsibility detail

**L1 `cli/`** owns root dispatch and shared argv parsing; the Click command of each subcommand is defined in its registered `workflows/`, `domain/`, or `io/` module. `app.py` holds the root `Click.Group` and the `_LAZY_SUBCOMMANDS` registry, whose entries use **absolute module paths** (§5.5). `preflight.py` (AmberTools / conda env / GPU preflight) lives here because it runs at CLI startup, before any L2 workflow.

**L2 `workflows/`** contains command modules plus shared workflow helpers. Modules registered in `cli/app.py:_LAZY_SUBCOMMANDS` own a `@click.command()` named `cli`; helper modules such as `_all_helpers.py`, `_opt_freq_common.py`, `_run_session.py`, `scan_common.py`, and `restraints.py` have no independent command.

**L3 `domain/`**. Chemistry-aware helper logic that may import `torch` / `numpy` / `pysisyphus.constants` (numeric back-ends), but **may not import** MLIP runtimes (`fairchem`, `orb_models`, `mace`, `aimnet`). Domain helpers are reusable by any L2 stage runner.

**L4a `backends/`**. The ML/MM ONIOM calculator core (`mlmm_calc.py`) lives here together with backend dispatch (`__init__.py`). The ML-region backends (UMA / ORB / MACE / AIMNet2) and OpenMM / hessian_ff coupling are dispatched from this layer. `mlmm_calc.py` hosts chemistry rules #1, #2, and #8 (§5.1).

**L4b `io/`**. Output-side I/O concerns include the per-stage summary writer, energy diagram, trajectory rendering, PDB/altloc handling, Hessian cache, numerical Hessian construction, and frequency/vibrational I/O (`hessian_calc.py`). `io/` never depends on `workflows/`; output format is owned here and consumed by stage runners.

**L5 `core/`**. The lowest layer. `defaults.py` is the **single source of truth** for shared defaults — grep here before adding a number elsewhere, then inspect justified command-local defaults. `utils.py` contains shared PDB / XYZ / plotting helpers.

### 2.4 Lazy-import mechanism (conceptual diagram)

```text
External consumer Package root Layer dir
------------------ ---------------- -----------

from mlmm.core.utils import x ────────────────────────────────────► mlmm/core/utils.py

import mlmm.io.trj2fig ──────────────────────────────────────────► mlmm/io/trj2fig.py

from mlmm.backends.mlmm_calc import ─────────────────────────────► mlmm/backends/mlmm_calc.py
 MLMMCore

from mlmm import MLMMCore ─────► mlmm/__init__.py
 __getattr__("MLMMCore")
 └─► _LAZY_IMPORTS["MLMMCore"]
 = "mlmm.backends.mlmm_calc"
 └─► importlib.import_module(...)
 └─► getattr(module, "MLMMCore")

mlmm myaction ─────────────────► mlmm/cli/app.py
 _LAZY_SUBCOMMANDS["myaction"]
 = ("mlmm.workflows.myaction", "cli", "...")
 └─► importlib.import_module(absolute path)
 └─► getattr(module, "cli") → Click command
```

Two import surfaces are supported:

1. **Layered import path**: external code imports directly from the layer directory, e.g. `from mlmm.backends.mlmm_calc import MLMMCore`.
2. **Root symbol attribute** (`from mlmm import MLMMCore`) — handled by `mlmm/__init__.py:_LAZY_IMPORTS` + PEP 562 `__getattr__`. The four re-exported symbols, `MLMMCore`, `MLMMASECalculator`, `mlmm`, and `mlmm_mm_only`, all resolve to `mlmm.backends.mlmm_calc` and are loaded on first access, so `import mlmm` stays cheap (only `__version__` is eager). Submodules are reached by their full path (`import mlmm.io.trj2fig`), not as attributes of the top-level package.

---

## 3. Fresh-eyes 5-step navigation (≈ 40 min total)

For a contributor opening the repo for the first time, follow this path top-to-bottom; each step closes one concern.

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md) | one-paragraph summary of the package + single-command usage |
| 2 | 5 | this file (`docs/architecture.md`) §2 + §4 | 6-layer dir tree, dependency direction, where each concern lives |
| 3 | 5 | [`mlmm/cli/app.py`](../mlmm/cli/app.py) | Click root group, `_LAZY_SUBCOMMANDS` registry (≈ 22 entries), absolute-path resolution |
| 4 | 20 | [`mlmm/workflows/all.py`](../mlmm/workflows/all.py) (skim) | one full subcommand top-to-bottom; trace `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) §3 + §4 | 5 add-a-X recipes + the "do not touch" hidden constraints |

After step 5 you can read any other file by following the file index in §4. The package is **flat within each layer**: there is no nested package below `mlmm/<layer>/`, so product modules are at most two directories below `mlmm/`.

---

## 4. File index — "where does this concern live?"

### 4.1 CLI / entry (L1 `cli/`)

| concern | file |
|---|---|
| Click root group + subcommand dispatch | `mlmm/cli/app.py` |
| Subcommand resolver (lazy import) | `mlmm/cli/default_group.py` |
| `python -m mlmm` shim | `mlmm/__main__.py` |
| Shared option decorator factories | `mlmm/cli/common_options.py` |
| Bool/YAML/exception CLI helpers | `mlmm/cli/decorators.py` |
| `--help-advanced` pager | `mlmm/cli/help_pages.py` |
| Bool flag parsing (`--flag` / `--no-flag` + value style) | `mlmm/cli/bool_compat.py` |
| AmberTools / conda env / GPU preflight | `mlmm/cli/preflight.py` |

### 4.2 Workflow stage runners (L2 `workflows/`)

Acronyms used below: GSM = growing-string method; COS = chain-of-states; RS-P-RFO = restricted-step partitioned rational-function optimization; RS-I-RFO = restricted-step image-function rational-function optimization; PHVA = partial Hessian vibrational analysis.

| concern | file |
|---|---|
| Full pipeline orchestrator | `mlmm/workflows/all.py` |
| Geometry optimization (ONIOM macro/micro pre-opt) | `mlmm/workflows/opt.py` |
| Scan and 2D/3D energy-landscape grids + shared | `mlmm/workflows/scan{,2d,3d,_common}.py` |
| MEP search (GSM) | `mlmm/workflows/path_search.py` |
| MEP optimizer core (pysisyphus COS) | `mlmm/workflows/path_opt.py` |
| TS optimization (RS-P-RFO / RS-I-RFO / TRIM + Bofill + macro/micro) | `mlmm/workflows/tsopt.py` |
| Vibrational analysis (PHVA + MLIP active block) | `mlmm/workflows/freq.py` |
| IRC integration (macro / micro) | `mlmm/workflows/irc.py` |
| Single-point DFT (ONIOM-embedded) | `mlmm/workflows/dft.py` |
| Active-site extraction (cluster cut-out + link-atom cap) | `mlmm/workflows/extract.py` |
| ML / Movable-MM / Frozen region assignment | `mlmm/workflows/define_layer.py` |
| AmberTools-driven MM parameter generation | `mlmm/workflows/mm_parm.py` |
| ONIOM input writer (Gaussian / ORCA) | `mlmm/workflows/oniom_export.py` |
| ONIOM input reader (sanity, atom-name diff) | `mlmm/workflows/oniom_import.py` |
| Kabsch / frozen-subset alignment | `mlmm/workflows/align_freeze.py` |

### 4.3 Chemistry helpers (L3 `domain/`)

| concern | file |
|---|---|
| R↔P bond change detection | `mlmm/domain/bond_changes.py` |
| Post-IRC bond summary | `mlmm/domain/bond_summary.py` |
| PDB element column normalizer | `mlmm/domain/add_elem_info.py` |

### 4.4 MLIP + ONIOM (L4a `backends/`)

| concern | file |
|---|---|
| ML/MM ONIOM calculator core + 4 inline MLIP backends + ONIOM coupling | `mlmm/backends/mlmm_calc.py` |
| `--precision` routing (`apply_precision_to_calc_cfg` / `_PRECISION_DISPATCH`) | `mlmm/backends/__init__.py` |
| Backend dispatch / factory (`_create_ml_backend`) | `mlmm/backends/mlmm_calc.py` |
| Optional DFT high-level adapter | `mlmm/backends/pyscf_dft.py` |

See [MLIP Backends](backends.md) for installation and runtime behavior. Backend
implementation changes currently touch `mlmm_calc.py` and the dispatcher.

### 4.5 I/O (L4b `io/`)

| concern | file |
|---|---|
| `summary.json` / `summary.log` writer | `mlmm/io/summary.py` |
| Plotly energy diagram | `mlmm/io/energy_diagram.py` |
| Trajectory → PNG / HTML / SVG / PDF | `mlmm/io/trj2fig.py` |
| PDB altloc resolution | `mlmm/io/pdb_fix.py` |
| parm7 atom indexing (CHEMISTRY-RULE:9) | `mlmm/io/pdb_indexing.py` |
| In-memory Hessian cache (per-run TTL) | `mlmm/io/hessian_cache.py` |
| Numerical Hessian build + frequency / vibrational I/O | `mlmm/io/hessian_calc.py` |
| Harmonic restraint setup | `mlmm/workflows/restraints.py` (L2 stage helper) |

### 4.6 Foundation (L5 `core/`)

| concern | file |
|---|---|
| Shared workflow and calculator defaults | `mlmm/core/defaults.py` |
| DFT settings (CHEMISTRY-RULE:4) | `mlmm/core/dft_settings.py` |
| PDB / XYZ / plot helpers | `mlmm/core/utils.py` |
| `-v/--verbose LEVEL` (0–3) logging wiring | `mlmm/core/logging.py` |
| Per-stage calc evaluation | `mlmm/core/calc_eval.py` |
| Output/result commit helpers | `mlmm/core/output.py`, `mlmm/core/result_commit.py` |
| Energy-component composition | `mlmm/core/pes_composition.py` |
| Residue tables | `mlmm/core/residue_data.py` |

### 4.7 Repo-internal bundled forks

| dir | role | divergent files (do NOT replace with upstream) |
|---|---|---|
| `pysisyphus/` | optimizer / TS / IRC engine | `irc/IRC.py`, `optimizers/hessian_updates.py`, `tsoptimizers/TSHessianOptimizer.py`, `calculators/*` |
| `thermoanalysis/` | thermochemistry (ΔG, ZPE, partition functions) | `QCData.py` (branding diff vs upstream) |
| `hessian_ff/` | analytical Hessian on MM force field — **PyPI 404, bundling is mandatory** | `analytical_hessian.py` (sole entry consumed by `mlmm/backends/mlmm_calc.py`) |

---

## 5. Scientific invariants

### 5.1 Nine chemistry rules (grep recipe)

Nine correctness-critical rules are spread across `backends/`, `workflows/`,
`core/`, and `io/`. Inline `# CHEMISTRY-RULE:N` markers identify their
implementation sites; `.github/scripts/check_engineering_markers.py` checks
marker completeness.

To find every chemistry rule before editing:

```bash
# List all 9 rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' mlmm/

# List every # DOMAIN_PURE marker (modules the CI check requires to carry it)
grep -rn '# DOMAIN_PURE' mlmm/
```

All 9 rules apply to `mlmm`:

| # | rule | host file |
|---|---|---|
| 1 | Subtractive ONIOM energy formula (`E = mm_real + ml_model − mm_model`) | `mlmm/backends/mlmm_calc.py` |
| 2 | Link-atom Hessian B-matrix projection | `mlmm/backends/mlmm_calc.py` |
| 3 | Macro / micro alternation for Hessian TS optimizers (RS-P-RFO default) | `mlmm/workflows/tsopt.py` |
| 4 | gpu4pyscf `rks_lowmem` closed-shell/GPU/lowmem guard | `mlmm/core/dft_settings.py` |
| 5 | def2 family auto-ECP injection | `mlmm/workflows/dft.py` |
| 6 | PHVA + MLIP active-block partial Hessian | `mlmm/workflows/freq.py` |
| 7 | `bofill_update` advanced-indexing scatter | `mlmm/workflows/tsopt.py` |
| 8 | 3-layer 5-pass partial Hessian assembly | `mlmm/backends/mlmm_calc.py` |
| 9 | parm7 atom indexing (1-based / serial gap handling) | `mlmm/io/pdb_indexing.py` |

Changes to these paths require focused regression tests and the relevant
scheduled numerical validation (see `CONTRIBUTING.md` §1.1).

**Recommended learning order (4 chemistry clusters)**:

| cluster | rules | shared concern | learn-first file |
|---|---|---|---|
| 5-pass Hessian set | #1, #2, #8, #9 | subtractive ONIOM + link-atom B-matrix + 3-layer assembly + parm7 indexing | `mlmm/backends/mlmm_calc.py` (host of 3 of the 9 rules: #1/#2/#8; #9 in `mlmm/io/pdb_indexing.py`) |
| TS optimization set | #3, #7 | macro / micro alternation + Bofill scatter | `mlmm/workflows/tsopt.py` |
| Vibrational set | #6 | PHVA + MLIP active-block partial Hessian | `mlmm/workflows/freq.py` |
| DFT set | #4, #5 | gpu4pyscf low-memory + def2 ECP injection | `mlmm/core/dft_settings.py` (#4), `mlmm/workflows/dft.py` (#5) |

For mlmm the practical curriculum is the 5-pass Hessian set first, then the TS set (#3, #7), then DFT (#4, #5), then vibrational (#6).

### 5.2 VRAM-management invariant (do not refactor `del` chains)

The IRC / TSopt / Freq stages explicitly `del` GPU-resident objects (`calc`, `geom`, `hess`) between stages to free CUDA memory; stage boundaries also use `gc.collect()` and, where CUDA allocations are present, `torch.cuda.empty_cache()`. **Do not refactor these release operations out** — long-running ML/MM `all` jobs on the full protein environment OOM without them.

### 5.3 Bundled forks: do NOT install upstream alongside

The bundled `pysisyphus/`, `thermoanalysis/`, and `hessian_ff/` packages are **forks**; `hessian_ff/` has no PyPI release at all. Reinstalling `pip install pysisyphus` or `pip install thermoanalysis` next to this package silently breaks:

- `pysisyphus/irc/IRC.py` — initial-displacement memory hygiene
- `pysisyphus/optimizers/hessian_updates.py` — GPU-resident in-place rank-two Bofill update; opt-in `PYSIS_BOFILL_CPU_OFFLOAD=1` fallback
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — Hessian TS optimizer kwargs
- `pysisyphus/calculators/...` — GPU-aware backend hooks
- `thermoanalysis/QCData.py` — branding / I/O diff vs upstream
- `hessian_ff/analytical_hessian.py` — sole entry consumed by `backends/mlmm_calc.py`; **no upstream alternative exists**

### 5.4 Package discovery and runtime dependencies

`[tool.setuptools.packages.find].include` uses the `mlmm*` glob to discover layer subpackages. A new top-level package layout must therefore be checked against wheel contents, and every imported runtime package must be declared in `dependencies`.

### 5.5 `_LAZY_SUBCOMMANDS` registry must use absolute paths

`mlmm/cli/app.py:_LAZY_SUBCOMMANDS` resolves every subcommand through an **absolute** module path. A relative dotted import (`".all"` etc.) would make resolution depend on the resolver module's `__package__` instead of the package root.

---

## 6. Bundled forks (repo-internal)

`mlmm_toolkit` ships **three** repo-internal modules at the repo top:

| dir | upstream PyPI? | purpose | scope of edits allowed |
|---|---|---|---|
| `pysisyphus/` | NO — fork, do not `pip install pysisyphus` alongside | optimizer, TS, IRC, COS, calculators | preserve the listed divergences; validate numerical changes with focused and scheduled numerical tests |
| `thermoanalysis/` | NO — fork (branding diff) | ΔG, ZPE, partition functions, `QCData` | preserve the `QCData` consumer contract; see its README |
| `hessian_ff/` | **NO — PyPI 404, bundling mandatory** | analytical Hessian on MM force field | preserve its derivative and public API contracts; see its README |

Each dir carries its own `README.md` listing the divergent files and the touch-restriction boundary.


---

## 7. Recommended deeper reading order (5–10 files)

After the Fresh-eyes 5-step navigation (§3), follow this depth-first reading order:

1. `mlmm/core/defaults.py` — internalise the default-value table; everything downstream reads from here.
2. `mlmm/cli/app.py` — Click root + `_LAZY_SUBCOMMANDS` registry.
3. `mlmm/workflows/all.py` — one full pipeline top-to-bottom.
4. `mlmm/workflows/extract.py` + `define_layer.py` — cluster cut-out + link-atom cap + ONIOM layer assignment.
5. `mlmm/workflows/mm_parm.py` — AmberTools parm7 generation.
6. `mlmm/backends/mlmm_calc.py` — the heart of ML/MM (CHEMISTRY-RULE:1, 2, 8).
7. `mlmm/workflows/tsopt.py` — Hessian TS optimizers + Bofill (CHEMISTRY-RULE:7) + macro / micro alternation (CHEMISTRY-RULE:3).
8. `mlmm/workflows/freq.py` — PHVA + MLIP active-block (CHEMISTRY-RULE:6).
9. `mlmm/workflows/irc.py` — VRAM hygiene + macro / micro IRC.
10. `mlmm/core/utils.py` — shared PDB / XYZ / plot helpers.

---

## 8. ML/MM (ONIOM) scope

`mlmm-toolkit` operates on the **full protein environment** via ONIOM:

- **ML region**: substrate + reaction-center residues, evaluated by one of the 4 MLIP backends (UMA / ORB / MACE / AIMNet2)
- **Movable-MM region**: a shell around the ML region, free to move under the AMBER force field
- **Frozen region**: the rest of the protein, held rigid

The split is encoded in B-factor channels of the input PDB and propagated through `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft`.

## Notes

- A few `core/` imports break the dependency direction today: `core.utils` imports `domain.add_elem_info` and `io.structure_formats`, and `core.calc_eval` imports `backends.mlmm_calc`; none of them forms a cycle.
- The `_check_domain_pure` gate only checks that the `# DOMAIN_PURE` marker is present on `backends/mlmm_calc.py`, `workflows/tsopt.py`, and `workflows/freq.py`; `workflows/sp.py` also carries it, and no `domain/` file does.
