# Architecture: mlmm-toolkit

## 1. Overview

`mlmm-toolkit` is a Python CLI that performs **ML/MM (ONIOM) enzymatic reaction-path analysis** on a complete protein environment. ML/MM here means a hybrid model in which a small reaction core is treated by a machine-learning interatomic potential (ML) and the surrounding protein by a molecular-mechanics (MM) force field, combined through the subtractive ONIOM (Our own N-layered Integrated molecular Orbital and molecular Mechanics) energy scheme.

The `all` workflow can generate a parm7 topology and encode the ONIOM region
split (ML / Movable-MM / Frozen) from structure and region inputs. Multi-input
runs search a path; single-input runs require either `--scan-lists` or
`--tsopt`.

Depending on its mode and flags, `all` composes stages from `extract`,
`mm-parm`, MEP search, TS optimization, IRC, frequency analysis, and
single-point DFT. TS, thermochemistry, and DFT stages are optional.

The package is laid out as **6 physical layer directories** (`cli/`, `workflows/`, `domain/`, `backends/`, `io/`, `core/`). The role and dependency direction of each are summarized in the §2.1 layer table below.

External code imports directly from the layer directory (`from mlmm.backends.mlmm_calc import MLMMCore`, `from mlmm.core.utils import …`, `import mlmm.io.trj2fig`, etc.). §2.4 details the two supported import surfaces.

Three bundled forks (`pysisyphus/`, `thermoanalysis/`, `hessian_ff/`) live at the repo top as repo-internal modules. They are deliberately **not** the upstream PyPI distributions (and `hessian_ff/` has no upstream at all — bundling is mandatory). See §6.

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
| **L5 Foundation** | `mlmm/core/` | shared defaults, PDB/XYZ/plot helpers, result commit/output support, and residue tables | `backends/`, `domain/`, `io/` (compatibility back-edges, see below) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/`, `<repo>/hessian_ff/` | repo-internal forks (optimizer / thermochemistry / analytical MM Hessian) | (sibling, layer-external) |

**Dependency direction**: the *design intent* is one-way `L1 → L2 → {L3, L4} → L5`. What is actually **enforced** — and by which checker — is:

- `.github/scripts/check_import_graph.py` (AST import-graph gate) enforces the parts that are true today: (a) the `mlmm` package has **no import cycle** (no strongly connected component among its modules); (b) `core` and `domain` **never import `workflows`**; (c) no bundled fork (`pysisyphus` / `hessian_ff` / `thermoanalysis`) imports `mlmm`.
- `.github/scripts/check_engineering_markers.py` covers a *different* concern — chemistry-rule / `DOMAIN_PURE` marker completeness and the MLIP-runtime import scope. It does **not** parse the layer import edges, so it does not enforce the direction.

The current cycle-free back-edges are limited to a few `core.utils` helpers importing `domain.add_elem_info` and `io.structure_formats`, plus `core.calc_eval` importing `backends.mlmm_calc`. `core.utils` does not import `workflows`. Shared charge/spin preparation and layer helpers live in `workflows/charge_prep.py` and `workflows/_opt_freq_common.py`. Bundled forks sit outside the layer graph and may be imported from any layer through their absolute package path (`from pysisyphus.X import Y`, `from hessian_ff.analytical_hessian import …`).

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
│ │ inline; CHEMISTRY-RULE:1 / 2 / 8 / 9 host)
│ │ ├── custom.py user ASE calculator loaded from --calc-file (custom backend)
│ │ ├── _determinism.py strict-determinism setup (--deterministic)
│ │ └── xtb_embedcharge_correction.py dormant compatibility implementation
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.log writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / HTML / SVG / PDF
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── hessian_cache.py in-memory Hessian cache
│ │ └── hessian_calc.py numerical-Hessian build + frequency / vibrational I/O helpers
│ │
│ ├── core/ # === L5 Foundation ===
│ │ ├── defaults.py shared workflow/calculator defaults
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

**L1 `cli/`**. This layer owns root dispatch and shared argv parsing; leaf Click commands are defined in the registered `workflows/`, `domain/`, and `io/` modules. `app.py` holds the root `Click.Group` plus the `_LAZY_SUBCOMMANDS` registry — every entry uses an **absolute module path** (`mlmm.workflows.all`, `mlmm.io.trj2fig`, …) so the resolver is independent of where `default_group.py` itself lives. The `mlmm`-specific `preflight.py` (AmberTools / conda env / GPU preflight) lives here because it runs during CLI startup before any L2 workflow is invoked.

**L2 `workflows/`** contains command modules plus shared workflow helpers. Modules registered in `cli/app.py:_LAZY_SUBCOMMANDS` own a `@click.command()` named `cli`; helper modules such as `_all_helpers.py`, `_opt_freq_common.py`, `_run_session.py`, `scan_common.py`, and `restraints.py` have no independent command. Large stage runners remain single modules in the current layout.

**L3 `domain/`**. Chemistry-aware helper logic that may import `torch` / `numpy` / `pysisyphus.constants` (numeric back-ends), but **may not import** machine-learning interatomic potential (MLIP) runtimes (`fairchem`, `orb_models`, `mace`, `aimnet`). Two distinct CI gates cover this, both in `.github/scripts/check_engineering_markers.py`:

- The MLIP-runtime deny list (`fairchem` / `orb_models` / `mace` / `aimnet`) is enforced repo-wide by `_check_external_library_scope`, which forbids those imports in any module outside `backends/`.
- The separate `# DOMAIN_PURE` module-docstring marker is a distinct CI gate (`_check_domain_pure`) that flags the specific backend-agnostic modules required to stay MLIP-free — `backends/mlmm_calc.py`, `workflows/tsopt.py`, `workflows/freq.py` (and present on `workflows/sp.py`). It is not itself the deny-list mechanism, and no `domain/` file carries it.

Domain helpers are reusable by any L2 stage runner.

**L4a `backends/`**. The ML/MM ONIOM calculator core (`mlmm_calc.py`) lives here together with backend dispatch (`__init__.py`) and a dormant electronic-embedding compatibility module (`xtb_embedcharge_correction.py`; public activation is rejected in v0.3.3). The ML-region backends (UMA / ORB / MACE / AIMNet2) and OpenMM / hessian_ff coupling are dispatched from this layer. `mlmm_calc.py` carries chemistry rules **#1 (subtractive ONIOM)**, **#2 (link-atom Hessian B-matrix)**, and **#8 (3-layer 5-pass partial Hessian)**; rule **#9 (parm7 atom indexing)** lives in `io/pdb_indexing.py` — see §5.1.

**L4b `io/`**. Output-side I/O concerns include the per-stage summary writer, energy diagram, trajectory rendering, PDB/altloc handling, Hessian cache, numerical Hessian construction, and frequency/vibrational I/O (`hessian_calc.py`). `io/` never depends on `workflows/`; output format is owned here and consumed by stage runners.

**L5 `core/`**. The lowest layer. `defaults.py` is the **single source of truth** for shared defaults — grep here before adding a number elsewhere, then inspect justified command-local defaults. `utils.py` contains shared PDB / XYZ / plotting helpers; `logging.py` (per-subcommand `-v/--verbose LEVEL`, 0–3), `calc_eval.py` (per-stage calculation evaluation), and `residue_data.py` (residue tables) also live here.

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

1. **Layered import path**: external code imports directly from the layer directory (see the §2.1 layer table; e.g. `from mlmm.backends.mlmm_calc import MLMMCore`).
2. **Root symbol attribute** (`from mlmm import MLMMCore`) — handled by `mlmm/__init__.py:_LAZY_IMPORTS` + PEP 562 `__getattr__`. The five re-exported symbols (`MLMMCore`, `MLMMASECalculator`, `mlmm`, `mlmm_ase`, `mlmm_mm_only`) all resolve to `mlmm.backends.mlmm_calc` and are loaded on first access, so `import mlmm` stays cheap (only `__version__` is eager). There is **no** root module-attribute surface — submodules are reached by their full path (`import mlmm.io.trj2fig`), not as attributes of the top-level package.

The CLI subcommand resolver (`cli/app.py:_LAZY_SUBCOMMANDS`) uses **absolute** module paths (e.g. `"mlmm.workflows.all"`) so subcommand discovery is independent of the resolver module's `__package__`.

---

## 3. Five-step navigation

For a contributor opening the repo for the first time, follow this path top-to-bottom; each step closes one concern.

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md) | one-paragraph elevator pitch + single-command usage |
| 2 | 5 | this file (`docs/architecture.md`) §2 + §4 | 6-layer dir tree, dependency direction, where each concern lives |
| 3 | 5 | [`mlmm/cli/app.py`](../mlmm/cli/app.py) | Click root group, `_LAZY_SUBCOMMANDS` registry (≈ 22 entries), absolute-path resolution |
| 4 | 20 | [`mlmm/workflows/all.py`](../mlmm/workflows/all.py) (skim) | one full subcommand top-to-bottom; trace `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) §3 + §4 | 5 add-a-X recipes + the "do not touch" hidden constraints |

After step 5 you can read any other file by following the file index in §4. The package is intentionally **flat within each layer**: there is no nested package below `mlmm/<layer>/`, so product modules are at most two directories below `mlmm/`.

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
| Bool flag compat (`--flag` / `--no-flag` + value style) | `mlmm/cli/bool_compat.py` |
| AmberTools / conda env / GPU preflight | `mlmm/cli/preflight.py` |

### 4.2 Workflow stage runners (L2 `workflows/`)

Acronyms used below: MEP = minimum-energy path; GSM = growing-string method; COS = chain-of-states; RSIRFO = restricted-step image-function rational-function optimization (also written RS-I-RFO); Bofill = the Bofill Hessian-update formula; PHVA = partial Hessian vibrational analysis; IRC = intrinsic reaction coordinate; Kabsch = the Kabsch rigid-body alignment algorithm.

| concern | file |
|---|---|
| Full pipeline orchestrator | `mlmm/workflows/all.py` |
| Geometry optimization (ONIOM macro/micro pre-opt) | `mlmm/workflows/opt.py` |
| 1D / 2D / 3D scans + shared | `mlmm/workflows/scan{,2d,3d,_common}.py` |
| MEP search (GSM) | `mlmm/workflows/path_search.py` |
| MEP optimizer core (pysisyphus COS) | `mlmm/workflows/path_opt.py` |
| TS optimization (RSIRFO + Bofill + macro/micro) | `mlmm/workflows/tsopt.py` |
| Vibrational analysis (PHVA + MLIP active block) | `mlmm/workflows/freq.py` |
| IRC integration (macro / micro) | `mlmm/workflows/irc.py` |
| Single-point DFT (gpu4pyscf subprocess, ONIOM-embedded) | `mlmm/workflows/dft.py` |
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
| Retired electronic-embedding compatibility module | `mlmm/backends/xtb_embedcharge_correction.py` |

See [MLIP Backends](backends.md) for installation and runtime behavior. Backend
implementation changes currently touch `mlmm_calc.py` and the dispatcher.

### 4.5 I/O (L4b `io/`)

| concern | file |
|---|---|
| `summary.json` / `summary.log` writer | `mlmm/io/summary.py` |
| Plotly energy diagram | `mlmm/io/energy_diagram.py` |
| Trajectory → PNG / HTML / SVG / PDF | `mlmm/io/trj2fig.py` |
| PDB altloc resolution | `mlmm/io/pdb_fix.py` |
| In-memory Hessian cache (per-run TTL) | `mlmm/io/hessian_cache.py` |
| Numerical Hessian build + frequency / vibrational I/O | `mlmm/io/hessian_calc.py` |
| Harmonic restraint setup | `mlmm/workflows/restraints.py` (L2 stage helper) |

### 4.6 Foundation (L5 `core/`)

| concern | file |
|---|---|
| Shared workflow and calculator defaults | `mlmm/core/defaults.py` |
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

See each dir's `README.md` for the touch-restriction boundary.

---

## 5. Scientific invariants

### 5.1 Nine chemistry rules (grep recipe)

Nine correctness-critical rules are spread across `backends/`, `workflows`,
and `core/defaults.py`. Inline `# CHEMISTRY-RULE:N` markers and
`# DOMAIN_PURE` module-docstring markers identify their implementation sites;
`.github/scripts/check_engineering_markers.py` checks marker completeness.

To find every chemistry rule before editing:

```bash
# List all 9 rule sites in the repo (host file + line)
grep -rnE '# CHEMISTRY-RULE:[0-9]+' mlmm/

# List every # DOMAIN_PURE marker (= chemistry-rule host modules)
grep -rn '# DOMAIN_PURE' mlmm/
```

All 9 rules apply to `mlmm`:

| # | rule | host file |
|---|---|---|
| 1 | Subtractive ONIOM energy formula (`E = mm_real + ml_model − mm_model`) | `mlmm/backends/mlmm_calc.py` |
| 2 | Link-atom Hessian B-matrix projection | `mlmm/backends/mlmm_calc.py` |
| 3 | Macro / micro alternation (RS-I-RFO hess mode microiteration) | `mlmm/workflows/tsopt.py` |
| 4 | gpu4pyscf `rks_lowmem` triple-guard | `mlmm/workflows/dft.py` |
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
| DFT set | #4, #5 | gpu4pyscf low-memory + def2 ECP injection | `mlmm/workflows/dft.py` |

For mlmm the practical curriculum is the 5-pass Hessian set first (#1, #2, #8 in `mlmm_calc.py`; #9 in `io/pdb_indexing.py`), then the TS set (#3, #7), then DFT (#4, #5), then vibrational (#6).

### 5.2 VRAM-management invariant (do not refactor `del` chains)

The IRC / TSopt / Freq stages explicitly `del` GPU-resident objects (`calc`, `geom`, `hess`) between stages to free CUDA memory; stage boundaries also use `gc.collect()` and, where CUDA allocations are present, `torch.cuda.empty_cache()`. **Do not refactor these release operations out** — long-running ML/MM `all` jobs on the full protein environment OOM without them.

### 5.3 Bundled forks: do NOT install upstream alongside

The bundled `pysisyphus/`, `thermoanalysis/`, and `hessian_ff/` packages are **forks** (and in the case of `hessian_ff/`, the **only** available distribution — PyPI returns 404). Reinstalling `pip install pysisyphus` or `pip install thermoanalysis` next to this package silently breaks:

- `pysisyphus/irc/IRC.py` — initial-displacement memory hygiene
- `pysisyphus/optimizers/hessian_updates.py` — GPU-resident in-place rank-two Bofill update; opt-in `PYSIS_BOFILL_CPU_OFFLOAD=1` fallback
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs
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

Each dir carries its own `README.md` listing the divergent files and the touch-restriction boundary. From the layer model these forks live **outside** the L1..L5 graph: any layer may import them via the absolute package path (`from pysisyphus.X import Y`, `from hessian_ff.analytical_hessian import …`) without breaking the `L1 → L2 → {L3, L4} → L5` direction.


---

## 7. Recommended deeper reading order (5–10 files)

After the Fresh-eyes tour (§3), follow this depth-first reading order:

1. `mlmm/core/defaults.py` — internalise the default-value table; everything downstream reads from here.
2. `mlmm/cli/app.py` — Click root + `_LAZY_SUBCOMMANDS` registry.
3. `mlmm/workflows/all.py` — one full pipeline top-to-bottom.
4. `mlmm/workflows/extract.py` + `define_layer.py` — cluster cut-out + link-atom cap + ONIOM layer assignment.
5. `mlmm/workflows/mm_parm.py` — AmberTools parm7 generation.
6. `mlmm/backends/mlmm_calc.py` — the heart of ML/MM (chemistry-rules #1, #2, #8 live here; #9 in `mlmm/io/pdb_indexing.py`).
7. `mlmm/workflows/tsopt.py` — RSIRFO + Bofill (CHEMISTRY-RULE:7) + macro / micro alternation (CHEMISTRY-RULE:3).
8. `mlmm/workflows/freq.py` — PHVA + MLIP active-block (CHEMISTRY-RULE:6).
9. `mlmm/workflows/irc.py` — VRAM hygiene + macro / micro IRC.
10. `mlmm/core/utils.py` — shared PDB / XYZ / plot helpers.

---

## 8. ML/MM (ONIOM) scope

`mlmm-toolkit` operates on the **full protein environment** via ONIOM:

- **ML region**: substrate + reaction-center residues, evaluated by one of 4 machine-learning interatomic potential (MLIP) backends (UMA / ORB / MACE / AIMNet2); v0.3.3 uses mechanical embedding
- **Movable-MM region**: a shell around the ML region, free to move under the AMBER force field
- **Frozen region**: the rest of the protein, held rigid

The split is encoded in B-factor channels of the input PDB and propagated through `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft`.
