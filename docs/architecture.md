# Architecture: mlmm-toolkit

> Audience: contributors and maintainers (not end users — see [`getting-started.md`](getting-started.md) for usage).
> Scope: 6-layer source tree map, layered responsibilities, hidden constraints, bundled-fork policy, and a Fresh-eyes 5-step navigation for new developers.

---

## 1. Overview

`mlmm-toolkit` is a Python CLI that performs **ML/MM (ONIOM) enzymatic reaction-path analysis** on a complete protein environment. From a PDB plus a substrate name, it automatically generates the parm7 topology, encodes the ONIOM region split (ML / Movable MM / Frozen MM) into B-factor channels, runs full-system Hessian-based TS search via a macro/micro alternation scheme, and produces the reaction path (extract → MM-param → ONIOM model → MEP → tsopt → IRC → freq → dft).

The package is laid out as **6 physical layer directories** (`cli/`, `workflows/`, `domain/`, `backends/`, `io/`, `core/`); the role and dependency direction of each are summarised in the §4 layer tables below. External code imports directly from the layer directory (`from mlmm.backends.mlmm_calc import MLMMCore`, `from mlmm.core.utils import …`, `import mlmm.io.trj2fig`, etc.); the previous flat-top shim layer has been retired in this release.

Three bundled forks (`pysisyphus/`, `thermoanalysis/`, `hessian_ff/`) live at the repo top as repo-internal modules. They are deliberately **not** the upstream PyPI distributions (and `hessian_ff/` has no upstream at all — bundling is mandatory). See §6.

---

## 2. Layered structure (6 physical directories)

### 2.1 Layer table

| layer | dir | responsibility | may depend on |
|---|---|---|---|
| **L1 Interface** | `mlmm/cli/` | Click root group, decorator factories, `--help-advanced`, bool flag normalisation, subcommand resolver, AmberTools preflight | `workflows/`, `core/` |
| **L2 Application** | `mlmm/workflows/` | per-subcommand orchestration; one file per stage runner (`all.py`, `path_search.py`, `tsopt.py`, `extract.py`, `oniom_export.py`, `mm_parm.py`, …) | `domain/`, `backends/`, `io/`, `core/` |
| **L3 Domain** | `mlmm/domain/` | chemistry-aware helper logic (bond change detection, bond summary, element-info propagation) | `core/` |
| **L4a Infra (MLIP + ONIOM)** | `mlmm/backends/` | MLIP backend dispatcher + per-backend adapter + ML/MM ONIOM calculator core | `core/` |
| **L4b Infra (I/O)** | `mlmm/io/` | output layout, summary, trajectory, PDB fix, energy diagram, Hessian cache, analytical-Hessian glue | `core/` |
| **L5 Foundation** | `mlmm/core/` | defaults (single source of truth), utils (PDB / XYZ / plot helpers), future `errors.py` / `types.py` | (none) |
| (bundle, not a layer) | `<repo>/pysisyphus/`, `<repo>/thermoanalysis/`, `<repo>/hessian_ff/` | repo-internal forks (optimiser / thermochemistry / analytical MM Hessian) | (sibling, layer-external) |

**Dependency direction (one-way)**: `L1 → L2 → {L3, L4} → L5` (per the §2.1 layer table). The directional rule is enforced by CI marker coverage (`.github/scripts/check_engineering_markers.py`). Bundled forks sit outside the layer graph and may be imported from any layer through their absolute package path (`from pysisyphus.X import Y`, `from hessian_ff.analytical_hessian import …`).

### 2.2 ASCII map of the package tree

```
mlmm_toolkit/ [GH: t-0hmura/mlmm_toolkit]
├── pyproject.toml packages.find = ["mlmm*",...] (glob, frozen)
├── README.md / CONTRIBUTING.md / CHANGELOG.md
├── docs/
│ ├── architecture.md ← this file
│ └──... (Sphinx site, unchanged)
├── mlmm/ ← package body, 6-layer physical dir
│ ├── __init__.py PEP 562 lazy: _LAZY_IMPORTS + __getattr__
│ ├── __main__.py `from mlmm.cli.app import cli`
│ ├── _version.py / py.typed
│ │
│ ├── cli/ # === L1 Interface ===
│ │ ├── app.py Click group + _LAZY_SUBCOMMANDS registry (absolute paths)
│ │ ├── decorators.py @add_chem_common_options / @add_mm_layer_common_options et al.
│ │ ├── help_pages.py --help-advanced pager
│ │ ├── bool_compat.py --flag / --no-flag normalisation
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
│ │ ├── define_layer.py ML / Movable MM / Frozen MM B-factor assignment
│ │ ├── mm_parm.py AmberTools-driven parm7 / rst7 generation
│ │ ├── oniom_export.py ONIOM input writer (Gaussian / ORCA)
│ │ ├── oniom_import.py ONIOM input reader (sanity / atom-name diff)
│ │ └── align_freeze.py Kabsch + frozen-subset rmsd
│ │
│ ├── domain/ # === L3 Domain ===
│ │ ├── bond_changes.py R↔P bond detection
│ │ ├── bond_summary.py post-IRC diagnostic
│ │ └── add_elem_info.py PDB element column normaliser
│ │
│ ├── backends/ # === L4a Infra (MLIP + ONIOM) ===
│ │ ├── __init__.py backend dispatch + registry
│ │ ├── mlmm_calc.py ML/MM ONIOM calculator core (4 MLIP backends UMA / Orb / MACE / AIMNet2
│ │ inline; CHEMISTRY-RULE:1 / 2 / 8 / 9 host)
│ │ │ Future: split into base.py + per-backend uma.py / orb.py
│ │ │ / mace.py / aimnet2.py + ONIOM subdir
│ │ └── xtb_embedcharge_correction.py xTB point-charge embedding correction (--embedcharge)
│ │
│ ├── io/ # === L4b Infra (I/O) ===
│ │ ├── summary.py summary.json / summary.md writer
│ │ ├── energy_diagram.py Plotly diagram
│ │ ├── trj2fig.py trajectory → PNG / HTML / SVG / PDF
│ │ ├── pdb_fix.py altloc resolution
│ │ ├── hessian_cache.py in-memory Hessian cache
│ │ └── hessian_calc.py numerical-Hessian build + frequency / vibrational I/O helpers
│ │
│ └── core/ # === L5 Foundation ===
│   ├── defaults.py C1 single source of truth for every default
│   ├── utils.py PDB / XYZ / plot helpers
│   ├── logging.py -v / -vv logging wiring
│   ├── calc_eval.py per-stage calc evaluation
│   └── residue_data.py residue tables
│
├── tests/ smoke / unit
├── .github/ workflows/ + scripts/ (docs-quality lint helpers; CI-only)
└── (repo-top sibling, layer-external bundled forks)
 pysisyphus/ ~90 file, repo-internal fork (slimmed; CLI driver + QM backends + wavefunction + dead optimisers / IRC / NEB variants removed)
 thermoanalysis/ 5 file, repo-internal fork
 hessian_ff/ 19 file / 4.2k LOC, NO upstream PyPI, mandatory bundling
```

### 2.3 Per-layer responsibility detail

**L1 `cli/`**. Only this layer constructs Click commands and parses argv. `app.py` holds the root `Click.Group` plus the `_LAZY_SUBCOMMANDS` registry — every entry uses an **absolute module path** (`mlmm.workflows.all`, `mlmm.io.trj2fig`, …) so the resolver is independent of where `default_group.py` itself lives. The `mlmm`-specific `preflight.py` (AmberTools / conda env / GPU preflight) lives here because it runs during CLI startup before any L2 workflow is invoked.

**L2 `workflows/`** (~18 files). One file per subcommand. Each file owns a single `@click.command()` named `cli` and its private helpers. Large stage runners (`all.py` = 4,147 LOC, `path_search.py` = 2,348 LOC, `tsopt.py` = 3,068 LOC, `extract.py` = 2,321 LOC, `oniom_export.py` = 2,002 LOC) remain as single files in the current layout; future work may split them into per-stage subdirectories, but this is **opt-in** and out of scope for this release line.

**L3 `domain/`**. Chemistry-aware helper logic that may import `torch` / `numpy` / `pysisyphus.constants` (numeric back-ends), but **may not import** MLIP runtimes (`fairchem`, `orb_models`, `mace`, `aimnet`). The MLIP-runtime deny list (`fairchem` / `orb_models` / `mace` / `aimnet`) is enforced repo-wide by `.github/scripts/check_engineering_markers.py` (`_check_external_library_scope`), which forbids those imports in any module outside `backends/`. The separate `# DOMAIN_PURE` module-docstring marker is a distinct CI gate (`_check_domain_pure`) that flags the specific backend-agnostic modules required to stay MLIP-free — `backends/mlmm_calc.py`, `workflows/tsopt.py`, `workflows/freq.py` (and present on `workflows/sp.py`) — and is not itself the deny-list mechanism; no `domain/` file carries it. Domain helpers are reusable by any L2 stage runner.

**L4a `backends/`**. The ML/MM ONIOM calculator core (`mlmm_calc.py` = 2,534 LOC) lives here together with the backend dispatch (`__init__.py`) and the standalone xTB point-charge embedding correction (`xtb_embedcharge_correction.py`, driven by `--embedcharge`). Today the 4 MLIP backends (UMA / Orb / MACE / AIMNet2) that evaluate the ML region and the OpenMM / hessian_ff coupling all sit inline inside `mlmm_calc.py`; future work may split this into `backends/{base, uma, orb, mace, aimnet2}.py` for the MLIP layer plus a `backends/mlmm_calc/` subdir for the ONIOM core (`core.py`, `ase_calc.py`, `embed_charge.py`, `hessianff_calc.py`, `openmm_calc.py`, `facade.py`). The current single-file `mlmm_calc.py` carries chemistry rules **#1 (subtractive ONIOM)**, **#2 (link-atom Hessian B-matrix)**, **#8 (3-layer 5-pass partial Hessian)**, and **#9 (parm7 atom indexing)** — see §5.1.

**L4b `io/`** (7 files). Output-side I/O concerns: per-stage summary writer, energy diagram, trajectory rendering, PDB altloc fix, Hessian cache, numerical Hessian construction + frequency / vibrational I/O (`hessian_calc.py`), harmonic-restraint setup. `io/` never depends on `workflows/`; output format is owned here and consumed by stage runners.

**L5 `core/`**. The lowest layer. `defaults.py` is the **single source of truth** for every CLI default — grep here before adding a number anywhere else. `utils.py` is a ~2,560-LOC grab-bag of PDB / XYZ / plotting helpers; future work may split it into `utils/{pdb,plot,coord,yaml,freeze,input_prep}.py`. `logging.py` (`-v` / `-vv` wiring), `calc_eval.py` (per-stage calc evaluation) and `residue_data.py` (residue tables) also live here. The internal-only modules `errors.py`, `types.py` / `_stage.py` are introduced here as they land.

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

Two import surfaces (the flat-top shim layer was retired in this
release; downstream code that used `from mlmm.<oldmod>` must migrate
to the layered path):

1. **Layered import path**: external code imports directly from the layer directory (`from mlmm.backends.mlmm_calc import MLMMCore`, `from mlmm.core.utils import …`, `import mlmm.io.trj2fig`, etc.).
2. **Root symbol attribute** (`from mlmm import MLMMCore`) — handled by `mlmm/__init__.py:_LAZY_IMPORTS` + PEP 562 `__getattr__`. The five re-exported symbols (`MLMMCore`, `MLMMASECalculator`, `mlmm`, `mlmm_ase`, `mlmm_mm_only`) all resolve to `mlmm.backends.mlmm_calc` and are loaded on first access, so `import mlmm` stays cheap (only `__version__` is eager). There is **no** root module-attribute surface — submodules are reached by their full path (`import mlmm.io.trj2fig`), not as attributes of the top-level package.

The CLI subcommand resolver (`cli/app.py:_LAZY_SUBCOMMANDS`) uses **absolute** module paths (e.g. `"mlmm.workflows.all"`) so that moving `default_group.py` into `cli/` does not silently break subcommand discovery (the registry no longer depends on `__package__`).

---

## 3. Fresh-eyes 5-step navigation (≈ 40 min total)

For a contributor opening the repo for the first time, follow this path top-to-bottom; each step closes one concern.

| step | minutes | open | what you learn |
|------|---------|------|-----------------|
| 1 | 3 | [`README.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/README.md) | one-paragraph elevator pitch + single-command usage |
| 2 | 5 | this file (`docs/architecture.md`) §2 + §4 | 6-layer dir tree, dependency direction, where each concern lives |
| 3 | 5 | [`mlmm/cli/app.py`](../mlmm/cli/app.py) | Click root group, `_LAZY_SUBCOMMANDS` registry (≈ 18 entries), absolute-path resolution |
| 4 | 20 | [`mlmm/workflows/all.py`](../mlmm/workflows/all.py) (4,147 LOC, skim) | one full subcommand top-to-bottom; trace `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft` |
| 5 | 7 | [`CONTRIBUTING.md`](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md) §3 + §4 | 5 add-a-X recipes + the "do not touch" hidden constraints |

After step 5 you can read any other file by following the file index in §4. The package is intentionally **flat-within-each-layer** — there is no nested package below `mlmm/<layer>/` (other than the future split of `backends/mlmm_calc/` into per-backend modules), so you never need to navigate more than two directories deep.

---

## 4. File index — "where does this concern live?"

### 4.1 CLI / entry (L1 `cli/`)

| concern | file |
|---|---|
| Click root group + subcommand dispatch | `mlmm/cli/app.py` |
| Subcommand resolver (lazy import) | `mlmm/cli/default_group.py` |
| `python -m mlmm` shim | `mlmm/__main__.py` |
| Shared option decorator factories | `mlmm/cli/decorators.py` |
| `--help-advanced` pager | `mlmm/cli/help_pages.py` |
| Bool flag compat (`--flag` / `--no-flag` + value style) | `mlmm/cli/bool_compat.py` |
| AmberTools / conda env / GPU preflight | `mlmm/cli/preflight.py` |

### 4.2 Workflow stage runners (L2 `workflows/`)

| concern | file |
|---|---|
| Full pipeline orchestrator | `mlmm/workflows/all.py` |
| Geometry optimisation (ONIOM macro/micro pre-opt) | `mlmm/workflows/opt.py` |
| 1D / 2D / 3D scans + shared | `mlmm/workflows/scan{,2d,3d,_common}.py` |
| MEP search (GSM) | `mlmm/workflows/path_search.py` |
| MEP optimiser core (pysisyphus COS) | `mlmm/workflows/path_opt.py` |
| TS optimisation (RSIRFO + Bofill + macro/micro) | `mlmm/workflows/tsopt.py` |
| Vibrational analysis (PHVA + UMA active block) | `mlmm/workflows/freq.py` |
| IRC integration (macro / micro) | `mlmm/workflows/irc.py` |
| Single-point DFT (gpu4pyscf subprocess, ONIOM-embedded) | `mlmm/workflows/dft.py` |
| Active-site extraction (cluster cap) | `mlmm/workflows/extract.py` |
| ML / Movable MM / Frozen MM region assignment | `mlmm/workflows/define_layer.py` |
| AmberTools-driven MM parameter generation | `mlmm/workflows/mm_parm.py` |
| ONIOM input writer (Gaussian / ORCA) | `mlmm/workflows/oniom_export.py` |
| ONIOM input reader (sanity, atom-name diff) | `mlmm/workflows/oniom_import.py` |
| Kabsch / frozen-subset alignment | `mlmm/workflows/align_freeze.py` |

### 4.3 Chemistry helpers (L3 `domain/`)

| concern | file |
|---|---|
| R↔P bond change detection | `mlmm/domain/bond_changes.py` |
| Post-IRC bond summary | `mlmm/domain/bond_summary.py` |
| PDB element column normaliser | `mlmm/domain/add_elem_info.py` |

### 4.4 MLIP + ONIOM (L4a `backends/`)

| concern | file |
|---|---|
| ML/MM ONIOM calculator core + 4 inline MLIP backends + ONIOM coupling | `mlmm/backends/mlmm_calc.py` |
| Backend dispatch + registry | `mlmm/backends/__init__.py` |
| xTB point-charge embedding correction (`--embedcharge`) | `mlmm/backends/xtb_embedcharge_correction.py` |
| per-backend adapter split (planned, not yet present) | `mlmm/backends/{base, uma, orb, mace, aimnet2}.py` |
| ONIOM core subdir (planned, not yet present) | `mlmm/backends/mlmm_calc/{core, ase_calc, embed_charge, hessianff_calc, openmm_calc, facade}.py` |

See [MLIP Backends](backends.md) for the add-a-backend recipe (currently scoped to the planned per-backend split; until that lands, backend additions touch `mlmm_calc.py` inline).

### 4.5 I/O (L4b `io/`)

| concern | file |
|---|---|
| `summary.json` / `summary.md` writer | `mlmm/io/summary.py` |
| Plotly energy diagram | `mlmm/io/energy_diagram.py` |
| Trajectory → PNG / HTML / SVG / PDF | `mlmm/io/trj2fig.py` |
| PDB altloc resolution | `mlmm/io/pdb_fix.py` |
| In-memory Hessian cache (per-run TTL) | `mlmm/io/hessian_cache.py` |
| Numerical Hessian build + frequency / vibrational I/O | `mlmm/io/hessian_calc.py` |
| Harmonic restraint setup | `mlmm/workflows/restraints.py` (L2 stage helper) |

### 4.6 Foundation (L5 `core/`)

| concern | file |
|---|---|
| **Every CLI default (single source of truth)** | `mlmm/core/defaults.py` |
| PDB / XYZ / plot helpers | `mlmm/core/utils.py` |
| `-v` / `-vv` logging wiring | `mlmm/core/logging.py` |
| Per-stage calc evaluation | `mlmm/core/calc_eval.py` |
| Residue tables | `mlmm/core/residue_data.py` |
| (future) internal Protocol / TypedDict | `mlmm/core/types.py` |

### 4.7 Repo-internal bundled forks

| dir | role | divergent files (do NOT replace with upstream) |
|---|---|---|
| `pysisyphus/` | optimiser / TS / IRC engine | `irc/IRC.py`, `optimizers/hessian_updates.py`, `run.py`, `tsoptimizers/TSHessianOptimizer.py`, `calculators/*` (5 files total) |
| `thermoanalysis/` | thermochemistry (ΔG, ZPE, partition functions) | `QCData.py` (branding diff vs upstream) |
| `hessian_ff/` | analytical Hessian on MM force field — **PyPI 404, bundling is mandatory** | `analytical_hessian.py` (sole entry consumed by `mlmm/backends/mlmm_calc.py`) |

See each dir's `README.md` for the touch-restriction boundary.

---

## 5. Hidden constraints (read this before any patch)

### 5.1 Nine chemistry rules (grep recipe)

Nine correctness-critical rules are spread across `backends/`, `workflows/`, and `core/defaults.py`. They are **not** detected by smoke tests — silent drift here breaks reaction-path accuracy. Inline `# CHEMISTRY-RULE:N` markers and `# DOMAIN_PURE` module-docstring markers identify the rules; `.github/scripts/check_engineering_markers.py` enforces marker completeness in CI.

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
| 6 | PHVA + UMA active-block partial Hessian | `mlmm/workflows/freq.py` |
| 7 | `bofill_update` advanced-indexing scatter | `mlmm/workflows/tsopt.py` |
| 8 | 3-layer 5-pass partial Hessian assembly | `mlmm/backends/mlmm_calc.py` |
| 9 | parm7 atom indexing (1-based / serial gap handling) | `mlmm/backends/mlmm_calc.py` |

Editing any of these requires a `[CHEMISTRY-RULE:N]` commit prefix and a HEAVY-tier numerical-golden gate pass (see `CONTRIBUTING.md` §1.1).

**Recommended learning order (4 chemistry clusters)**:

| cluster | rules | shared concern | learn-first file |
|---|---|---|---|
| 5-pass Hessian set | #1, #2, #8, #9 | subtractive ONIOM + link-atom B-matrix + 3-layer assembly + parm7 indexing | `mlmm/backends/mlmm_calc.py` (host of 4 of the 9 rules) |
| TS optimisation set | #3, #7 | macro / micro alternation + Bofill scatter | `mlmm/workflows/tsopt.py` |
| Vibrational set | #6 | PHVA + UMA active-block partial Hessian | `mlmm/workflows/freq.py` |
| DFT set | #4, #5 | gpu4pyscf low-memory + def2 ECP injection | `mlmm/workflows/dft.py` |

For mlmm the practical curriculum is the 5-pass Hessian set first (#1, #2, #8, #9 — all in `mlmm_calc.py`), then the TS set (#3, #7), then DFT (#4, #5), then vibrational (#6).

### 5.2 VRAM-management invariant (do not refactor `del` chains)

The IRC / TSopt / Freq stages explicitly `del` GPU-resident objects (`calc`, `geom`, `hess`) between stages to free CUDA memory; the `all` workflow additionally runs `gc.collect()` at stage boundaries. **Do not refactor these `del` / `gc.collect()` statements out** — long-running ML/MM `all` jobs with the full protein environment OOM without them.

### 5.3 Bundled forks: do NOT install upstream alongside

The bundled `pysisyphus/`, `thermoanalysis/`, and `hessian_ff/` packages are **forks** (and in the case of `hessian_ff/`, the **only** available distribution — PyPI returns 404). Reinstalling `pip install pysisyphus` or `pip install thermoanalysis` next to this package silently breaks:

- `pysisyphus/irc/IRC.py` — initial-displacement memory hygiene
- `pysisyphus/optimizers/hessian_updates.py` — Bofill scatter on advanced indices, CPU-only `bofill_update` path for GPU OOM avoidance
- `pysisyphus/run.py` — stage wiring
- `pysisyphus/tsoptimizers/TSHessianOptimizer.py` — RSIRFO kwargs
- `pysisyphus/calculators/...` — GPU-aware backend hooks
- `thermoanalysis/QCData.py` — branding / I/O diff vs upstream
- `hessian_ff/analytical_hessian.py` — sole entry consumed by `backends/mlmm_calc.py`; **no upstream alternative exists**

### 5.4 `pyproject.toml` arrays are 0-diff

`[tool.setuptools.packages.find].include` and `dependencies` are treated as **0-diff arrays** during this release. The `include` glob (`mlmm*`) already auto-discovers any new layer subpackage; adding a `vendor/` or `internal/` container directory, or pinning a new runtime dependency, breaks the install contract and is forbidden by the release scope. Reflow / comment edits are fine; **array contents** are frozen.

### 5.5 `_LAZY_SUBCOMMANDS` registry must use absolute paths

`mlmm/cli/app.py:_LAZY_SUBCOMMANDS` resolves every subcommand through an **absolute** module path. Switching any entry back to a relative dotted import (`".all"` etc.) silently breaks subcommand discovery whenever `default_group.py` moves, because the resolver's `__package__` then drifts away from the package root. See internal design notes.

---

## 6. Bundled forks (repo-internal, two-repo independent evolution)

`mlmm_toolkit` ships **three** repo-internal modules at the repo top:

| dir | upstream PyPI? | purpose | scope of edits allowed |
|---|---|---|---|
| `pysisyphus/` | NO — fork, do not `pip install pysisyphus` alongside | optimiser, TS, IRC, COS, calculators | annotation-only in this release (docstring + type hints); logic edits forbidden |
| `thermoanalysis/` | NO — fork (branding diff) | ΔG, ZPE, partition functions, `QCData` | same as `pysisyphus/` |
| `hessian_ff/` | **NO — PyPI 404, bundling mandatory** | analytical Hessian on MM force field | same as `pysisyphus/` |

Each dir carries its own `README.md` listing the divergent files and the touch-restriction boundary. From the layer model these forks live **outside** the L1..L5 graph: any layer may import them via the absolute package path (`from pysisyphus.X import Y`, `from hessian_ff.analytical_hessian import …`) without breaking the `L1 → L2 → {L3, L4} → L5` direction.


---

## 7. Recommended deeper reading order (5–10 files)

After the Fresh-eyes tour (§3), follow this depth-first reading order:

1. `mlmm/core/defaults.py` — internalise the default-value table; everything downstream reads from here.
2. `mlmm/cli/app.py` — Click root + `_LAZY_SUBCOMMANDS` registry.
3. `mlmm/workflows/all.py` — one full pipeline top-to-bottom.
4. `mlmm/workflows/extract.py` + `define_layer.py` — cluster cap + ONIOM layer assignment.
5. `mlmm/workflows/mm_parm.py` — AmberTools parm7 generation.
6. `mlmm/backends/mlmm_calc.py` — the heart of ML/MM (chemistry-rules #1, #2, #8, #9 all live here).
7. `mlmm/workflows/tsopt.py` — RSIRFO + Bofill (CHEMISTRY-RULE:7) + macro / micro alternation (CHEMISTRY-RULE:3).
8. `mlmm/workflows/freq.py` — PHVA + UMA active-block (CHEMISTRY-RULE:6).
9. `mlmm/workflows/irc.py` — VRAM hygiene + macro / micro IRC.
10. `mlmm/core/utils.py` — shared PDB / XYZ / plot helpers.

---

## 8. ML/MM (ONIOM) scope

`mlmm-toolkit` operates on the **full protein environment** via ONIOM:

- **ML region**: substrate + reaction-centre residues, evaluated by one of 4 MLIP backends (UMA / Orb / MACE / AIMNet2); an optional xTB point-charge embedding correction (`--embedcharge`) adds MM→ML environmental effects
- **Movable MM region**: a shell around the ML region, free to move under the AMBER force field
- **Frozen MM region**: the rest of the protein, held rigid

The split is encoded in B-factor channels of the input PDB and propagated through `extract → mm-parm → ONIOM model → MEP → tsopt → IRC → freq → dft`.
