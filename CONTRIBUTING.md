# Contributing to mlmm-toolkit

Thank you for your interest in contributing to **mlmm-toolkit**.

This document is for **contributors and maintainers**. For end-user usage, see [`README.md`](README.md) and [`docs/getting-started.md`](docs/getting-started.md).

---

## 1. Before you start

Run the validation relevant to the files and numerical paths changed.

### 1.1 Required validation

| stage | what runs | how to invoke locally | failure means |
|---|---|---|---|
| 1. Unit tests | all configured pytest roots | `pytest -q` | logic regression |
| 2. Engineering markers | `# CHEMISTRY-RULE:N` coverage, `# DOMAIN_PURE` coverage, external-library import scope | `python .github/scripts/check_engineering_markers.py` | a required marker is missing, or an MLIP SDK is imported outside `backends/` |
| 3. Help registry drift | CLI `--help` and `--help-advanced` compliance with registry | `python .github/scripts/check_help_registry.py` | CLI option mismatch — re-run after CLI changes |
| 4. Smoke | `tests/smoke/run.sh` exercises the canonical ONIOM CLI surface (`mm-parm` → `define-layer` → `extract` → `path-search` → `tsopt` → `irc` → `freq` → `all`) on a representative system | copy `tests/smoke/` to scratch, then invoke `bash run.sh` from a site-specific scheduler wrapper | functional regression |

### 1.2 Change preparation

1. Read [`docs/architecture.md`](docs/architecture.md) §5 "Scientific invariants" before changing the VRAM-release sequence, chemistry rules, bundled forks, package discovery, dependencies, or `_LAZY_SUBCOMMANDS`.
2. Check [`mlmm/core/defaults.py`](mlmm/core/defaults.py) for shared runtime
   defaults and the relevant Click option definition for command-local CLI
   defaults.
3. Before editing `pysisyphus/`, `thermoanalysis/`, or `hessian_ff/`, read that directory's `README.md`. Validate logic changes with focused regression tests and the relevant scheduled numerical tests.
4. Identify which layer your change belongs to (`cli/`, `workflows/`, `domain/`, `backends/`, `io/`, `core/`). Stay inside one layer per commit when possible; the dependency direction is `L1 → L2 → {L3, L4} → L5` and must not be inverted.

### 1.3 Dev setup (lint / type-check tooling)

`ruff` and `pyright` are recommended local checks (run them before pushing;
they are not CI-enforced gates). Both are dev-only tools not installed by the
runtime `pip install`. Add them once:

```bash
pip install -e ".[dev]" ruff pyright
```

(or install ruff + pyright from your package manager of choice; pin via
`pip install ruff==0.6.* pyright==1.1.*` for reproducible optional local checks.)

### 1.4 Diagnostic dump examples

```bash
# Default run — INFO-level logging, no diagnostic dump
mlmm all -i R.pdb P.pdb -c 'SAM' -l 'SAM:1' -b uma --out-dir ./result_all

# --dump on freq: write thermoanalysis.yaml alongside the standard outputs
mlmm freq -i opt.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 --dump

# --dump-hess <path>: dump the Hessian for downstream IRC restart
mlmm freq -i opt.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
    --dump-hess /scratch/hess.npz
```

Use `--dump` when reproducing a resource-intensive regression or attaching artefacts to a bug report. Use `-v 3` when diagnosing an import-time or stage-bridge issue (e.g. AmberTools preflight failure, parm7 mismatch); the additional log volume is acceptable for short runs.

### 1.5 Downstream parser freeze rule

Treat `summary.json` as a versioned public schema and `summary.log` as
human-readable output. Document intentional schema changes, add a migration
test for structural breaks, and preserve machine-readable keys unless the
schema version changes.

---

## 2. Project layout

See [`docs/architecture.md`](docs/architecture.md) for the full 6-layer dir tree, file index, dependency direction, and recommended reading order. Short version:

- `mlmm/cli/` — L1 Interface (Click group, decorator factories, help, bool compat, subcommand resolver, AmberTools preflight).
- `mlmm/workflows/` — L2 Application (one file per subcommand stage runner, including the ONIOM-specific `define_layer`, `mm_parm`, `oniom_export`, `oniom_import`).
- `mlmm/domain/` — L3 Domain (chemistry-aware helpers: bond changes, bond summary, element info).
- `mlmm/backends/` — L4a Infra (MLIP dispatcher + ML/MM ONIOM calculator core; future splits into per-backend adapter + ONIOM subdir).
- `mlmm/io/` — L4b Infra (summary, trajectory, diagram, PDB fix, Hessian cache, analytical-Hessian glue).
- `mlmm/core/` — L5 Foundation (defaults, utils, future errors / types / logging).
- `pysisyphus/`, `thermoanalysis/`, `hessian_ff/` — bundled forks at the repo top; **not** upstream PyPI (and `hessian_ff/` has no upstream — bundling is mandatory).
- `tests/` — unit and regression tests.
- `tests/smoke/` — short representative job covering the canonical ONIOM CLI surface.

---

## 3. Recipes

The following "add-a-X" recipes name the files to touch and the checks that
cover common contributor changes.

### 3.1 Add a subcommand

**Goal**: expose a new CLI subcommand `mlmm myaction --opt1 X --opt2 Y`.

| step | action | file |
|---|---|---|
| 1 | Add a Python module `mlmm/workflows/myaction.py` with a top-level `@click.command(...)` named `cli` | new file in L2 |
| 2 | Register shared runtime defaults in `mlmm/core/defaults.py`; keep command-local Click defaults with the option definition | `mlmm/core/defaults.py` (L5), command module or `mlmm/cli/common_options.py` |
| 3 | Wire the command into the lazy registry — add `"myaction": ("mlmm.workflows.myaction", "cli", "<short description>")` to `_LAZY_SUBCOMMANDS` | `mlmm/cli/app.py` (L1) |
| 4 | Declare booleans as canonical Click toggle pairs (`--flag/--no-flag`). Runtime parameter introspection discovers ordinary decorators automatically; add a manual pre-import hint in `mlmm/cli/app.py` only when a lazy/parser-wrapper path requires it. | `mlmm/cli/default_group.py`, `mlmm/cli/app.py` |
| 5 | Add a docs page `docs/myaction.md` (and `docs/ja/myaction.md` if you maintain the JP set); add a unit test in `tests/test_myaction.py` | new files |

**Validation**: the unit test exercises the behavior, the reference check
detects registry/help drift, and the smoke suite covers commands on the
canonical smoke surface.

**Note on absolute paths**: `_LAZY_SUBCOMMANDS` entries MUST use **absolute** module paths (`mlmm.workflows.myaction`). Relative dotted strings (`".myaction"`) silently break subcommand discovery if `default_group.py` ever moves; see `docs/architecture.md` §5.5.

### 3.2 Add an MLIP backend

**Goal**: introduce a new MLIP backend `XYZModel` consumable as `--backend xyz`.

| step | action | file |
|---|---|---|
| 1 | Add the backend adapter and availability check to the existing ML backend abstraction | `mlmm/backends/mlmm_calc.py` |
| 2 | Add backend-specific model, precision, and accepted-option defaults | `mlmm/core/defaults.py` |
| 3 | Route the backend through calculator construction and provenance | `mlmm/backends/mlmm_calc.py`, `mlmm/core/utils.py` |
| 4 | Add the backend token to Click choices and generated help | `mlmm/cli/common_options.py`, relevant workflows |
| 5 | Document model identifiers, install command, accepted kwargs in `docs/backends.md`; add a smoke entry in `tests/smoke/run.sh` | `docs/backends.md`, `tests/smoke/run.sh` |

**Validation**: add focused factory/provenance tests, then exercise the backend
end to end in the smoke suite.

### 3.3 Add an output format

**Goal**: emit a new artefact `summary_v2.csv` alongside the existing `summary.json` / `summary.log`.

| step | action | file |
|---|---|---|
| 1 | Add a writer function in `mlmm/io/summary.py` that consumes the same in-memory summary dict | `mlmm/io/summary.py` (L4b) |
| 2 | Default emit path / on-or-off flag lives in `mlmm/core/defaults.py` | `mlmm/core/defaults.py` (L5) |
| 3 | Attach a per-subcommand `@click.option("--dump-<artefact>",...)` to the L2 stage runner when users can opt out | the stage runner |
| 4 | Advertise the new artefact in the output-layout documentation and the `summary.json` schema so downstream consumers can discover it | `docs/output-layout.md`, `mlmm/io/summary.py` |
| 5 | Add docs in `docs/json-output.md` + a unit test for round-trip serialisation | `docs/json-output.md`, new test |

**Validation**: the unit test in step 5 checks the writer; add smoke coverage
when the artefact belongs to the canonical smoke surface. Machine-readable
changes are governed by §1.5.

### 3.4 Add a workflow stage

**Goal**: insert a new stage (e.g. an intermediate `validate` step between TSOpt and Freq) into the `all` workflow.

| step | action | file |
|---|---|---|
| 1 | Implement the stage as a standalone subcommand first (Recipe 3.1) | `mlmm/workflows/validate.py` |
| 2 | Add an internal entry to the `all` pipeline orchestrator, preserving the VRAM `del` + `gc.collect()` pattern between stages | `mlmm/workflows/all.py` |
| 3 | Pass the stage result explicitly to the aggregate producer and any later consumer | `mlmm/workflows/all.py` |
| 4 | Update `mlmm/io/summary.py` to record the new stage's entry in `summary.json` | `mlmm/io/summary.py` |
| 5 | Update `tests/smoke/run.sh` to include the new stage in the representative run | `tests/smoke/run.sh` |

**Validation**: the smoke run in step 5 exercises a new `all` stage end to end before merge.

### 3.5 Add a test

**Goal**: add a unit test for new behaviour or a regression test for a fixed bug.

| step | action | file |
|---|---|---|
| 1 | Pick the right tier: pure-Python or chemistry-rule logic → `tests/test_<feature>.py`; multi-stage smoke → `tests/smoke/` | as appropriate |
| 2 | Use `pytest` style: one assertion per logical thing; name the test for the symptom (`test_irc_initial_displacement_does_not_oom`) | new test |
| 3 | If the test consumes a fixture, prefer the `tests/data/` directory; do **not** add large binary fixtures (> 100 KB) — use generators | `tests/data/`, `tests/conftest.py` |
| 4 | Run `pytest tests/test_<feature>.py -q -x` until green, then `pytest -q` to confirm no cross-test breakage across all configured roots | local |
| 5 | If the test depends on a new public Click command or symbol, land Recipe 3.1 / 3.3 first so the suite remains green | sequencing |

**Validation**: run `pytest`; CI blocks a failing merge.

---

## 4. Scientific and compatibility constraints

These are **hard constraints** enforced by the release process. Violating them either breaks correctness (chemistry rules), behaviour-level guarantees (`pyproject.toml` arrays, downstream-parser log lines), or upstream-fork compatibility.

### 4.1 Nine chemistry rules

The reaction-path correctness rules listed in [`docs/architecture.md`](docs/architecture.md) §5.1 must not be reordered, simplified, or factored out. They are marked with `# CHEMISTRY-RULE:N` inline comments and `# DOMAIN_PURE` module-docstring markers. The CI gate `.github/scripts/check_engineering_markers.py` enforces marker completeness and confines MLIP-only SDK imports (`fairchem`, `orb_models`, `mace`, `aimnet`) to the `backends/` layer. For **mlmm specifically** all 9 rules apply: #1 (subtractive ONIOM energy), #2 (link-atom Hessian B-matrix), #8 (3-layer 5-pass partial Hessian) in `backends/mlmm_calc.py`; #9 (parm7 atom indexing) in `io/pdb_indexing.py`; #3 (macro/micro alternation), #7 (`bofill_update` advanced-indexing) in `workflows/tsopt.py`; #6 (PHVA + MLIP active block) in `workflows/freq.py`; #4 (gpu4pyscf `rks_lowmem`), #5 (def2 auto-ECP) in `workflows/dft.py`.

To locate the markers:

```bash
grep -rnE '# CHEMISTRY-RULE:[0-9]+' mlmm/
grep -rn '# DOMAIN_PURE' mlmm/
```

### 4.2 VRAM-management invariant (`del` chains)

The IRC / TSopt / Freq stages explicitly `del calc`, `del geom`, `del hess` between stages and the `all` workflow runs `gc.collect()` at stage boundaries. **Do not refactor those `del` / `gc.collect()` statements out** — long-running ML/MM jobs with the full protein environment OOM without them.

### 4.3 Divergent files in bundled forks

The per-directory README tables are the live inventory of divergent files.
A logic change requires focused regression tests and the relevant scheduled
numerical tests.

Do not `pip install pysisyphus` or `pip install thermoanalysis` from PyPI alongside this package — silent runtime breakage. `hessian_ff/` has no upstream package; only the bundled copy works.

### 4.4 Package discovery and runtime dependencies

`pyproject.toml [tool.setuptools.packages.find].include` uses the `mlmm*` glob,
so files added under an existing layer need no discovery change. Check wheel
contents when introducing a new top-level package layout, and declare every
imported runtime package in `dependencies`.

### 4.5 `_LAZY_SUBCOMMANDS` absolute-path rule

Entries in `mlmm/cli/app.py:_LAZY_SUBCOMMANDS` MUST use absolute module paths (`"mlmm.workflows.all"`, never `".all"`). Relative dotted strings silently break the resolver when `default_group.py` moves; see [`docs/architecture.md`](docs/architecture.md) §5.5.

### 4.6 Chemistry default choices

Changes to the default basis set (def2-TZVPD), functional (ωB97M-V),
convergence thresholds, ECP handling, solvent models, or ONIOM region shell
radii require a documented numerical comparison and maintainer review. Inspect
`mlmm/core/defaults.py` and the command-local Click option before proposing a
change.

### 4.7 Downstream-parser-visible log lines

Any `summary.log` or `summary.json` line that downstream parsers consume is **frozen byte-for-byte**. See §1.5 above (Downstream parser freeze rule).

---

## 5. Where to ask

| forum | best for |
|---|---|
| GitHub issue | reproducible bug, feature request, design question with a concrete proposal |
| GitHub discussion | open-ended design / chemistry-method discussion, "what is the right way to..." |

---

## License

By contributing, you agree that your contributions will be licensed under the
[GPL-3.0 License](LICENSE).
