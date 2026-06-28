# `thermoanalysis/` (bundled fork)

> **This is a repo-internal fork of [thermoanalysis](https://github.com/eljost/thermoanalysis), NOT the upstream PyPI package. Do not `pip install thermoanalysis` alongside this package — it will silently overwrite the bundled copy.**

The fork is shipped inside this repository; treat it as part of `mlmm_toolkit` rather than a swappable upstream.

## Why a fork?

The bundled fork carries a small but important divergence from upstream that breaks side-by-side installation:

- **`QCData.py` branding / IO diff** — the bundled file emits output identifiers and quantum-chemistry data format expected by `mlmm/workflows/freq.py` (PHVA + UMA active-block partial Hessian consumer, chemistry-rule #6). Upstream's signatures do not match.

The rest of the package (`thermo.py`, `config.py`, `constants.py`) is close to upstream but rebuilt against the bundled `pysisyphus/` fork's unit conventions and the ONIOM partition-function conventions used by `mlmm`.

## Divergent files (do NOT replace with upstream)

| file | divergence | rule |
|------|------------|------|
| `QCData.py` | branding + IO signature differences | freq stage consumer contract, partial-Hessian path |

## release scope

During this release, **only annotation edits are allowed** on this directory:

- docstring additions / improvements
- type hints
- section banners (`# ===... ===`)
- per-file module docstring

**Forbidden** during polish:

- any change to numerical behaviour, control flow, or function signatures of `QCData.py`
- any new external dependency
- any rename of public symbols (would break `mlmm/workflows/freq.py` callers)

Logic edits to `QCData.py` would only be justified to track an upstream improvement; they must be explicitly requested via a `[CHEMISTRY-RULE]` commit and verified against the existing thermochemistry golden tests.

## Upstream compatibility

If you `pip install thermoanalysis` into the same environment as `mlmm_toolkit`, Python's import machinery may resolve to the bundled copy or the upstream copy depending on `sys.path` order. The flat-top placement at `<repo-top>/thermoanalysis/` is required for the bundled copy to take precedence in the editable-install path; do not move it under `mlmm/`.

## See also

- [`../docs/architecture.md`](../docs/architecture.md) §5.3, §6 — repo-internal fork policy
- [`../CONTRIBUTING.md`](../CONTRIBUTING.md) §4.3 — do-not-touch list (4 divergent fork files)
- `THIRD_PARTY_NOTICES.txt` — thermoanalysis upstream attribution
