# Changelog

All notable changes to **mlmm_toolkit** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- `-b/--backend` option added to all subcommands for consistent backend selection.
- `-o/--out-dir` option added to all subcommands for output directory control.
- `--embedcharge-cutoff` option for embed-charge distance threshold.

### Changed
- `--dry-run` moved to `--help-advanced` across all subcommands.
- `-s/--scan-lists` unified with `--spec` alias for scan specification.
- `--embedcharge` promoted to `--help` primary display.
- Bool options now documented as `--flag/--no-flag` alongside `--flag True/False`.

---

## [Previous]

### Added
- Multi-backend MLIP support: UMA, ORB, MACE, AIMNet2 via `--backend` option.
- `define-layer` subcommand for ML/MM layer assignment.
- `oniom-export` / `oniom-import` subcommands for Gaussian/ORCA ONIOM interoperability.
- `mm-parm` subcommand for Amber parm7/rst7 topology generation.
- Progressive help: `--help` (primary options) / `--help-advanced` (full options).
- `logging.getLogger(__name__)` in all CLI modules for structured logging.
- `CONTRIBUTING.md` and `CHANGELOG.md`.
- Bilingual documentation (English + Japanese) under `docs/` and `docs/ja/`.
- Coverage reporting (`pytest-cov`) with `--cov-fail-under=50` in CI.
- `--timeout=120` in pytest configuration.

### Fixed
- `freq.py`: Added float64 enforcement and explicit `(H+H^T)/2` symmetrization in `_mw_projected_hessian`.
- `irc.py`: Changed `return_partial_hessian` default from `True` to `False` for IRC calculations.
- `irc.py`: Consolidated duplicate `IRC_KW_DEFAULT` with `defaults.py` `IRC_KW`.
- Silent `except Exception: pass` blocks converted to `logger.debug(...)` across `oniom_export.py`, `utils.py`, `mm_parm.py`, `xtb_embedcharge_correction.py`, `opt.py`, `irc.py`.
- Documentation: YAML `bias.k` corrected from `100` to `300` in scan2d/scan3d docs.
- Documentation: `--opt-mode` aliases updated from `light`/`heavy` to `grad`/`hess` as primary names.
- Documentation: `--scan-lists` repeatable flag corrected in Japanese docs.
- Documentation: Missing `--convert-files` added to Japanese `dft.md` CLI table.
- Scoped pysisyphus log suppression (replaced global `logging.disable`).
- Version banner suppressed during `--help` tab completion (`ctx.resilient_parsing` guard).

### Changed
- `scan.py`: `--opt-mode` choices extended to include `grad`/`hess` for consistency with other subcommands.
