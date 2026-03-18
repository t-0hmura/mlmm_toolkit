# Changelog

All notable changes to **mlmm-toolkit** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.2.2] — 2026-03-18

### Added
- GitHub Pages deployment and PyPI release workflows.
- Bidirectional scan (4-tuple) documentation (EN/JA).

### Fixed
- Regenerated CLI reference docs to match current `--help` output.
- Documentation accuracy: defaults, YAML examples, EN/JA alignment.
- Generalized UMA-specific wording to MLIP where applicable.
- Improved hessian_ff JIT build error message.

## [0.2.0] — 2026-03-16

**Complete rewrite from `mlmm` (v0.1.1) to `mlmm-toolkit` (v0.2.0).**
This release replaces the previous pysisyphus-wrapper architecture with a unified Click-based CLI toolkit. The package name, CLI interface, dependency stack, and internal architecture have all changed.

### Breaking Changes

- **Package name**: `mlmm-toolkit` on PyPI. Install with `pip install mlmm-toolkit`.
- **CLI completely redesigned**: The old entry points (`mlmm`, `def_ml_region`, `bond_scan`, `ts_search`, `energy_summary`, `trj2fig`, `add_elem_info`, `get_freeze_indices`, `xyz_geom2pdb`) are replaced by a single `mlmm <subcommand>` interface with 21 subcommands.
- **OpenMM removed**: MM calculations now use `hessian_ff`, a bundled C++ native extension for Amber force fields. OpenMM and OpenMM-CUDA-12 are no longer dependencies.
- **RDKit removed**: No longer required.
- **pysisyphus bundled**: No longer installed from a separate git repository; a modified fork is included in the package.
- **fairchem-core from PyPI**: No longer installed from a custom git fork.
- **numpy constraint relaxed**: `numpy<2.0` → `numpy>=1.24` (NumPy 2.x compatible).
- **Configuration**: Inline kwargs replaced by a centralized `defaults.py` as the single source of truth for all default values.

### Added — CLI & Workflow

- **`mlmm all`**: End-to-end workflow command. Given PDB files (R → P), automatically extracts active-site pockets, generates MM parameters, assigns ONIOM layers, runs MEP search, and optionally performs TS optimization, IRC, vibrational analysis, and single-point DFT — all in one invocation.
- **`mlmm opt`**: Single-structure geometry optimization with LBFGS (grad mode) or RFO (hess mode) with microiteration support.
- **`mlmm scan` / `scan2d` / `scan3d`**: 1D, 2D, and 3D constrained distance scans along user-specified atom pairs.
- **`mlmm path-search`**: Recursive minimum-energy path search with GSM (Growing String Method) and DMF (Direct Max Flux).
- **`mlmm path-opt`**: Single-pass MEP optimization (GSM or DMF).
- **`mlmm tsopt`**: Transition-state optimization (dimer method with partial Hessian, or RS-I-RFO with full Hessian).
- **`mlmm irc`**: Intrinsic reaction coordinate calculation from a TS geometry. Outputs `forward_last.pdb` and `backward_last.pdb` endpoint structures.
- **`mlmm freq`**: Vibrational analysis and thermochemistry (partial or full Hessian).
- **`mlmm dft`**: GPU-accelerated single-point DFT via PySCF / gpu4pyscf (`pip install "mlmm-toolkit[dft]"`).
- **`mlmm extract`**: Active-site pocket extraction from full protein-ligand PDB structures.
- **`mlmm mm-parm`**: Automatic Amber parm7/rst7 generation via AmberTools (tleap + GAFF2/AM1-BCC).
- **`mlmm define-layer`**: Assign 3-layer ML/MM partitioning via B-factor encoding (ML=0, MovableMM=10, FrozenMM=20).
- **`mlmm oniom-export` / `oniom-import`**: Gaussian/ORCA ONIOM input/output interoperability.
- **`mlmm add-elem-info`**: Fix missing or incorrect PDB element columns.
- **`mlmm fix-altloc`**: Resolve alternate location indicators in PDB files.
- **`mlmm energy-diagram`**: Plot energy profiles with Plotly (interactive HTML + static PNG).
- **`mlmm trj2fig`**: Trajectory visualization (PNG snapshots from XYZ trajectories).

### Added — Core Features

- **Multi-backend MLIP support**: UMA (default), ORB, MACE, AIMNet2 — selectable via `-b/--backend` on all subcommands.
- **ONIOM-like ML/MM decomposition**: `E = E_MM_real + E_ML_model − E_MM_model` with link-atom Jacobian transformation.
- **hessian_ff**: Bundled C++ native extension for analytical Amber force field energies, forces, and Hessians (replaces OpenMM).
- **Microiteration scheme**: Efficient optimization of large ML/MM systems (~10,000 atoms) by separating ML and MM degrees of freedom.
- **xTB embed-charge correction**: Optional point-charge embedding for the ML region (`--embedcharge`).
- **Partial Hessian approach**: Compute Hessians only for the ML region + boundary atoms, enabling TS optimization and frequency analysis on large systems.
- **Automatic versioning**: setuptools-scm replaces hard-coded `__version__`.
- **Progressive help**: `--help` shows primary options; `--help-advanced` shows the full option set.
- **DefaultGroup**: Lazy-loading Click subcommand architecture with automatic boolean normalization (`--flag/--no-flag` and `--flag True/False` both supported).

### Added — Dependencies & Extras

- **Optional backends**: `pip install "mlmm-toolkit[orb]"`, `"mlmm-toolkit[mace]"`, `"mlmm-toolkit[aimnet2]"`.
- **Optional DFT**: `pip install "mlmm-toolkit[dft]"` for PySCF + gpu4pyscf-cuda12x + CuPy.
- **Optional PDBFixer**: `pip install "mlmm-toolkit[pdbfixer]"` for hydrogen addition.
- **New core dependencies**: `click`, `torch_geometric`, `torch_scatter`, `pyparsing`, `tabulate`.
- **Bundled packages**: `pysisyphus` (modified fork), `thermoanalysis`, `hessian_ff`.

### Added — Documentation & Testing

- Bilingual documentation (English + Japanese) under `docs/` and `docs/ja/`.
- Per-subcommand reference docs with CLI tables and YAML configuration examples.
- `CONTRIBUTING.md` with contributor guidelines.
- GitHub Actions CI: `pytest.yml`, `smoke_test.yml`, `docs_quality.yml`.
- 198 unit tests with `pytest --timeout=120`.
- Smoke test suite: 34 end-to-end tests covering all subcommands.
- Working examples in `examples/toy_system/` and `examples/methyltransferase/`.

### Fixed (relative to v0.1.1)

- IRC energy-based initial displacement: Added `step_length` clamp (max 0.5 au) to prevent divergence when `min_eigval ≈ 0`.
- Hessian cache: Fixed `set_calculator()` clearing `within_partial_hessian` before `cart_hessian` assignment.
- HessianOptimizer: Fixed GPU/CPU device mismatch in `hessian_recalc`.
- Frequency analysis: Added float64 enforcement and explicit `(H+Hᵀ)/2` symmetrization.
- PDB element inference: Fixed `_infer_element_from_pdb_atom_name()` misidentifying protein hydrogens (e.g., `HG2` → `Hg`). Now uses residue-context-aware `guess_element()`.
- Silent `except Exception: pass` blocks converted to `logger.debug(...)` across all modules.
- Charge derivation: Uses `--model-pdb` instead of full input when provided.
- TypeError in `_pseudo_irc_and_match` when mapping value is `None`.

### Architecture — v0.1.1 → v0.2.0

| Aspect | v0.1.1 (`mlmm`) | v0.2.0 (`mlmm-toolkit`) |
|--------|------------------|--------------------------|
| Package name | `mlmm` | `mlmm-toolkit` |
| Codebase size | ~3,700 lines (13 files) | ~35,000 lines (40+ files) |
| CLI framework | 8 separate entry points (argparse) | 1 entry point, 21 Click subcommands |
| MM backend | OpenMM (finite difference) | hessian_ff (analytical C++) |
| ML backends | UMA, AIMNet2 | UMA, ORB, MACE, AIMNet2 |
| pysisyphus | Git dependency | Bundled (modified fork) |
| fairchem-core | Git fork | PyPI release |
| Configuration | Inline kwargs | Centralized `defaults.py` |
| ONIOM support | Basic link atoms | Full Jacobian + 3-layer B-factor |
| Documentation | README only | Full bilingual docs (EN/JA) |
| CI/CD | None | GitHub Actions (3 workflows) |
| Tests | None | 198 unit + 34 smoke tests |
