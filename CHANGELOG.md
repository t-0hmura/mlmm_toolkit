# Changelog

All notable changes to **mlmm-toolkit** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## Unreleased

## [0.3.0] — 2026-06-15

### Changed
- `--dft-func-basis` is now surfaced in the primary `mlmm <subcmd> --help`
  (previously only under `--help-advanced`), so the DFT//MLIP functional/basis
  is discoverable without the advanced listing.
- Standalone `path-opt` / `path-search` now default to `--preopt` (each MEP
  endpoint is pre-optimized before the search); the previous default was
  `--no-preopt`. The `all` pipeline forwards the flag explicitly and is
  unaffected.
- `--precision fp64` now also forces the Hessian to fp64 (`H_double`) so the
  optimiser / eigen linear algebra cannot silently run in a lower precision
  than the model; a config that set `H_double=False` under fp64 is overridden
  with a warning.
- AIMNet2 now rejects both `--precision fp64` (its model inputs are cast to
  float32 upstream) and `--deterministic` (its forces come from a custom CUDA
  kernel outside torch's deterministic-algorithms control), with clear errors
  instead of running misleadingly.
- The `all`-pipeline determinism check in the smoke suite is now an
  informational monitor (reports drift, does not gate); bit-exactness is an
  opt-in via `--deterministic`, not a default guarantee.
- `--help` of `mlmm` now groups subcommands under semantic sections
  ("Pipelines" / "Pipeline stages" / "Inputs & topology" / "Analysis")
  in a configurable, deterministic order; subcommands not listed in any
  section fall through to a trailing "Other" bucket so we never hide an
  entry silently.
- CLI exception renderer appends `Try 'mlmm <subcmd> -h' for help.` to
  every user-input-style error so first-time users see a recovery path,
  and routes the full traceback through `logging.getLogger(...).exception`
  so log scrapers / `-v` users get the structured record alongside the
  human-readable terminal echo.
- `_calc_energy` deduplicated into `mlmm.core.calc_eval`; both
  `workflows/opt.py` and `workflows/tsopt.py` now re-export the helper.
- Repo-wide ruff `F401` sweep: removed 66 unused imports and two
  orphaned helpers (`_get_masses`, `_build_tr_basis`) in
  `mlmm/io/hessian_calc.py`.
- 4 mis-typed parameters annotated `Optional[...]` where the default is
  `None` (`MLMMCore.__init__: input_pdb / real_parm7 / model_pdb` and
  `build_model_pdb_from_bfactors: tolerance`).
- `docs/cli-conventions.md` now spells out the four permanent boolean
  forms (`--flag` / `--no-flag` / `--flag True/yes/1/on` /
  `--flag False/no/0/off`) and adds a "Contributing a new bool flag"
  section pointing at the `add_*_option` factory + `_COMMAND_BOOL_*`
  registries.

### Fixed
- mlmm-toolkit's default MM backend now declares `ninja` as a dependency. Its C++ kernels are
  JIT-compiled through `torch.utils.cpp_extension`, which on modern torch requires
  Ninja and no longer has a distutils fallback, so a clean-env install (e.g. a
  fresh HPC conda env) could fail at runtime with "native bonded extension is
  unavailable". Ninja ships as a pip wheel on every platform (incl. linux-aarch64),
  so the extensions now build out of the box; the dead `use_ninja=False` fallback
  was removed.
- OPC / TIP4P 4-point water with a virtual site (Amber `EPW`, element `EP`) is
  now read and parameterised correctly through the PDB/ASE input layer and the
  OpenMM MM backend (`computeVirtualSites`), instead of mishandling the massless
  extra point. (mlmm-toolkit's default MM backend still needs a 3-point model —
  see `docs/mm-parm.md`.)
- `orb_precision` now reaches `_OrbBackend` (the kwarg was being silently
  dropped via `**_kwargs`, so `--precision fp64` on the ORB backend
  always ran at the default `float32-high`). The legacy alias
  `"float32"` is rewritten to `"float32-high"` for backward compatibility.
- Degenerate ML/MM link distances (|r(MM) − r(ML)| < 1e-6 Å) now raise
  `ValueError` instead of silently dropping the link H. The previous
  `continue` left the link slot at its template (0,0,0) position,
  corrupting every downstream energy / force / Hessian. Same fix
  applied in `dft._append_link_hydrogens` for the ONIOM export path.
- `uma_precision` is now hidden from the echoed config for non-`uma`
  backends (it was missing from `_BACKEND_KEY_PREFIXES`).
- `_COMMAND_BOOL_SINGLE_FLAG_OPTIONS` is now wired through to the
  DefaultGroup, so `--auto-mm-keep-temp` is correctly classified as a
  single-flag bool by the `bool_compat` shim.
- `dataset_list` log filter is scoped to the `fairchem` logger subtree
  instead of the root logger (the previous code silently suppressed any
  log record across the process whose message contained `dataset_list`).
- `workflows/tsopt.py` heavy- / light-mode fallback no longer swallows
  energy + imaginary-mode failures with `except Exception: pass`. Narrow
  to the actual classes raised by `_calc_energy` /
  `_frequencies_cm_and_modes`, log via `logger.warning`, and emit a NaN
  sentinel for energy so downstream consumers can distinguish a missing
  value from a zero.
- `convert_and_annotate_xyz_to_pdb` docstring corrected: defaults are
  0 / 10 / 20 (ML / movable / frozen) per `annotate_pdb_bfactors_inplace`,
  not the 100 / 50 / 150 values the docstring previously claimed.
- ML-region charge/spin parity check hoisted to `MLMMCore.__init__`
  preflight, with a corrected off-by-one in `selection_indices` (the
  hoist initially treated them as 1-based, but `_mk_model_parm7`
  returns parmed 0-based atom indices; the off-by-one produced sum_Z
  values for a shifted atom slice, so e.g. `-q -1 -m 1` for a layered
  ARG-bearing pocket was rejected even though the actual ML electron
  parity was valid). A bad `--charge` / `--multiplicity` combination
  now fails in O(ms) instead of after the multi-second ML model load.

### Added
- `--deterministic` flag on every compute subcommand (`opt`, `tsopt`,
  `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`,
  `all`, `sp`) for bit-reproducible GPU runs (deterministic algorithms +
  an `index_reduce_` shim). Process-global, slower, and fails loud if the
  build cannot honour it; `MLMM_STRICT_DETERMINISTIC=1` is the env-var
  equivalent. Verified bit-identical energy and forces on uma / orb / mace.
- `docs/reproducibility.md` documenting the determinism / precision model
  and the per-backend reproducibility guarantees.
- `tests/test_help_grouping.py` locks the four-bucket `--help` section
  rendering + order.
- `MLMMCore.compute` / `mlmm` Calculator class / MLMMCore.__init__
  workspace setup all gained docstrings (previously empty).
- `result.json` / `summary.json` envelope now carries
  `schema_version: "1.0"` (bumped when the structure changes) and
  `write_result_json` mirrors every per-stage `result.json` payload to
  a sibling `summary.json` so MCP clients and agents can converge on a
  single filename across every subcommand. `result.json` is preserved
  for back-compat. `RESULT_JSON_STATUS_VALUES` enumerates the allowed
  `status` strings (`success` / `partial` / `error` / `unknown`).
- Structured error envelope when a subcommand fails: the JSON envelope
  now carries `error_class_chain` (MRO names so agents can match the
  exception hierarchy without parsing text), `error_module`, and
  `error_label` alongside the legacy `error` / `error_type` / `status`
  keys.
- `mlmm.mcp._runner` exposes `SubcmdResultDict` (TypedDict matching the
  runtime tool payload), `MCP_SUBCMD_RESULT_SCHEMA_VERSION = "1.0"`,
  and `MCP_SUBCMD_RESULT_STATUSES` enum. `SubcmdResult.to_dict` now
  emits `schema_version` so MCP clients can pin the contract.
- `docs/output-layout.md` (new): single-page reference for the filename
  conventions per subcommand + agent recipe for reading `summary.json`
  with class-chain pattern matching. `docs/json-output.md` and
  `docs/mcp_server.md` updated with the schema_version, summary.json
  mirror, and error envelope semantics.
- `tests/test_write_result_json.py`, `tests/test_error_envelope.py`,
  `tests/test_mcp_runner.py` (~9 tier-1 assertions) lock in the new
  envelope contracts so regressions show up at pytest time.
- `mlmm.workflows._all_helpers` (new module) hosts module-level helpers
  extracted from `workflows/all.py:cli()`:
    * `AllContext` (frozen dataclass) bundling the 72 `mlmm all` CLI
      parameters in declaration order.
    * `copy_path_outputs_to_root` / `promote_diag_for_root` (replace
      nested closures).
    * `build_energy_level_dict` factors the 5-way R/TS/P energy-level
      dict pattern (UMA / Gibbs / DFT / Gibbs-DFT) into one helper.
    * `build_pipeline_summary_payload` factors the summary-log payload
      assembly so the dict construction is unit-testable separately
      from the I/O wrapper.
    * `build_tsopt_overrides` / `build_freq_overrides` /
      `build_dft_overrides` factor the inline override-dict assembly
      (`if x is not None: cfg[k] = ...` ladder) for the post-MEP
      TSOPT / FREQ / DFT stage calls.
    * `append_backend_forwarding_args` consolidates the 8 copies of
      the backend / charge-embedding / link / mm / cmap / args-yaml
      argv-build ladder shared by every cli()-internal subcommand
      forwarder. Preserves the "only emit `--no-embedcharge` when
      the user explicitly typed it" subtlety so a YAML
      `calc.embedcharge: true` is not silently overridden.
- `mlmm.core.utils.resolve_ml_layer_assignment` (new helper) collapses
  the ~30-LOC layer-resolution block that was duplicated verbatim
  across `scan.py`, `scan2d.py`, and `scan3d.py` (movable_cutoff
  gating, detect_layer fallback to explicit ML region, model_pdb /
  model_indices resolution, calc_cfg mutation). Single source of
  truth for the layer-assignment behavior.
- `tests/test_all_helpers.py` (9 cases) pins the extracted helper
  contracts (copy-outputs no-op, copy-outputs canonical artefact set,
  promote-diag none/rewrite, energy-level kcal projection +
  no-input-mutation, AllContext frozen + 72 fields, pipeline summary
  payload shape, AllContext signature drift guard).
- TS-opt microiteration now supports internal coordinate macro steps
  (`--coord-type dlc|redund|tric`) by running the MM-only micro relax
  on a fresh cart twin Geometry and projecting the converged positions
  back onto the macro geometry via `Geometry.reset_coords()`. Cart-only
  macro path is byte-equivalent to the previous behaviour. Removes the
  pre-existing pysisyphus cart<->internal roundtrip drift that aborted
  microiter runs with `assert_allclose` after a few macro cycles.
- Smoke `tests/smoke/run.sh` expanded with per-stage `--coord-type
  {dlc,redund,tric}` + `--precision fp64` + `--mm-backend openmm` +
  `--link-atom-method fixed` test coverage (test50a/d/g/j/k/m/n/o);
  test50 itself capped at `--max-cycles 5 --no-tsopt/thermo/dft` so
  the DLC code path is exercised without the multi-hour GSM
  convergence the uncapped run requires.
- `--precision fp32|fp64` accepted on every calculator-constructing
  subcommand. The flag was previously available only on `tsopt / freq /
  irc / sp`; it now also covers `opt / all / path-opt / path-search /
  scan / scan2d / scan3d`. For `all`, the value propagates to every
  child stage through the shared args YAML, so a single top-level
  switch covers the full pipeline.
- `--irc-pos-def` (IRC convergence guard requiring PSD mass-weighted
  Hessian) is opt-in on `irc`; blocks the IRC "shoulder" false
  convergence where the rms-only criterion calls success before
  reaching the local minimum.
- `mlmm.core.residue_data` (new module) hosts the `AMINO_ACIDS` / `ION`
  / `WATER_RES` tables shared between `workflows/extract` and
  `domain/add_elem_info`; removes the L3 -> L2 import inversion that
  previously had `domain/add_elem_info` reaching back into
  `workflows/extract`. The legacy `from mlmm.workflows.extract import
  AMINO_ACIDS, ...` form is preserved via re-export.

### Removed
- **BREAKING:** `mlmm pysis` subcommand and the `mlmm.pysis_runner` module.
  This was a thin wrapper that registered the `mlmm` calculator into
  pysisyphus's `CALC_DICT` and shelled out to the pysisyphus YAML runner,
  providing v0.1.x YAML-workflow compatibility (`mlmm opt.yaml`). That
  compatibility surface is dropped: drive runs through the `mlmm`
  subcommands (with `--config` for YAML-supplied defaults) instead. Using
  the `mlmm` calculator directly from Python (`from mlmm import mlmm`)
  is unaffected.
- **BREAKING:** Flat-top compatibility shim layer removed. The package now
  lives under 6 layer directories (`cli/`, `workflows/`, `domain/`,
  `backends/`, `io/`, `core/`); the shims at `mlmm/<file>.py` that
  re-exported the new locations have been deleted in this release. External
  code must migrate dotted imports to the layered paths:

  | Old (removed)              | New                                |
  |----------------------------|------------------------------------|
  | `mlmm.{all,opt,tsopt,freq,irc,scan,scan2d,scan3d,path_opt,path_search,extract,dft,mm_parm,oniom_export,oniom_import,define_layer}` | `mlmm.workflows.<same>` |
  | `mlmm.align_freeze_atoms`  | `mlmm.workflows.align_freeze`      |
  | `mlmm.scan_common`         | `mlmm.workflows.scan_common`       |
  | `mlmm.{defaults,utils}`    | `mlmm.core.<same>`                 |
  | `mlmm.{mlmm_calc,xtb_embedcharge_correction}` | `mlmm.backends.<same>` |
  | `mlmm.{bond_changes,bond_summary,add_elem_info}` | `mlmm.domain.<same>` |
  | `mlmm.{energy_diagram,trj2fig,hessian_cache,hessian_calc}` | `mlmm.io.<same>` |
  | `mlmm.harmonic_constraints` | `mlmm.workflows.restraints`       |
  | `mlmm.fix_altloc`          | `mlmm.io.pdb_fix`                  |
  | `mlmm.summary_log`         | `mlmm.io.summary`                  |
  | `mlmm.cli_utils`           | `mlmm.cli.decorators`              |
  | `mlmm.{bool_compat,default_group,preflight}` | `mlmm.cli.<same>`   |
  | `mlmm.advanced_help`       | `mlmm.cli.help_pages`              |

  All console-script subcommands are preserved across this move **except**
  `mlmm pysis` (removed, see below); `mlmm sp` was added. Other Python
  imports change as tabulated above.
- `--trust-band` / `--hessian-window` / `--weighted-trust` CLI flags
  (and their `add_*_option` factories). The trust-radius / multistep
  Hessian-update knobs they exposed showed no benefit on small TS
  benchmarks and actively slowed convergence (rho-band trust update
  −33 %, hessian_window > 1 −47 % cycles on a 20-atom TS), with no
  evidence of speed-up on production-scale systems. The vendored
  pysisyphus `HessianOptimizer` kwargs are left dormant; no
  behaviour change since defaults were always legacy.

### Documentation
- Documented that the `[orb]` extra's `torch_scatter` has no PyPI binary wheel
  (sdist only, fails under PEP517 build isolation): install from PyG's
  prebuilt-wheel index, e.g.
  `pip install "mlmm-toolkit[orb]" -f https://data.pyg.org/whl/torch-2.8.0+cu129.html`.
- Documented the default MM backend's 4-point-water (OPC/TIP4P) limitation and the
  3-point (OPC3/TIP3P) / `--mm-backend openmm` workarounds in `mm-parm.md`.
- Removed stray internal authoring notes from in-source comments and
  docs; technical rationale preserved verbatim, no runtime change.

## [0.2.9] — unreleased

Consolidated changes since v0.2.4, pending the `v0.2.9` tag. Highlights: end-to-end
JSON output, MM-only optimization mode, GPU memory headroom for 16 GB
consumer cards via Hessian-update CPU offload, 1-based atom indices in
all user-facing I/O, energy-plateau convergence fallback, and a
comprehensive documentation overhaul (EN/JA).

### Added — CLI features
- `mlmm bond-summary`: detect bond changes between two structures
  (positional args supported, R/P sanity check companion to `mlmm all`).
- `mlmm opt --mm-only`: skip the MLIP component and minimize on the MM
  force field only. Layers are still honored; `--opt-mode hess` rejected
  (MM-only calculator does not provide a Hessian); microiteration is
  auto-disabled. Useful as a cheap MM pre-relaxation before ML/MM
  ONIOM optimization.
- `--out-json` across all MLIP subcommands (`opt`, `tsopt`, `freq`,
  `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`): emits a
  machine-readable `result.json` per run, including backend, charge,
  spin, and timing. `mlmm all` migrates `summary.yaml` → `summary.json`.
- `--link-atom-method scaled|fixed` on all computation subcommands.
- `--cmap / --no-cmap` to exclude CMAP from the model parm7 (Gaussian
  ONIOM compatibility).
- `--hess-device`, `--read-hess`, `--dump-hess`, `--skip-final-freq`
  for explicit Hessian device control and serialization.
- `--engine gpu|cpu` for `mlmm dft`; `--lowmem` selects
  `gpu4pyscf.rks_lowmem.RKS` as the closed-shell default.
- `--modified-residue` for `mlmm extract` / `mlmm all` (non-standard
  amino acid handling).
- `--freeze-atoms` for `mlmm irc`.
- ML-region charge/multiplicity sanity validator that runs before the
  first MLIP evaluation, surfacing wavefunction-charge mismatches early.

### Added — Convergence and stability
- Energy-plateau convergence fallback for optimizers (range-based,
  `1e-4 au` threshold; skipped for chain-of-states optimizers).
- Pre-MEP global alignment of reactant / product structures.
- Auto ECP for def2 basis sets in `mlmm dft`.

### Added — Documentation & infrastructure
- Agent skill set under `skills/` covering install, structure
  I/O, CLI subcommands, workflows / outputs, and HPC submission.
- ChemRxiv preprint added as `preferred-citation` in `CITATION.cff`
  (DOI placeholder until Zenodo release).
- xTB install documentation and improved not-found error message.
- JSON Output Reference page (EN + JA).
- GPU4PySCF Blackwell-class OOM workaround note in `dft.md` (EN + JA).
- `mm-parm` antechamber sqm odd-electron pitfall.

### Changed
- Atom indices unified to **1-based** for all user-facing I/O
  (`--freeze-atoms`, layer atoms, scan target display, etc.) with
  corresponding EN/JA documentation updates.
- `--refine-path` default `True` (path-opt is the default MEP mode).
- `--exclude-backbone` default `False`.
- RFO / RS-I-RFO `trust_max` reduced from `0.20` to `0.10` for MLIP
  stability.
- `--max-nodes` default `20` for `path-search` / `path-opt`.
- `--thresh` default `gau_loose` (was `gau`).
- Unified `ml_hessian_mode` and `hessian_calc_mode` under a single
  `hessian_calc_mode` setting.
- Imaginary mode filenames shortened (`final_imag_mode` → `imag`).
- Shorter / cleaner CLI logging with deduplicated blank lines, relative
  paths in echo, and suppressed empty config blocks.

### Fixed — GPU memory
- `pysisyphus/optimizers/hessian_updates.py`: `bofill_update` runs on
  CPU for `torch.Tensor` input, avoiding a ~5 GB GPU peak (sr1 / psb /
  mix temporaries each ~1.35 GB for ~4000-atom active DOF). Result is
  transferred back to the Hessian's device. Pairs with the IRC fix
  below; together they bring active-DOF 4000–5000 IRC into the 16 GB
  consumer-GPU envelope.
- `pysisyphus/irc/IRC.py`: stash `forward_mw_hessian` on CPU during
  backward integration (consumer transfers it back to numpy at the end
  anyway). Frees ~`N_dof² × 8 B` of VRAM for the entire backward run.
- All workflow CLIs (`opt`, `tsopt`, `freq`, `irc`, `scan`,
  `path-search`): explicit `del` of heavy locals (`calc`, `optimizer`,
  `H_t`, `modes`) in `finally` before `gc.collect()` and
  `torch.cuda.empty_cache()`. Fixes the post-tsopt → freq OOM observed
  in `mlmm all` on bezA-class systems with 16 GB GPU (cyclic-gc could
  not break torch.nn.Module hook / closure cycles).

### Fixed — Correctness
- ML-region `model_charge` / `model_mult` now follow the CLI-resolved
  value (`-q` or `-l` derivation) across `opt`, `tsopt`, `freq`,
  `irc`, `scan`, `path-search`, `path-opt`. Previously the
  `MLMM_CALC_KW` default `0` silently overrode user input because
  `dict.get(key, fallback)` returned the present default rather than
  the CLI value; wavefunctions were computed against charge 0
  regardless of `--ligand-charge` (UMA's robustness masked downstream
  symptoms). Single-line fix replaces the get-with-fallback dance with
  direct assignment from CLI-resolved `charge` / `spin`.
- `mlmm/all.py`: honor input PDB B-factor layers when active-site
  extraction is skipped (`-c` not supplied with pre-layered input).
- `mlmm/mlmm_calc.py`: use 1-based ATOM/HETATM file position as `idx`,
  not raw PDB serial. PDBs with serial gaps (e.g. 3411→3418) no longer
  trigger `IndexError` in `_mk_model_parm7`. Regression test:
  `tests/test_mlmm_calc_serial_gap.py`.
- Microiteration oscillation with scaled link atoms.
- IRC initial-displacement bisection (in-place mutation + missing
  numpy import).
- EulerPC corrector Euler integration safety guards (zero / NaN
  gradient norms).
- Tangent normalization and SVD alignment NaN / zero guards.
- DFT failure path: graceful skip of energy diagrams + status in
  summary; error `result.json` written on subcommand failure for all
  CLI commands.
- `bond-summary` PDB loading (`geom_from_pdb_str` → `geom_from_pdb`).
- `_resolve_device`: handle `'auto'` → `cuda` / `cpu`.
- PDB element inference: text-based coord replacement + robust atom
  name parsing for protein hydrogens.
- `mlmm/all.py` modified-residue handling: `AMINO_ACIDS` restore
  wrapped in `try / finally`.
- MODEL / ENDMDL missing on first frame in PDB trajectory conversion.

### Fixed — Documentation hygiene
- `skills/`: factual cleanup of trajectory file names
  (`opt_trj.xyz` → `optimization_trj.xyz`), default-dict references
  (`UMA_CALC_KW` → `MLMM_CALC_KW`), TS Hessian flag (`--hessian-init`
  → `--hessian-calc-mode`), output-tree paths, status enum values, and
  removal of phantom `.log` rows. Cosmetic-only changes to the AI-side
  cheatsheets; the docs/source are the canonical references.

### Documentation
- Bilingual documentation overhaul (EN/JA) covering terminology,
  accuracy, and consistency: net charge (not total), pocket → ML
  region determination, total → net charge wording, link-atom and
  microiteration concepts, EN/JA alignment for 1-based atom indices.
- "At a glance" 5-bullet block (Use when / Method / Outputs
  / Defaults / Next step) added to the 10 core subcommand pages
  (`opt`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`,
  `path-opt`, `path-search`, `dft`).
- `docs/yaml-reference.md` (EN+JA): document `mlmm:` block as alias
  for `calc:` (10 callsites verified); add `dft.engine`, `dft.ecp`,
  `freq.active_dof_mode`; correct lbfgs/rfo Used-by columns; drop
  dead `stopt.rfo` example.
- `docs/cli-conventions.md` (EN+JA): drop `< override-yaml` from
  precedence chain (slot kept in `cli_utils.py` for legacy compat but
  no CLI subcommand exposes the flag); `--opt-mode` alias coverage
  (`grad / hess / light / heavy / lbfgs / rfo`); em-dash convention.
- `docs/json-output.md` (EN+JA): replace "every MLIP-based subcommand"
  `--out-json` claim with explicit list (excludes `path-search`,
  `all`, `bond-summary`, `define-layer`, `mm-parm`); enumerate the
  6-alias `opt_mode` Choice; add `n_segments_reactive` field for
  `mlmm all`.
- Reference docs regenerated from current CLI; markdown filenames
  normalized from underscores to hyphens; theme toggle simplified to
  light ↔ dark.

## [0.2.4] — 2026-03-18

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

- **Optional backends**: `pip install "mlmm-toolkit[orb]"`, `"mlmm-toolkit[aimnet]"`. MACE (`mace-torch`) conflicts with `fairchem-core` (UMA) due to incompatible `e3nn` versions; use separate conda environments if both are needed.
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
