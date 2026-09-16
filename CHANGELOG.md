# Changelog

All notable changes to **mlmm-toolkit** will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

_No changes yet._

## [0.3.5] — 2026-09-16

### Breaking changes

- JSON schema 3.0: replace IRC `forward_converged` / `backward_converged` with `forward_status` / `backward_status` (`stopped`, `failed`, `disabled`). Integration and downhill-departure diagnostics remain available.
- `all` requires converged, topology-validated endpoints. Missing or unusable required stages now report `partial` or `failed` rather than `success`.
- `search.max_depth` now counts subdivision levels: `0` disables subdivision; the same value permits one level fewer than before.
- TSOPT reports `energy_missing` when final energy evaluation fails; `optimization_status` retains the optimizer verdict.

### Added

- Expose `--max-depth` for recursive path search and report executed optimizers, convergence and stop reasons in JSON and logs.

### Changed

- Cartesian RS-P-RFO defaults to TS-BFGS updates and a 0.1 Å maximum-atom initial/maximum trust radius; explicit settings take precedence.
- Support fairchem-core 2.22. DMF reports IPOPT convergence consistently (success codes 0 and 1), and its default iteration limit increases from 300 to 3000.
- Retain all signed physical frequencies and mode vectors, including soft positive modes, in exports and thermochemistry. Thermal corrections may change; the 5 cm⁻¹ reporting threshold remains separate from raw curvature checks.

### Fixed

- Improve Cartesian OPT/TSOPT Hessian updates, trust-region steps and final curvature checks, including no-flatten continuation from higher-order saddles. Reduce constrained-Hessian memory use and correct stale frequency-cache reuse.
- Fix the UMA first-call CUDA mismatch reported in [pdb2reaction #298](https://github.com/t-0hmura/pdb2reaction/pull/298).
- Fix ORB installation on Python 3.13 Colab runtimes; prefer PDB results and trajectories, and correct scan command display and status messages.
- Use standard optimization with embedding; preserve microiteration failure diagnostics and correct final mode-index handling. Restore AIMNet2 evaluation and Plotly image export.

### Documentation

- Update workflow and installation guidance.

## [0.3.4] — 2026-09-02

### Added

- Add full-system COMT and refined-path BezA Endpoint examples to the Colab GUI.

### Changed

- Improve the Colab workflow guidance, built-in example walkthroughs, and
  linked setup and runtime help.

## [0.3.3] — 2026-08-29

> Upgrade warning: unchanged inputs can produce different geometries, energies/barriers,
> vibrational classifications, thermochemistry, and scientific/terminal status. Users of
> `result.json`/`summary.json` must review the Breaking changes and output-schema updates below.

### Breaking changes

- Remove the public `--tr-projection` option and the `legacy-active` treatment.
  Frozen-boundary PHVA now always uses the constrained treatment; stale YAML
  values fail explicitly.
- Restrict `path-opt` / `path-search` / `all` `--thresh` to single-structure
  optimizations. It no longer sets the GSM string-optimizer preset, which now
  has its own `--thresh-gsm`. A command that relied on `--thresh` to tighten or
  loosen the string optimizer must pass `--thresh-gsm` as well.
- **CMAP now defaults to the force-field-faithful policy in both MM layers.**
  Parm7 CMAP terms are preserved in both REAL and MODEL calculations, so complete
  model-internal terms cancel in `E_real_low + E_high - E_model_low` while boundary
  terms remain part of the low-level coupling. `--no-cmap` is an explicit
  modified-force-field opt-out and removes CMAP from both layers. Earlier releases
  removed CMAP only from MODEL by default, which mixed two low-level Hamiltonians.
- **`add-elem-info` no longer overwrites its input by default.** Omitting `-o`
  now writes `<input>_add_elem.pdb`. Pass `--inplace` to replace the input;
  `--overwrite` continues to mean re-infer existing element fields.
- **JSON schema 2.0 (breaking).** UMA-specific summary keys became backend-neutral MLIP keys
  (`post_segments[].uma` → `.mlip`, `gibbs_uma` → `gibbs_mlip`, `gibbs_dft_uma` → `gibbs_dft_mlip`);
  the old keys were removed. This is `schema_version: "2.0"`. Update parsers before upgrading.
- **Removed the `add_uma_precision_option` Python alias** (`mlmm.cli.common_options`). It previously
  aliased `add_precision_option`; external scripts importing the old name break with `ImportError`.
  Migration: import `add_precision_option`.
- **MCP tools now reject `extra_args` that override a managed output option.** Passing `-o`/`--out-dir`
  (or another MCP-managed output switch) through `extra_args` now raises `ValueError` instead of being
  forwarded verbatim. Migration: use the tool's own output parameters.
- **MCP `search_paths` now requires `product_pdb` (keyword-only).** The tool previously
  sent one structure to a CLI that requires at least two. Migration: pass the product
  endpoint explicitly; optional `intermediate_pdbs` are inserted in reaction order.
- **Hessian handoffs now verify electronic state.** `freq --dump-hess` writes
  schema 2 with model charge and multiplicity. `irc --read-hess` rejects a
  state mismatch; schema-1 files require explicit
  `--allow-unverified-hess-state` after independent verification.
- **ONIOM reference imports now verify atom order.** Exported inputs carry an
  atom-order digest. `oniom-import --ref-pdb` rejects malformed or mismatched
  markers; markerless files with repeated elements require
  `--allow-unverified-ref-order`, which cannot bypass a known mismatch.
- **The energy-plateau stop is now opt-in and off by default**
  (`opt.energy_plateau: false`), and it never applies to the MM micro
  iterations of a `--microiter` run. A flat energy is not evidence of a
  stationary point, so a run that stopped there reported `stalled` while
  leaving the geometry unconverged; in the micro step a flat MM energy with
  forces still above threshold ended the macro/micro alternation with the
  environment unrelaxed, which left extra imaginary modes in the TS search.
  `--max-cycles` (macro) and `microiter.micro_max_cycles` (micro) remain the
  real bounds. Turn the macro stop back on with `--stop-plateau` (see Added);
  YAML `opt.energy_plateau: true` also still enables it.

### Added

- Add `--gsm-param {equi,energy}` to `all`, `path-opt`, and `path-search` as an advanced GSM node-parameterization control; `equi` remains the default.
- Add `--stop-plateau/--no-stop-plateau`, `--stop-plateau-thresh`, and
  `--stop-plateau-window` to `opt`, `tsopt`, and `all`, exposing the
  energy-plateau stop and its two tuning values on the command line.
- Add `--thresh-gsm` and `--thresh-dmf` to `all`, `path-opt`, and
  `path-search`. The former selects the GSM string-optimizer preset; the latter
  selects the DMF IPOPT dual-infeasibility tolerance (`tight`, `middle`,
  `loose`, or a positive float).
- Add `all --tsopt-from-mep-tan/--no-tsopt-from-mep-tan`. The MEP tangent guides the
  initial TS root by default; when disabled, TSOPT selects from the initial-structure
  Hessian modes.
- Add a runnable full-system BezA example with reactant, intermediate, and
  product structures plus endpoint-MEP and staged-scan workflows.
- Report citations for the methods actually used at the end of `summary.log`
  and final stdout, and expose the same `{method, citation, doi}` records as
  `summary.json.references`.
- Write directly inspectable ML-region structures both before and after link
  hydrogens generated from parm7 ML/MM boundary bonds, including PDB companions
  for PDB input.
- Install AmberTools in the Colab notebook through an isolated Miniforge, so
  `mm-parm` and the DMF path mode run there; the Colab Python environment is
  left untouched.
- Cancel a running Colab job from the interface, and send an occupied output
  directory to `result(1)`, `result(2)`, … instead of writing into it.
- Add an mmCIF/large-PDB bridge (atom-identity–preserving; multi-model input keeps the first model,
  with a warning), exact selectors, safe duplicate atom names, and root/segment CIF companions with
  original identifiers.
- Add `tsopt --ref-mode`, opt-in IRC never-stop traversal, and analytical
  Hessians for ORB, MACE, and AIMNet2.
- Add a release-pinned, keyboard-accessible Colab GUI that preserves full-system
  coordinate/topology identity, accepts a separately prepared model, builds workflow-aware
  commands, reaps interrupted jobs, and limits result views/downloads to the current run.
  Its compact Input, Viewer, Options, and Results workspace co-locates workflow
  selection with a 4:3 molecular view. Exact 3D picks use Mol*'s native focus and
  selection behavior without rebuilding the viewer. The GUI keeps primary-input
  selections separate from view-only secondary structures and rolls back failed structure
  switches. Center selectors retain chain and insertion-code identity;
  incompatible secondary views suppress primary selection and measurement overlays.
  Validation follows the exact command and content hashes for every existing-file input
  declared by the selected Click command, while its transcript remains separate from
  current-run diagnostics. Private runtime files avoid upload collisions; Results links one
  trajectory control to the 3D frame and energy-profile cursor, while keeping structures,
  diagrams, SVG, HTML, CSV, PDF, and other generated files in a separate preview. Missing
  energies, failures, cancellations, and diagnostic bundles remain explicit.
  The input queue appends uploads as removable rows, keeps water visible on request,
  and derives searchable per-option controls and click-to-open help from the selected
  live CLI. Key and advanced controls remain collapsed until needed, while the
  generated command line stays visible.
- Add `opt`/`all --reject-uphill/--no-reject-uphill` (default off). On `all`,
  it is forwarded to the post-IRC endpoint re-optimization child only.
- Add `all --irc-step-size` so the end-to-end workflow can forward a smaller
  EulerPC step to every post-TS IRC branch.
- Detect the molecular point group and external rotational symmetry number for
  every `freq` structure and include the rotational `1/sigma` correction
  automatically. `thermo.symmetry_number` remains an advanced YAML override.
- Add selectable UMA/ORB/MACE/AIMNet2 frame rescoring to `trj2fig`, with model,
  precision, and provenance controls. Comment-energy mode remains
  calculator-free, and rescoring is a pure MLIP calculation rather than ONIOM.
- Announce the first load of each ML backend model with
  `[backend] Preparing MLIP model (...)...` and `[backend] Done.`.

### Changed

- Remove the repeated status/image-count/backend footer from the Colab summary table; the compact run context owns status and artifact details.
- Remove the redundant `Path with N moving images.` startup line; tagged GSM sections and `String=...` records identify progress.
- Rename the prepared Toy inputs to `r_toy.pdb`, `p_toy.pdb`, and `p_toy.parm7` and label the notebook route as MEP mode.
- Enable repeated trajectory playback by default, select stitched `finished_irc_trj` for IRC profiles, and place the trajectory/energy view before result status and generated-file details in Colab.
- Show Hessian cache-reuse notices at `-v 2`; cache identity and rejection details remain at `-v 3`.
- Map the `hess` TS-optimizer preset to RS-P-RFO. Standalone `tsopt` retains RS-I-RFO through `--opt-mode rsirfo`.
- Remove the legacy `light` and `heavy` optimizer aliases; use `grad`/`hess` or the algorithm names exposed by each subcommand.
- Keep finite product cycle defaults: 100000 for ordinary optimization,
  300 for GSM/DMF, 125 for IRC, 100000 for ML/MM microiterations, and 100 for
  DFT SCF. An explicit YAML `null` remains the uncapped engine value.
- Detect B-factor ML/MM layers by default and expose
  `--detect-layer/--no-detect-layer` consistently across commands and `all`
  child stages. Explicit `--model-pdb` / `--model-indices` membership remains
  supported.
- Treat `--embedcharge` as an experimental, opt-in, computationally expensive
  feature. MLIP/MM commands apply the xTB point-charge delta correction; `dft`
  adds MM point charges directly to the PySCF QM Hamiltonian.
- Report every option's effective default. Options whose real default lives in
  a config block are declared `None` so an explicit value stays distinguishable
  from an omission; each now carries that default as a display string, so
  `--help`, `--help-advanced`, the generated reference and the Colab Options
  pane stop reading as unset. The Colab controls -- dropdowns included -- label
  and select `default: <value>`, or `default: None` when there genuinely is none.
- Drop the notice that replaced the preparation panel for workflows that
  extract internally; the panel is simply hidden.
- Name the custom-calculator entry point `--calc-file-func-name`. It names a
  callable inside `--calc-file`, which the previous spelling `--calc-factory`
  left unsaid. `--calc-factory` is removed rather than aliased.
- Head the stdout citation block `====== Citations & References ======` like
  every other console section. `summary.log` keeps its numbered
  `[6] Methods and citations`; the two share one renderer, so the log file's
  section index was leaking into the console. Each method and its numbered
  reference are now printed on separate lines.
- Reconcile a completed Colab run from the browser heartbeat when its final
  widget update is lost, restoring the controls and Results without rerunning
  the GUI cell.
- Keep one execution path in the Colab GUI. An unreachable async twin of the
  whole run path (`_start_async_task`, `_do_validate_async`, `_do_run_guarded`,
  `_do_run_async`, `_stream_async`, `_stop_async_process`,
  `_validate_command_async`, `_async_task_done`) was never called, so a reader
  had two implementations to reason about and only one ever ran. Removed with
  three unreferenced one-liners; the notebook cell is 264 lines shorter.
- Correct the `baker` description in the `opt` skill page: convergence requires
  all four force/step thresholds AND the energy change, not the published
  force-and-either form. The contract gate pinned the stale sentence.
- Rename the Colab preparation button to "Extract ML region & use it", and put a
  one-line hint where the panel used to vanish: `all` extracts internally and
  `extract` is itself the extraction command, so the panel is hidden for both
  and the hint names the route to a standalone `--model-pdb`.
- Make uphill-trial rejection opt-in for L-BFGS and RFO minimization. When
  explicitly enabled, its energy tolerance is `1e-4` Hartree. TS optimization
  always leaves the rejection disabled.
- Restore ordinary IRC's strict energy-rise stop: any positive one-step rise
  stops that direction. `--irc-never-stop` bypasses this and the other physical
  endpoint criteria until the cycle cap.
- Derive the ML-region charge in the Colab notebook from `--ligand-charge`
  instead of requiring a confirmed `-q`; a ticked charge box is now an explicit
  override.
- Install the DFT extra by default in the Colab notebook and verify plot export
  by rendering a PNG. A first run takes several minutes and needs a GPU
  runtime.
- Label thermochemistry as `E + G_corr = G`, force uphill rejection off for
  transition-state optimization, and keep its toggle limited to minimum and
  post-IRC endpoint optimization.
- Apply the `baker` preset as a deliberately tightened variant of the published
  criterion (Bakken and Helgaker, J. Chem. Phys. 117, 9160 (2002)): maximum
  force, RMS force, maximum step, RMS step and the energy change must all hold,
  including the final check of a retained lower-energy geometry at the
  uphill-rejection trust floor. A zero-length step settles the energy criterion,
  because the geometry cannot move.
- Stop forcing one reparametrization pass per growing-string cycle. An image
  already inside the parametrization tolerance is no longer displaced, and
  coincident parameter densities are rejected instead of divided.
- Pin backend setup recipes to the official PyTorch 2.8 wheel matrix, install
  dedicated-environment MACE only after removing `fairchem-core`, and make HPC
  templates fail fast unless the `hessian_ff` JIT compiler prerequisites are
  present.
- Remove the unused internal `AllContext` parameter mirror and break the product
  import cycles (`core.utils`↔`extract`, `freq`↔`opt`) by relocating the shared
  charge/spin preparation and layer helpers; the relocation itself makes no CLI or
  JSON contract change (see the `sp` ML-region resolution change below, which does
  move numbers).
- Reject explicit analytical Hessians with `workers > 1`.
- Project only rigid modes that are an actual null space of the frozen system for
  PHVA, IRC, Dimer, and TS validation, and record the effective mode. The former
  active-fragment projection could hide a real imaginary mode, so `n_imag`, ZPE
  and ΔG‡ move on frozen-boundary systems. The superseded public
  `--tr-projection` option and `legacy-active` treatment are removed.
- Require `n_imag = 1` for TS success, preserve rejected optimizer state, and
  report resolved backend/model/precision and the highest common rate-limiting method.
- Keep the default Hessian TS search on standard restricted-step root
  following. Trial mode-loss rejection, intermediate eigenvalue-structure
  gating, automatic saddle recovery, and automatic displaced multistarts are
  disabled by default. Exact PHVA remains the terminal authority and requires
  `n_imag = 1`; explicit `--flatten` remains available for extra modes.
- Double the default segment path resolution (`max_nodes_segment` 10 → 20), which
  changes the MEP, its highest-energy image, and therefore the reported barrier.
- Unify the residue/ion/water catalog so charge inference recognizes the same
  ions as element inference; charge summaries now include monatomic ions that
  were previously unrecognized.
- Populate missing element columns in the topology-matched PDB exported by
  `mm-parm`, preserving LEaP atom records and order; extraction recovery text
  now names the non-destructive output or explicit `--inplace` route.
- Keep finite-difference Hessian assembly, low-rank Bofill updates, and mass
  scaling device-resident on GPU runs, avoiding per-step host round-trips.
- Make the `tsopt` Hessian-radius contract explicit: an omitted
  `--radius-hessian` includes the movable MM atoms required by the default
  `partial` final-frequency basis, while an explicit narrower cutoff is rejected
  before optimization unless `--active-dof-mode` requests only evaluated atoms.
- Make periodic IRC HDF5 trajectory checkpointing opt-in
  (`irc.dump_every: null` by default). Enabled checkpoints contain coordinates,
  energies, and gradients, but never a dense Hessian.
- Honor YAML `opt.thresh`, the `dft` method, and `--func-basis` in `all`: these are
  forwarded to the post-MEP TSOPT and R/P endpoint-optimization children only when set
  explicitly, so a YAML value is no longer clobbered by a hardcoded default. Endpoint
  geometry/energy and DFT numbers change when configured.
- Build the align/refine and `--no-tsopt` HEI-probe calculators from the resolved
  calculator template, so `uma_precision`/`uma_model`/`workers`/`cutoff` and
  supported YAML `calc` settings now reach them; aligned geometry and HEI
  energy/barrier move when these are set.
- Skip the BFGS Hessian update when the curvature `s·y ≤ 0` (was applied), and force
  `use_active=False` on internal-coordinate frequency analysis, fixing a partial-Hessian
  index mismatch reachable via `--coord-type` with frozen atoms.
- Decide IRC endpoint minimality on an exact Cartesian Hessian under the opt-in
  `--irc-pos-def`, replacing the integrator's quasi-Newton Hessian.
- Mask the harmonic-restraint Hessian on a separate `hessian_constrained_atoms` set
  from the forces, changing frequencies where restraints and frozen atoms interact.
- Seed `scan2d`/`scan3d` grid points only from explicitly converged relaxations,
  reference relative and minimum energies to seed-eligible points only, and honor YAML
  `opt.thresh`/`opt.max_cycles` (previously clobbered by the CLI default); relative
  energies shift when any grid point fails.
- Isolate Direct-Max-Flux configuration per invocation (`fresh_dmf_config` deep-copy).
- Load `path-opt --ref-pdb` geometry from the original unrounded `.xyz` (was: a PDB
  rounded to 0.001 Å, then reloaded), removing a coordinate round-trip.
- Resolve the `sp` ML region through the shared validator instead of a bare
  `except`-and-use-the-whole-system fallback: a resolution failure now raises rather
  than silently treating the full input as the ML layer (which changed the ONIOM energy).
- Take `dft` QM-region atom indices from the one-based file ordinal (was the deposited
  PDB serial and `name[0]` element), so the QM region and link-hydrogen pairs shift on
  gapped serials or two-letter elements.
- Drop the ±1 terminal-charge correction for Amber-capped C-/N-terminal residues in the
  `extract` charge summary, changing the net protein/active-site charge for capped termini.
- Attach the built wheel and source distribution to each GitHub Release before
  trusted PyPI publication.
- `result.json`/`summary.json` gained additive field families — a `run_id` (from
  `MLMM_RUN_ID`; a conflicting id raises `RunIdentityError`); `execution_status` and
  `scientific_status` with reasons and item/expected/observed ids; per-stage and
  per-point outcomes; per-segment `converged`, `irc`, `endpoint_opt`, `ts_imag`, and
  `dft_status`; scan `energy_reference`/`n_points_usable`; a serialized
  `microiteration` object in `opt`/`tsopt` output; a `thermo_policy` block in the
  frequency YAML/`result.json`; and resolved provenance. These are additive for
  consumers that tolerate unknown fields. `scientific_status` also participates in
  usability/promotion decisions, not only provenance. Aggregate success requires
  every applicable producer convergence signal; direct-TS segments do not invent
  an MEP gate.
- `key_output_files` now lists only the artifacts claimed by the current
  invocation's run manifest rather than files discovered under the output tree, so a
  reused `-o/--out-dir` no longer reports stale files from an earlier run.
- Frequency JSON/YAML records `symmetry_number` and
  `symmetry_number_source`; `all` copies complete child provenance into each
  post-segment's `thermo_symmetry` map for R/TS/P. IRC result JSON records
  `electronic_state_verified` for a reused frequency Hessian.
- The MCP tool-return envelope moved from `schema_version` `1.0` to `1.1`, adding a
  per-invocation `run_id`, a `summary_run_mismatch` status, and a run-id byte check
  (distinct from the summary `schema_version: "2.0"`).
- Commit `result.json`/`summary.json` and converted structures by staged atomic
  replace and raise on a write failure that was previously swallowed; the mmCIF/PDB
  bridge hard-fails on out-of-range or non-finite coordinates and unresolved elements
  rather than emitting a corrupted fixed-column record. Valid-input output bytes are
  unchanged.

### Fixed

- Route the first Colab `Load results` action through the native browser bridge.
- Organize Colab Results into `Energy profile & Trajectory` and
  `Imaginary frequency`, with a contextual `View` or `Mode` selector.
- Install declared wheel dependencies during release validation and allow a
  published release tag to be retried safely after a workflow failure.
- Restore the intended Mol* expanded-view defaults, clear the manual-command
  notice on Rebuild, and keep combined-IRC endpoint labels below their markers.
- Initialize generated Colab numeric controls from each option's effective CLI
  default instead of its Click range minimum, and omit the flag again when the
  control is restored to that default.
- Keep energy-weighted GSM parametrization finite across minimum-energy and flat intervals while rejecting invalid energies.
- Correct CLI/YAML/default documentation, runnable output recipes, TS-only guidance, and canonical toggle names.
- Pin Torch 2.8.0 in every minimal no-deps CI lane and match the primary citation title.
- Bound Plotly/Kaleido static-image export to ten minutes in an isolated
  renderer process so a stalled Chrome startup or shutdown cannot block a
  scientific workflow indefinitely.
- Disable polynomial line search inside the Hessian Dimer and reject attempts
  to enable it, because the effective Dimer force is not conjugate to the
  reported physical energy.
- Keep IRC branches diagnostic when energy-based displacement cannot find a
  downhill departure; later small-gradient termination can no longer certify
  or cache such a branch.
- Preserve distance-derived frozen atoms when a calculator is attached, and
  constrain default Hessian targets to atoms inside `movable_cutoff`.
- Propagate requested endpoint-preoptimization convergence into standalone and
  `all` path scientific status while retaining diagnostic path continuation.
- Make `all --dry-run` validate extraction/scan syntax and the ML-region
  charge/multiplicity state before reporting success.
- Forward the resolved embedding state and cutoff together when YAML enables
  embedding and CLI overrides only the cutoff.
- Add `all --freeze-atoms`, merge it with YAML `geom.freeze_atoms` and the
  Frozen-MM layer across every child stage, expose it in Colab `all`, and keep
  the unsupported selection hidden only for `dft`.
- Honor an explicit `all --opt-mode` for both TSOPT and post-IRC endpoint
  optimization when `--opt-mode-post` is omitted.
- Separate the Colab workspace-path and example actions, restore `opt --dump`
  trajectories from the `[command]` entry in `run.log`, and persist `run.log`
  for ordinary output-directory CLI runs.
- Keep `dft --embedcharge` as direct PySCF electrostatic embedding and do not
  add the optional MLIP/xTB point-charge correction to the embedded DFT energy.
- Number TS candidates, validated TS structures, and intermediates across linked
  Colab MEP/IRC trajectories, and label combined IRC endpoints by direction.
- Label the Colab advanced-control disclosure as `Show All options`.
- Load an existing Colab-workspace file through the validated input queue and
  restore `./result_all/` whenever a built-in example is loaded.
- Show workflow warnings from `scientific_status_reasons` above Colab result tables while
  keeping Run details limited to execution metadata.
- Keep the complete MEP summary visible at default pipeline verbosity,
  including no-change results and every segment's barrier and reaction energy.
- Reconcile reported frequency modes against active DOF, rigid modes, and the
  near-zero window, and render reproduced commands shell-safely.
- Give generated text and log previews a bounded, visible scroll region in the
  Colab Results panel.
- Report `all` summary modes as `MEP` / `Scan` / `TS-only` and show absolute
  root and internal path-module output directories in `summary.log`.
- Match the interactive Colab energy-level segment and connector proportions
  to the CLI diagram and shrink overlong level annotations to fit their bars.
- Keep the scan2d/scan3d starting reference as the `-1` row in `surface.csv`
  while excluding it from energy baselines, interpolation, and plots.
- Compact the Colab run/output and Results status surfaces and enlarge
  scan-grid pick markers while retaining the scan3d color scale.
- Offset the scan2d coloured base plane above the z-axis floor to prevent
  coplanar rendering flicker.
- Activate the Colab Results tab before rendering completed results, keep the
  result selectors reusable, and reserve plot space for MEP/IRC labels.
- Keep cancelled Colab runs on the active tab and out of completed-result views.
- Report the exact MLIP model and UMA task in logs, JSON, and Colab Results, and
  include model-specific UMA, Orb-v3, MACE, and OMol25 citations.
- Use pause semantics for Colab trajectory playback.
- Restore the reviewed Colab selection UI: persistent removable chips, concise
  workflow labels and center summaries, explicit frozen-atom completion, and a
  compact staged-scan editor. The upload panel now treats parm7 as optional for
  `all`, which generates it when omitted.
- Fail closed on unknown execution, failed segment DFT, and non-converged path bridges, and clear stale scientific-status reasons.
- Preserve internal result identity while presenting aggregate structure-linked MEP/IRC profiles, reset prior result state before reload, and recover only declared or current-run artifacts.
- Preserve the TSOPT stop cause, restrict the output tree to current-run artifacts, and report the resolved frequency zero cutoff in `summary.log`.
- Label subtractive energies as ML/MM in human output and include compact three-layer counts.
- Keep optimizer verbosity monotonic and show Hessian cache reuse at `-v 2` without exposing raw DFT child diagnostics below `-v 3`.
- Translate internal partial-result codes into concise, actionable warnings in `summary.log` and final stdout.
- Preserve the blank line before the first MLIP model-load announcement even when an stderr warning immediately precedes it.
- Preserve the MM micro-optimizer until its non-convergence diagnostics have been collected, preventing a post-macro microiteration failure from raising `UnboundLocalError`.
- Record terminal PHVA as `skipped`, rather than `unavailable`, when a non-converged Hessian TS optimizer does not authorize the analysis.
- Print a section heading before each tagged recursive GSM segment and before standalone Growing String optimization.
- Bound the complete minimum-RFO line-search/GDIIS displacement by scaling it
  to the trust radius while preserving its accelerated direction.
- Enforce the same trust radius on hard-case Newton and reference RFO steps,
  preventing near-degenerate solvers from returning oversized displacements.
- Keep an explicit finite cycle cap when resuming a checkpoint written by an
  uncapped optimizer.
- Retain endpoint-optimization diagnostics when `--dump` is enabled.
- Gate MEP `--ref-mode` handoff to Hessian TS optimizers; Dimer records the handoff as not applicable instead of receiving an option it rejects.
- Preserve a non-converged or stalled TS final structure and stop before terminal PHVA.
- Use the configurable `freq.zero_cutoff_cm` value for standalone frequency analysis, flattening, and TS saddle classification.
- Restrict reference-aligned reaction-mode selection to negative exact-PHVA modes and validate the selected frequency before IRC. Invalid or missing selections use an explicit lowest-imaginary root-0 fallback with unverified reaction identity.
- Keep `--skip-final-freq` artifact-preserving but stop composite `all` before IRC because the negative reaction direction was not validated.
- Separate numerical optimization status from saddle order. A converged higher-order stationary point is not a first-order TS, but `all` may perform warning-labelled diagnostic IRC when a valid negative root exists.
- Synchronize EN/JA docs, skills, live help, generated command references, and checked-in contract tests with the current behavior.
- Remove the replaced path-tangent/single-mode helpers and one unreferenced mass-weighted-frequency wrapper; larger workflow and Notebook refactors remain deferred.
- Keep Hessian status and timing lines adjacent to the surrounding optimizer
  cycle rows instead of inserting blank lines before and after each Hessian
  evaluation.
- Run exact PHVA after every configured Hessian-TS convergence criterion passes
  for the proposed step, or after an opt-in energy-plateau stop. A run that
  instead exhausts `max_cycles` skips final PHVA and mode export, reports both
  imaginary-mode fields as `null`, and remains `not_converged`. A plateau runs
  PHVA but retains `stalled`.
- Keep RS-P-RFO line searches off by default while honoring explicit YAML
  values. RS-I-RFO and TRIM discard these unused kwargs, and Hessian TS searches
  require exactly one root.
- Apply the shared `opt` block to the `tsopt` macro optimizer. `rms_force`,
  `rms_force_only`, `max_force_only`, `force_only`, `overachieve_factor`,
  `min_step_norm`, `assert_min_step`, `converge_to_geom_rms_thresh`,
  `check_eigval_structure` and `line_search` were read from `opt` and then
  dropped, so a YAML convergence control was silently ignored by TS
  optimization while `opt` honored it. Both the microiteration and ordinary
  paths now build their macro kwargs through the same helper, and only values
  you changed are forwarded, so `rsirfo.*` stays authoritative for untouched
  keys and a default configuration is byte-identical to before.
- Read the `lbfgs` section in `tsopt`, so the microiteration MM relaxation's
  L-BFGS controls (`keep_last`, `max_step`, `line_search`, ...) are reachable
  from YAML. They were settable from neither YAML nor CLI.
- Report why a micro relaxation stopped. The micro runs with its stdout
  redirected, so the log said only `status=not_converged`; it now names the
  executed cycles against the bound, the last max force and step, and — when
  the bound is what ended it — `microiter.micro_max_cycles` as the setting to
  raise.
- Raise `microiter.micro_max_cycles` from 10000 to 100000 so the MM relaxation
  is bounded by its convergence criterion rather than by the cycle cap. The
  micro step is held to the same preset as the macro step, and under `baker`
  the few transient relaxations that follow a large rearrangement need more
  than 10000 L-BFGS cycles: measured on the release-smoke ML/MM TS lane, 791
  relaxations with a median of 56 cycles and three transients at 16815 / 13167
  / 10624. At the old cap those three reported `not_converged`, whose
  fail-closed handling aborted the TS search after 45 of 791 macro steps and
  left four imaginary modes; with the new bound every relaxation converges and
  the search reaches a certified first-order saddle. No convergence threshold
  changed, and the 788 relaxations that already converged are unaffected.
- Apply the whole shared `opt` block to the macro step of an `opt --microiter`
  run. It previously forwarded only `thresh`, `max_cycles`, `dump` and
  `out_dir`, so YAML convergence controls (`rms_force`, `rms_force_only`,
  `max_force_only`, `force_only`, `overachieve_factor`, `min_step_norm`,
  `assert_min_step`, `converge_to_geom_rms_thresh`, `check_eigval_structure`)
  and `line_search` / `print_every` were silently ignored in the default `hess`
  path while `--no-microiter` honored them. The macro step now uses the same
  merge rule as an ordinary run, so only values you actually changed are passed
  and a default configuration reaches the optimizer unchanged.
- Reject single-class B-factor metadata such as an all-zero PDB as an ML/MM
  layer partition instead of silently treating the full system as ML.
- Preserve valid B-factor movable/frozen MM layers when `sp` uses an explicit
  `--model-pdb` or `--model-indices`; invalid layer metadata no longer remains enabled.
- Honor `dmf.ipopt_options.dual_inf_tol` when it is set in YAML instead of
  replacing it with a fixed preset.
- Reject dependent Amber virtual sites before allocating the `hessian_ff`
  backend, and direct users to OpenMM or a three-point-water topology.
- Apply the inclusive B-factor tolerance consistently when reading ML,
  movable-MM, and frozen-MM layers; leave missing or malformed B-factors
  unassigned.
- Resolve the legacy `MACE-OFF23_small`, `_medium`, and `_large` aliases to the
  corresponding upstream MACE-OFF model sizes.
- Preserve selected-root diagnostics in non-Cartesian coordinates and convert
  torch eigenvectors explicitly for reference-mode overlap checks.
- Apply the tightened `baker` criterion on every evaluable cycle, including the
  first retained geometry.
- Convert Cartesian coordinates from Bohr to Angstrom before the
  thermochemistry hand-off; the moments of inertia, rotational entropy and
  absolute Gibbs energies were computed from Bohr values in an Angstrom contract.
- Require every internal-coordinate component to satisfy the back-transformation
  tolerance before the iteration is accepted; a single converged component used
  to end it.
- Apply the thermochemistry policy's `zpe_scale` (passed on as
  `zpe_scale_factor`) exactly once to the reported zero-point energy and to the
  vibrational internal energy. Values other than the default `1.0` were scaled
  twice.
- Collect invalid primitive indices from the actual dihedral and bend index
  lists instead of assuming contiguous ranges, keeping both sets.
- Correct the diagnostic QRRHO free-rotor partition function (frequency instead
  of wavenumber, plus the missing factor of pi); reported QRRHO entropy and Gibbs
  energy come from the Grimme interpolation and are unchanged.
- Evaluate the vibrational heat capacity through a series expansion for small
  exponents and `expm1` elsewhere, removing the cancellation error near zero.
- Accept equivalent leading-digit PDB and trailing-digit Amber hydrogen names
  while preserving strict atom-order checks for non-hydrogen atoms.
- Route `dft` through `MLMMCore` for topology preparation, MM calculators, and
  subtractive recombination. DFT now replaces only the high-level model energy
  and requests no forces; there is no separate model-parm7 builder or MM path.
- Keep the IRC running when the EulerPC corrector oscillates. The corrector
  descends the two-point interpolated surface rather than the real potential,
  so a reversal there is an interpolation artifact; it now warns and keeps the
  last non-oscillating point instead of aborting the whole run.
- Stop `scan2d` after writing `surface.csv` with a clear diagnostic when fewer
  than three non-collinear converged grid points remain, instead of passing an
  underdetermined data set to SciPy's RBF interpolator.
- Read the bundled pysisyphus optimizer's legacy comma-delimited Hartree XYZ comments in
  `trj2fig` and other strict trajectory consumers without treating unrelated
  numeric comments as energies.
- Expose the shared `--allow-charge-mult-mismatch` escape hatch on `tsopt`,
  matching the other ML/MM compute commands that run the same ML-region
  electron-parity validation.
- Require `scan2d --scan-lists` during Click parsing, reject unknown
  `energy-diagram` options instead of silently accepting them, and keep
  detail-only pipeline messages out of verbosity level 1.
- Reject unknown options and orphan arguments in legacy grouped-value commands,
  while accepting grouped or repeated path-search inputs and topology references
  consistently in shell and in-process Click invocations.
- Reject `all --scan-lists` with multiple inputs instead of silently selecting
  the endpoint-path route and ignoring the staged scan.
- Honor `tsopt --dump` when the initial microiteration MM relaxation fails closed
  before the first macro step; the combined trajectory now contains that executed
  micro relaxation instead of being absent.
- Write `thermoanalysis.yaml` from the `all` pipeline's `freq` stages by default, so
  `--thermo` actually yields thermochemistry. The child `freq` inherited its own
  `--dump` default (off), so `all --thermo` computed and printed the thermochemistry
  but never persisted the file the Gibbs assembly reads — every run reported
  `thermochemistry result is missing` and produced no Gibbs diagram. An explicit
  `--no-dump` still suppresses the file.
- Keep Hessian-Dimer orientations and off-center images on the frozen Cartesian
  constraint manifold, refreshing constraint-compatible rigid null modes at
  each central image.
- Keep scan energies on the unbiased PES and prevent stale Hessian reuse across
  TS, IRC, frequency, and endpoint optimization.
- Roll back rejected RFO/L-BFGS state, recover TS searches from `n_imag=0`, and
  keep path-guided flattening explicitly opt-in and mode-safe.
- Keep interpolated/GDIIS RFO displacements on the reference step's tensor
  device, preventing a CUDA-to-NumPy conversion failure during Hessian scans.
- Report a microiterated optimization's real terminal convergence. The macro
  loop's verdict was computed and then dropped, so a converged TS — including one
  whose exact-PHVA validation found a single imaginary mode — was written as
  `not_converged`, and `all` refused to start its IRC.
- Release the transition-state probe calculator before returning, so the IRC
  phase's leased core is the only heavy ML/MM core alive across the handoff.
- Run the `opt --flatten` loop when the optimizer stalls on an energy plateau.
  The loop rebuilds the Hessian and displaces along the remaining imaginary
  modes, so a stall is exactly when it is wanted; only a flatten *retry* that
  stalls again stops the loop.
- Keep optimizing when `opt.dump_restart` is set on an optimizer class whose
  restart state is not declared: the unsupported checkpoint is refused once,
  further dumping is disabled, and the run continues instead of aborting.
- Resolve the rigid-mode default from one place, so the `freq`/IRC/TS-optimizer
  fallbacks can no longer drift apart from the documented default.
- Correct IRC endpoint labels, `all` worker propagation, YAML custom-factory
  provenance, and CIF publication at pipeline root and segment level.
- Keep custom `--calc-file` selection consistent across child and in-process
  `all` stages when `--backend` is also supplied.
- Select one coherent altLoc per residue and identity-check frequency-to-IRC
  Hessian handoff against geometry, atom order, and active-DOF basis.
- Honor EulerPC's normalized IRC filename prefix in conversion, endpoint checks,
  and JSON; ship the third-party/OpenMM CMAP notices in built distributions.
- Apply the ZPE scale factor exactly once in thermochemistry, so a non-unity
  factor no longer enters the reported ZPE and `U`/`H`/`G` quadratically.
- Isolate each Dimer's random state from the process-global NumPy RNG, and raise
  a clear error for an invalid rotation method.
- Fingerprint the native `hessian_ff` build and runtime identity before loading a
  prebuilt extension, refusing a stale or host-incompatible binary.
- Stop reporting an electronic energy as a Gibbs free energy: `all` builds the
  per-segment MLIP and DFT//MLIP/MM Gibbs diagrams only when every state's frequency
  free energy (and DFT thermal correction) is finite, otherwise it skips the diagram
  and warns, instead of substituting the MLIP/DFT electronic energy or a `0.0`
  thermal correction. Reported ΔG changes wherever a thermochemistry value was
  missing.
- Compare the full per-atom signature when guarding `extract`'s multi-input atom
  order, so a swapped mid-chain atom now raises instead of passing a first/last-atom
  spot check.

### Documentation

- Rebuild the CLI references and skills for the current ML/MM input contract,
  validate required parm7 and XYZ topology references in runnable skill examples,
  and document JSON 2.0 truth/provenance plus the Colab workflow.
- Clarify that recognized monatomic ions use the internal charge table and that
  `-l` entries such as `MG:3` are ignored rather than overriding that table.

## [0.3.2] — 2026-07-10

### Added

- **`--dry-run` on `scan2d` and `scan3d`.** Validates options and prints the execution
  plan without running the scan, pairing with `scan` and pdb2reaction.
- **`--workers` / `--workers-per-node` on every MLIP subcommand** (`sp`, `opt`, `tsopt`, `freq`,
  `irc`, `scan` / `scan2d` / `scan3d`, `path-opt`, `path-search`, `all`), pairing with pdb2reaction.
  `--workers > 1` routes the UMA backend through fairchem's `ParallelMLIPPredictUnit` (needs
  `fairchem-core[extras]`); the parallel predictor exposes no autograd model, so analytical Hessians
  are unavailable and an `Analytical` request is rejected with an error. The default
  `--workers 1` keeps the in-process predictor and is byte-for-byte the previous behavior.
- **Microiteration now works with every Hessian TS optimizer.** `--microiter` (default on)
  previously engaged only with RS-I-RFO (`--opt-mode hess` / `rsirfo`); it now also drives the
  RS-P-RFO (`--opt-mode rsprfo`) and TRIM (`--opt-mode trim`) macro step, alternating a 1-step
  macro TS move with MM-only L-BFGS relaxation. All three are `TSHessianOptimizer` subclasses that
  share the `optimize()`/`prepare_opt()`/Bofill-update contract the macro loop drives. The default
  TS optimizer (RS-I-RFO) is unchanged.

### Changed

- **`--precision` now defaults per backend instead of globally to fp32: ORB runs fp64
  when no precision is given (MACE already did), UMA keeps fp32.** ORB's fp32 is the
  reduced `float32-high` (TF32) matmul mode. Pass `--precision fp32`
  explicitly to select the reduced-precision screening configuration.
- **Behavior change (default): the default UMA model is now `uma-s-1p2`** (was
  `uma-s-1p1`). Other models (`uma-s-1p1`,
  `uma-m-1p1`, MACE-OMOL, Orb-v3-omol) remain selectable via `-b` / `--backend-model` / config.
- **Centralized the default UMA model** in a single constant `DEFAULT_UMA_MODEL`
  (`mlmm/core/defaults.py`); the `uma-s-1p1` defaults previously hardcoded across backends and
  `io/trj2fig.py` now all reference it. CLI help and docs were updated for consistency.

### Fixed

- **MACE backend could not load its own default model.** The `MACE-OMOL-0` default was
  routed to `mace_off()`, which treats any non-preset, non-URL string as a local file
  path (raising `FileNotFoundError`). It now uses the dedicated `mace_omol` factory.
- **An analytical-Hessian request the backend cannot honour now raises instead of
  degrading silently.** A backend build without an autograd model rejects
  `--hessian-calc-mode analytical` with an error naming `FiniteDifference`, rather than
  quietly changing the Hessian method.
- **`make_is_param_explicit` logs on a failed parameter-source query.** Behaviour is
  unchanged (still treats the param as not explicit); the debug log makes a genuine
  typo diagnosable, matching pdb2reaction's `cli_param_overridden`.
- **A `--config` YAML `calc.precision` is now dispatched on every subcommand.** The
  per-subcommand call sites guarded the dispatch helper on the CLI value being
  non-`None`, so a YAML-only precision was never applied (ORB silently ran at its
  default precision); `all` was already correct. The guard is dropped on all 9
  subcommands, re-enabling the helper's YAML path and its invalid-token / aimnet2-fp64
  validators.
- **`--backend-model` supplied via a `--config` YAML is now honoured** on every
  subcommand (same guard as above; only `all` was correct).
- **The deterministic-compute comparison now reuses one MM topology.** Topology
  generation is outside `--deterministic`; pass a fixed `--parm` when comparing
  exact artifacts.
- **Run summary recorded the default UMA model and `mlip_backend: "unknown"`, ignoring `--backend-model` / `-b`.**
  The `all` workflow's summary payload never populated `uma_model`, so `summary.py` fell back to the default
  `MLMM_CALC_KW` model (now `uma-s-1p2`): a run launched with e.g. `--backend-model uma-s-1p1` recorded
  `UMA model: uma-s-1p2` in `summary.log` and `"uma_model": null` / `"mlip_backend": "unknown"` in `summary.json`.
  The computation honored `--backend-model`; the summary now records the
  resolved model (`--backend-model` or `DEFAULT_UMA_MODEL`) and backend (`-b`
  or the `uma` default).

## [0.3.1] — 2026-07-05

### Changed

- Documentation corrections to match the code (`--precision` per-backend default; `pdbfixer`
  is optional; backend/architecture notes).

### Fixed

- **Charge/spin were silently dropped on the UMA MLIP path.** `AtomicData.from_ase` was
  called without `r_data_keys`, so fairchem-core ≥2.x ran the UMA/ORB/MACE backends at
  charge=0/spin=0 regardless of the requested charge/multiplicity. Now passes
  `r_data_keys=["spin","charge"]`, and `atoms.info["spin"]` is the spin multiplicity (2S+1)
  for the OMol backends — it previously used the unpaired-electron count (`mult-1`), which sent
  closed-shell singlets to the null spin token.
- `mm-parm` used a local copy of the residue/ion charge tables that had drifted from
  `mlmm.core.residue_data` (it lacked the phosphorylated S1P/T1P/Y1P and several His tautomers,
  assigning them charge 0, and used a different disulfide cutoff). It now imports the canonical
  tables.

## [0.3.0] — 2026-06-28

### Added

- `--calc-file PATH` (with `--calc-factory NAME`): drive the ML region with an
  arbitrary ASE Calculator loaded from a user Python file (a `custom` backend) —
  usable on every subcommand and forwarded through the `all` pipeline. Couple
  GFN-xTB, DFTB+, ORCA, or any ASE-compatible engine without modifying the
  package; the MM side and ONIOM coupling are unchanged and Hessians use the
  finite-difference path. See `docs/backends.md`.
- `--backend-model NAME` flag on every backend-using subcommand (`opt`,
  `tsopt`, `freq`, `irc`, `scan` / `scan2d` / `scan3d`, `path-opt`,
  `path-search`, `sp`, `all`) to override the model variant for the selected
  `--backend` (e.g. `--backend uma --backend-model uma-m-1p1`), routed to the
  backend's model kwarg via `apply_backend_model_to_calc_cfg`. Previously the
  model variant was settable only through `--config` YAML.
- `--deterministic` flag on every compute subcommand (`opt`, `tsopt`,
  `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`,
  `all`, `sp`) requests deterministic algorithms and an `index_reduce_` shim.
  `MLMM_STRICT_DETERMINISTIC=1` is the environment-variable equivalent;
  backend/model/SDK and target-stack repeatability must be verified.
- `docs/reproducibility.md` documents backend support and verification guidance.
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
    - `AllContext` (frozen dataclass) bundling the 72 `mlmm all` CLI
      parameters in declaration order.
    - `copy_path_outputs_to_root` / `promote_diag_for_root` (replace
      nested closures).
    - `build_energy_level_dict` factors the 5-way R/TS/P energy-level
      dict pattern (UMA / Gibbs / DFT / Gibbs-DFT) into one helper.
    - `build_pipeline_summary_payload` factors the summary-log payload
      assembly so the dict construction is unit-testable separately
      from the I/O wrapper.
    - `build_tsopt_overrides` / `build_freq_overrides` /
      `build_dft_overrides` factor the inline override-dict assembly
      (`if x is not None: cfg[k] = ...` ladder) for the post-MEP
      TSOPT / FREQ / DFT stage calls.
    - `append_backend_forwarding_args` consolidates the 8 copies of
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

### Changed

- **Behavior change (default):** `all --refine-path/--no-refine-path` now
  defaults to `--no-refine-path`. The `all` pipeline's MEP stage runs a
  single-pass `path-opt` by default; pass `--refine-path` to run the recursive
  `path-search` (automatic multi-step bond-change segmentation), which was the
  previous default. The default MEP work directory is now `_work/path_opt/`
  (was `_work/path_search/`). The standalone `path-search` subcommand is
  unchanged. Docs/skills updated throughout.
- `--dft-func-basis` is now surfaced in the primary `mlmm <subcmd> --help`
  (previously only under `--help-advanced`), so the DFT//MLIP/MM functional/basis
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
- The `all`-pipeline determinism comparison uses a fixed smoke stack and input;
  exact repeatability for other stacks requires separate verification.
- `--help` of `mlmm` now groups subcommands under semantic sections
  ("Pipelines" / "Pipeline stages" / "Inputs & topology" / "Analysis")
  in a configurable, deterministic order; subcommands not listed in any
  section fall through to a trailing "Other" bucket so we never hide an
  entry silently.
- CLI exception renderer appends `Try 'mlmm <subcmd> -h' for help.` to
  every user-input-style error so first-time users see a recovery path,
  and routes the full traceback through `logging.getLogger(...).exception`
  so log scrapers / `-v` users get the structured record alongside the
  terminal echo.
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
  (and their `add_*_option` factories). The vendored
  pysisyphus `HessianOptimizer` kwargs are left dormant; no
  behaviour change since defaults were always legacy.

### Fixed

- Bond-change detection (`domain/bond_changes.compare_structures`) is now
  row-chunked instead of building dense N×N distance matrices, removing a CUDA
  out-of-memory failure on large solvated clusters (~20k+ atoms) during
  `path-search` / `scan` kink detection on 16–24 GB GPUs.
- mlmm-toolkit's default MM backend now declares `ninja` as a dependency. Its C++ kernels are
  JIT-compiled through `torch.utils.cpp_extension`, which on modern torch requires
  Ninja and no longer has a distutils fallback, so a clean-env install (e.g. a
  fresh HPC conda env) could fail at runtime with "native bonded extension is
  unavailable". Ninja ships as a pip wheel on every platform (incl. linux-aarch64),
  so the extensions now build out of the box; the dead `use_ninja=False` fallback
  was removed.
- The `hessian_ff` native extensions now build into a local-filesystem
  directory by default (system temp dir, or `TORCH_EXTENSIONS_DIR`), since
  `torch`'s cpp_extension build lock deadlocks on network filesystems and
  previously hung the first build on NFS/Lustre. Requires GCC ≥ 9
  (`conda install -c conda-forge gxx_linux-64`).
- Dropped `torch_geometric` / `torch_scatter` from the package dependencies:
  neither is imported by mlmm and current `fairchem-core` does not require them,
  while `torch_scatter` (sdist only) broke a clean `pip install mlmm-toolkit`
  under PEP 517 build isolation.
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

### Documentation

- Documented that the `[orb]` extra's `torch_scatter` has no PyPI binary wheel
  (sdist only, fails under PEP517 build isolation): install from PyG's
  prebuilt-wheel index, e.g.
  `pip install "mlmm-toolkit[orb]" -f https://data.pyg.org/whl/torch-2.8.0+cu129.html`.
- Documented automatic vs manual ML-region definition (concepts) and clarified
  that `extract --add-linkh` is distance-based and unnecessary when preparing a
  `--model-pdb` for mlmm (the ML/MM calculator caps the boundary from the `--parm`
  topology). Genericized the `json-output` environment example to placeholders,
  and switched the README overview image to an absolute URL so it renders on PyPI.
- Documented the default MM backend's 4-point-water (OPC/TIP4P) limitation and the
  3-point (OPC3/TIP3P) / `--mm-backend openmm` workarounds in `mm-parm.md`.

## [0.2.9] — Unreleased

Consolidated changes since v0.2.4, pending the `v0.2.9` tag. Highlights: end-to-end
JSON output, MM-only optimization mode, GPU memory headroom for 16 GB
consumer cards via Hessian-update CPU offload, 1-based atom indices in
all user-facing I/O, energy-plateau convergence fallback, and a
comprehensive documentation overhaul (EN/JA).

### Added

#### CLI features

- `mlmm bond-summary`: detect bond changes between two structures
  (positional args supported, R/P sanity check companion to `mlmm all`).
- `mlmm opt --mm-only`: skip the MLIP component and minimize on the MM
  force field only. Layers are still honored; `--opt-mode hess` rejected
  (MM-only calculator does not provide a Hessian); microiteration is
  auto-disabled. Useful as a cheap MM pre-relaxation before ML/MM
  ONIOM optimization.
- `--out-json` across all MLIP subcommands (`opt`, `tsopt`, `freq`,
  `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`): emits a
  `result.json` record per run, including backend, charge,
  spin, and timing. `mlmm all` migrates `summary.yaml` → `summary.json`.
- `--link-atom-method scaled|fixed` on all computation subcommands.
- `--cmap / --no-cmap` to control one CMAP policy across both MM layers.
- `--hess-device`, `--read-hess`, `--dump-hess`, `--skip-final-freq`
  for explicit Hessian device control and serialization.
- `--engine gpu|cpu` for `mlmm dft`; `--lowmem` selects
  `gpu4pyscf.rks_lowmem.RKS` as the closed-shell default.
- `--modified-residue` for `mlmm extract` / `mlmm all` (non-standard
  amino acid handling).
- `--freeze-atoms` for `mlmm irc`.
- ML-region charge/multiplicity sanity validator that runs before the
  first MLIP evaluation, surfacing wavefunction-charge mismatches early.

#### Convergence and stability

- Energy-plateau convergence fallback for optimizers (range-based,
  `1e-4 au` threshold; skipped for chain-of-states optimizers).
- Pre-MEP global alignment of reactant / product structures.
- Auto ECP for def2 basis sets in `mlmm dft`.

#### Documentation & infrastructure

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
- `--refine-path` introduced (default `True`, recursive `path-search`; the
  default was later flipped to `--no-refine-path` / single-pass `path-opt` in 0.3.0).
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

### Fixed

#### GPU memory

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

#### Correctness

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

#### Documentation hygiene

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

### Breaking changes

- **Package name**: `mlmm-toolkit` on PyPI. Install with `pip install mlmm-toolkit`.
- **CLI completely redesigned**: The old entry points (`mlmm`, `def_ml_region`, `bond_scan`, `ts_search`, `energy_summary`, `trj2fig`, `add_elem_info`, `get_freeze_indices`, `xyz_geom2pdb`) are replaced by a single `mlmm <subcommand>` interface with 21 subcommands.
- **OpenMM removed**: MM calculations now use `hessian_ff`, a bundled C++ native extension for Amber force fields. OpenMM and OpenMM-CUDA-12 are no longer dependencies.
- **RDKit removed**: No longer required.
- **pysisyphus bundled**: No longer installed from a separate git repository; a modified fork is included in the package.
- **fairchem-core from PyPI**: No longer installed from a custom git fork.
- **numpy constraint relaxed**: `numpy<2.0` → `numpy>=1.24` (NumPy 2.x compatible).
- **Configuration**: Shared runtime defaults were centralized in `defaults.py`; command-local CLI defaults remain with their Click options.

### Added

#### CLI & Workflow

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

#### Core Features

- **Multi-backend MLIP support**: UMA (default), ORB, MACE, AIMNet2 — selectable via `-b/--backend` on all subcommands.
- **ONIOM-like ML/MM decomposition**: `E = E_MM_real + E_ML_model − E_MM_model` with link-atom Jacobian transformation.
- **hessian_ff**: Bundled C++ native extension for analytical Amber force field energies, forces, and Hessians (replaces OpenMM).
- **Microiteration scheme**: Efficient optimization of large ML/MM systems (~10,000 atoms) by separating ML and MM degrees of freedom.
- **xTB embed-charge correction**: Optional point-charge embedding for the ML region (`--embedcharge`).
- **Partial Hessian approach**: Compute Hessians only for the ML region + boundary atoms, enabling TS optimization and frequency analysis on large systems.
- **Automatic versioning**: setuptools-scm replaces hard-coded `__version__`.
- **Progressive help**: `--help` shows primary options; `--help-advanced` shows the full option set.
- **DefaultGroup**: Lazy-loading Click subcommand architecture with automatic boolean normalization (`--flag/--no-flag` and `--flag True/False` both supported).

#### Dependencies & Extras

- **Optional backends**: `pip install "mlmm-toolkit[orb]"`, `"mlmm-toolkit[aimnet]"`. MACE (`mace-torch`) conflicts with `fairchem-core` (UMA) due to incompatible `e3nn` versions; use separate conda environments if both are needed.
- **Optional DFT**: `pip install "mlmm-toolkit[dft]"` for PySCF + gpu4pyscf-cuda12x + CuPy.
- **Optional PDBFixer**: `pip install "mlmm-toolkit[pdbfixer]"` for hydrogen addition.
- **New core dependencies**: `click`, `torch_geometric`, `torch_scatter`, `pyparsing`, `tabulate`.
- **Bundled packages**: `pysisyphus` (modified fork), `thermoanalysis`, `hessian_ff`.

#### Documentation & Testing

- Bilingual documentation (English + Japanese) under `docs/` and `docs/ja/`.
- Per-subcommand reference docs with CLI tables and YAML configuration examples.
- `CONTRIBUTING.md` with contributor guidelines.
- GitHub Actions CI: `pytest.yml`, `smoke_test.yml`, `docs_quality.yml`.
- 198 unit tests with `pytest --timeout=120`.
- Smoke test suite: 34 end-to-end tests covering all subcommands.
- Working examples in `examples/toy_system/` and `examples/methyltransferase/`.

### Changed

#### Architecture — v0.1.1 → v0.2.0

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

### Fixed

#### (relative to v0.1.1)

- IRC energy-based initial displacement: Added `step_length` clamp (max 0.5 au) to prevent divergence when `min_eigval ≈ 0`.
- Hessian cache: Fixed `set_calculator()` clearing `within_partial_hessian` before `cart_hessian` assignment.
- HessianOptimizer: Fixed GPU/CPU device mismatch in `hessian_recalc`.
- Frequency analysis: Added float64 enforcement and explicit `(H+Hᵀ)/2` symmetrization.
- PDB element inference: Fixed `_infer_element_from_pdb_atom_name()` misidentifying protein hydrogens (e.g., `HG2` → `Hg`). Now uses residue-context-aware `guess_element()`.
- Silent `except Exception: pass` blocks converted to `logger.debug(...)` across all modules.
- Charge derivation: Uses `--model-pdb` instead of full input when provided.
- TypeError in `_pseudo_irc_and_match` when mapping value is `None`.
