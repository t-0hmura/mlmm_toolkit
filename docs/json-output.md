# JSON Output Reference

mlmm provides machine-readable JSON output for programmatic consumption by AI agents, scripts, and downstream tools.

## `--out-json` flag

Most MLIP-based and reporting subcommands (`opt`, `sp`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`, `extract`, `trj2fig`, and `energy-diagram`) support `--out-json / --no-out-json` (default: off).
When enabled, authoritative `result.json` and its identical `summary.json` compatibility mirror are written beside the normal outputs.

```bash
mlmm opt -i r_complex_layered.pdb --parm real.parm7 -q 0 -m 1 \
  --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

The `all` and `path-search` commands write `summary.json` without an `--out-json` flag once execution reaches their summary writer. Early CLI or input validation can fail before the file exists.

### `summary.json` mirror

`write_result_json` stages the same bytes for both names, publishes the
`summary.json` compatibility mirror first, and publishes authoritative
`result.json` last. A successful return guarantees identical bytes. If
publication is interrupted, treat `result.json` as authoritative and require
successful process/writer completion; when an orchestrator assigned `run_id`,
validate it as well rather than assuming both names identify one generation.

## Common envelope

The shared writer and aggregate summary producers supply the fields below.
Rows marked optional are present only when the producer supplies that data:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Envelope schema version; current value comes from `mlmm.core.utils.RESULT_JSON_SCHEMA_VERSION` — pin against that constant rather than the literal in this doc. Bumps signal a structural change. |
| `command` | string | Leaf envelopes use the subcommand name (e.g. `"opt"`); aggregate `all` / `path-search` summaries record the full invocation. |
| `mlmm_version` / `mlmm_toolkit_version` | string | Package version (`mlmm_version` in leaf envelopes; `mlmm_toolkit_version` in aggregate summaries). |
| `status` | string | Command-specific: `all` uses `success`/`partial`/`failed`; `path-search` uses `success`/`partial`; `opt` and `tsopt` use the numerical outcomes `converged`/`not_converged`/`stalled`; completed analysis/integration stages use `completed`; exception envelopes use `error`. TS saddle order is reported separately in `saddle_validation` / `hessian_status`. |
| `elapsed_seconds` | float | Optional wall-clock time; omitted when the producer does not pass timing to the shared writer. |
| `environment` | object | Hardware info (see below) |
| `run_id` | string | Optional. Present when an orchestrator (including MCP) assigns a current invocation identity; conflicting caller values are rejected. |

MLIP/ML/MM calculator stages additionally record:

| Field | Type | Description |
|-------|------|-------------|
| `mlip_backend` | string \| null | Backend identifier (`uma`, `orb`, `mace`, `aimnet2`, `dft`, or `custom`); DFT leaves emit `dft`, and plot-only commands that did not evaluate a calculator emit null |
| `mlip_model` | string \| null | Exact model/checkpoint; `filename:factory` for `--calc-file` |
| `mlip_precision` | string \| null | Effective public precision (`fp32` or `fp64`); null for custom calculators |
| `mm_backend` | string \| null | MM Hessian/energy backend (`hessian_ff` or `openmm`); null when a plot-only command did not evaluate a calculator |
| `link_atom_method` | string \| null | Link-atom placement (`scaled` or `fixed`); null for plot-only output |
| `use_cmap` | bool \| null | Whether CMAP terms were enabled; null for plot-only output |

### Execution and scientific truth

Multi-stage and scan producers add the fields below when they can evaluate constituent work. These fields are additive and producer-dependent; the command-specific `status` remains in place. Consumers should gate scientific use on `scientific_status` and the leaf outcomes. Missing required acceptance signals are fail-closed. IRC endpoint stationarity is diagnostic rather than an acceptance signal; IRC usability is reported separately from propagation validity.

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` | string | Normally `completed` or `failed`; reports whether required constituent commands executed. |
| `scientific_status` | string | `success`, `partial`, or `failed`; reports whether the produced scientific result is complete and usable. |
| `scientific_status_reasons` | string[] | Reasons for unusable or missing leaves; omitted on clean success. This is distinct from an aggregate workflow's legacy `status_reasons`. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed leaf identifiers used to detect missing aggregate work. |
| `stage_outcomes` | object[] | Stage leaves with `stage`, `item_id`, `required`, `executed`, `converged`, `usable`, `reason`, and `artifacts`. |
| `point_outcomes` | object[] | Scan points with `point_id`, `executed`, `converged`, `energy_valid`, `artifact_written`, `seed_eligible`, and `reason`. |

When present, `run_id` identifies the current invocation. The `all` aggregate
rebuilds `current_output_paths` and `key_output_files` from that invocation's
manifest, so stale files in a reused output tree are excluded.

### Error envelope (when `status == "error"`)

| Field | Type | Description |
|-------|------|-------------|
| `error` | string | `str(exc)` of the original exception |
| `error_type` | string | Exception class name (e.g. `"OptimizationError"`) |
| `error_class_chain` | list[string] | Full MRO class names (e.g. `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`) so agents can match the hierarchy without parsing text |
| `error_module` | string | Module the exception class was defined in |
| `error_label` | string | High-level CLI stage label (e.g. `"opt"`, `"tsopt-stage"`) |

**`environment`**:

| Field | Type | Example |
|-------|------|---------|
| `device` | string | `"cuda"` or `"cpu"` |
| `gpu_name` | string | `"<gpu model>"` |
| `gpu_vram_gb` | float | `<vram in GB>` |
| `cuda_version` | string | `"<cuda version>"` |
| `cpu` | string | `"<cpu model>"` |
| `n_cpus` | int | `<int>` |
| `ram_gb` | float | `<ram in GB>` |

An optimizer may also report `"status": "stalled"`: the energy stopped decreasing over the configured window (an energy plateau) while the configured force/step convergence criteria remained unmet. A stall is a distinct, non-converged outcome — it is never reported as `converged` — and it stops further flatten/retry work rather than repeating a non-progressing optimization. When present, a `stop_reason` string records the energy range, window, and the failed criteria. A stall may be retried (e.g. from a perturbed geometry or with tighter step control); it is not an alias for `max_cycles` exhaustion or a generic failure. In microiteration, a stalled macro step or latest micro (MM) relaxation remains a stalled result and cannot satisfy macro convergence.

## Subcommand schemas

### `sp`

| Field | Type | Description |
|-------|------|-------------|
| `status` / `stage` | string / string | `"ok"` / `"sp"` |
| `input` | string | Input structure path |
| `real_parm7` | string | Full-system Amber topology path |
| `charge` / `spin` | int / int | Model-region charge and multiplicity |
| `energy_au` | float | ONIOM single-point energy (Hartree) |
| `forces_path` | string | Path to `forces.npy` |
| `hessian_path` | string \| null | Path to `hessian.npy`, or null without `--hess` |
| `elapsed` | string | Human-readable elapsed time |

### `opt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"converged"`, `"not_converged"`, or `"stalled"` (energy plateau; see above) |
| `stop_reason` | string | Present only for a non-converged stop (stalled/stopped); records the energy plateau range/window and failed criteria |
| `energy_hartree` | float | Final ONIOM energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"lbfgs"`, or `"rfo"` |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `n_atoms` | int | Total atoms (all layers) |
| `n_freeze_atoms` | int | Frozen atoms |
| `thresh` | string | Convergence threshold preset |
| `max_cycles` | int | Maximum allowed cycles |
| `input_file` | string | Input filename |
| `final_max_force` | float | Last max gradient (Hartree/Bohr) |
| `final_rms_force` | float | Last RMS gradient |
| `final_max_step` | float | Last max displacement (Bohr) |
| `final_rms_step` | float | Last RMS displacement |
| `convergence_thresholds` | object | Numeric thresholds for the named preset |
| `rigid_projection` | object\|null | Present when `--flatten` performs PHVA; frozen-boundary TR provenance |
| `files` | object | Output file map |

### `tsopt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | Backward-compatible numerical outcome; use `optimization_status` and `saddle_validation` separately |
| `optimization_status` | string | Numerical optimizer outcome: `"converged"`, `"not_converged"`, or `"stalled"`; independent of saddle order |
| `saddle_validation` | string | `"first_order"`, `"higher_order"`, `"no_imaginary"`, or `"unavailable"` from terminal exact PHVA |
| `saddle_order_verified` | bool | `true` only for `saddle_validation: "first_order"` |
| `hessian_status` | string | `"completed"`, `"failed"`, `"skipped"`, or `"unavailable"`; `hessian_error` gives the failure reason |
| `reaction_mode_index` | int\|null | Selected negative exact-PHVA root for downstream IRC; fallback root 0 is explicitly labelled and does not verify reaction identity |
| `reaction_mode_frequency_cm` | float\|null | Frequency of the selected negative root |
| `reaction_mode_source` | string\|null | Reference-aligned or explicit fallback source used for root selection |
| `energy_hartree` | float | TS energy (Hartree) |
| `n_imaginary_modes` | int\|null | Number of imaginary frequencies; `null` if PHVA was not run |
| `imaginary_frequencies_cm` | float[]\|null | Imaginary frequencies (cm⁻¹, negative); no PHVA: `[]` with `--skip-final-freq`, otherwise `null` |
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"dimer"`, `"rsprfo"`, `"rsirfo"`, or `"trim"`; `hess` selects RS-P-RFO |
| `opt_mode_requested` | string | Requested CLI preset |
| `optimizer` | string | Effective optimizer algorithm used by the run |
| `n_atoms` | int | Total atoms |
| `n_opt_cycles` | int | Optimization cycles |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `reference_mode_file` | string\|null | Advanced path-derived mode supplied with `--ref-mode`; Hessian-family only |
| `safeguards` | object | Hessian-family rejection/recovery, exact-saddle, and target-mode diagnostics |
| `rigid_projection` | object | Frozen-boundary TR provenance for Dimer/flatten/final saddle analysis |
| `files` | object | Final geometry + vib mode files |

Terminal exact PHVA runs only after numerical convergence. A non-converged or
stalled run retains the terminal geometry and records PHVA as skipped. A PHVA failure is recorded as
`hessian_status: "failed"` without discarding the structure or fabricating
frequencies. Numerical status and saddle order are separate: a converged
higher-order stationary point remains `optimization_status: "converged"` with
`saddle_validation: "higher_order"`, and is not a certified first-order TS.
`all` may continue warning-labelled diagnostic IRC only with a validated
negative root. Numerical non-convergence, zero imaginary modes, failed/skipped
PHVA, or no valid negative root stops after TS artifact registration and before
IRC. Explicit `--skip-final-freq` retains the final structure with
`n_imaginary_modes: null` and `imaginary_frequencies_cm: []`.

### `freq`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_modes` | int | Total normal modes |
| `n_imaginary` | int | Imaginary frequency count |
| `frequencies_cm` | float[] | All frequencies (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | Negative frequencies only |
| `thermochemistry` | object\|null | Thermodynamic data (see below) |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `n_atoms` | int | Total atoms |
| `n_freeze_atoms` | int | Frozen atoms |
| `rigid_projection` | object | Frozen-boundary TR provenance used for frequencies and thermochemistry |
| `files` | object | Output map; includes `hessian_npz` when `--dump-hess` is used |

**`thermochemistry`** (null if thermoanalysis unavailable):

| Field | Type | Unit |
|-------|------|------|
| `temperature_K` | float | K |
| `pressure_atm` | float | atm |
| `point_group` | string | Automatically detected molecular point group |
| `point_group_source` | string | `"auto"` or conservative `"auto-fallback"` |
| `symmetry_number` | int | External rotational symmetry number |
| `symmetry_number_source` | string | `"auto"`, `"auto-fallback"`, `"config"`, or `"override"` |
| `electronic_energy_ha` | float | Hartree — the `E` of the reported `E + G_corr = G` |
| `zpe_ha` | float | Hartree |
| `thermal_correction_energy_ha` | float | Hartree |
| `thermal_correction_enthalpy_ha` | float | Hartree |
| `thermal_correction_free_energy_ha` | float | Hartree |
| `sum_EE_and_ZPE_ha` | float | Hartree |
| `sum_EE_and_thermal_energy_ha` | float | Hartree |
| `sum_EE_and_thermal_free_energy_ha` | float | Hartree |
| `E_thermal_cal_per_mol` | float | cal/mol |
| `Cv_cal_per_mol_K` | float | cal/(mol K) |
| `S_cal_per_mol_K` | float | cal/(mol K) |

### `irc`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_frames_forward` / `n_frames_backward` / `n_frames_total` | int | IRC frames |
| `energy_first_hartree` | float | First stitched-path endpoint; standalone IRC assigns no chemical identity |
| `energy_ts_hartree` | float | TS energy |
| `energy_last_hartree` | float | Last stitched-path endpoint; standalone IRC assigns no chemical identity |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `energy_reactant_hartree` / `energy_product_hartree` | float | Compatibility aliases for first/last; do not infer R/P identity from the names |
| `forward_requested` / `backward_requested` | bool | Whether each direction was requested |
| `forward_status` / `backward_status` | string | `stopped`, `failed`, or `disabled`; use these for directional propagation status |
| `forward_endpoint_stationary` / `backward_endpoint_stationary` | bool\|null | Whether the raw endpoint met the stationary-point threshold; diagnostic only |
| `forward_converged` / `backward_converged` | bool\|null | Compatibility aliases for `*_endpoint_stationary`; not the IRC usability gate |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool\|null | Whether the branch established a downhill departure from the TS |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string\|null | Non-empty only for a numerical propagation failure |
| `never_stop` | bool | Whether opt-in physical endpoint-stop bypass mode was enabled |
| `never_stop_energy_bypasses` | int | Number of energy-rise or one-step energy-change stops actually bypassed |
| `rigid_projection` | object | Frozen-boundary TR provenance for the initial/updated Hessian |
| `rigid_projection.electronic_state_verified` | bool | For a file-seeded Hessian, whether model charge and multiplicity were identity-verified |
| `bond_changes` | object | Directed first→last `{formed: [...], broken: [...]}`; omitted if comparison was unavailable |
| `bond_changes_direction` | string | `"finished_first_to_finished_last"` when bond changes are present |
| `files` | object | Trajectory and endpoint files (XYZ plus available PDB/CIF companions) |

**`rigid_projection` provenance:** the object records the selected treatment
(`treatment`), `effective_rank`, active/frozen atom counts and indices, and the
Hessian source/shape used by that workflow. The treatment is always
`constrained`; a stale non-constrained configuration fails explicitly. A `freq --dump`
run writes the same object to `thermoanalysis.yaml`. Field names for the final
two values follow the producing workflow (`hessian_source` / `hessian_shape`,
or `source` / `raw_hessian_shape`).

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `scan_opt_mode` | string | Fixed `grad` preset used by the L-BFGS constrained relaxations |
| `scan_optimizer` | string | Effective optimizer identity (`lbfgs`) |
| `n_stages` | int | Number of scan stages |
| `stages` | object[] | Per-stage data |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `files` | object | Output files |

**`stages[]`**: `n_steps`, `converged`, `pairs_1based`, `energies_hartree`, `final_energy_hartree`, `bond_changes`, and (additive) `optimizer_status` (`converged`/`not_converged`/`stalled`) plus `stop_reason` when the stage's last optimizer stopped without convergence

### `scan2d` / `scan3d`

| Field | Type | Description |
|-------|------|-------------|
| `n_grid_points` | int | Total grid points |
| `n_points_attempted` | int | Fresh-run grid points attempted (preoptimization row excluded) |
| `n_points_usable` | int | Fresh-run points with explicit convergence, finite energy/coordinates, and a written geometry artifact |
| `grid_points` | object[] | Explicit grid-index, distances, energy, convergence, and `geometry_file` mapping used by interactive Results viewers |
| `current_output_paths` | string[] | Current-run CSV/HTML/PNG and grid geometries; consumers can exclude stale files without filename inference |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | Surface minimum energy |
| `charge` | int \| null | Model-region charge; null for plot-only `scan3d --csv` |
| `spin` | int \| null | Model-region multiplicity; null for plot-only `scan3d --csv` |
| `files` | object | CSV + plot files |

Fresh `scan2d`/`scan3d` results include the common MLIP/ML/MM calculator
provenance. Plot-only `scan3d --csv` keeps the same keys but writes null because
the imported energy grid does not identify the calculator that produced it.
Plot-only results omit `n_points_attempted`; they include `n_points_usable`
only when the imported CSV has complete convergence and artifact provenance.

### `path-opt`

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool | Convergence flag |
| `mep_mode` | string | `"dmf"` or `"gsm"` |
| `image_energies_hartree` | float[] | All image energies |
| `n_images` | int | Image count |
| `hei_index` | int | Highest-energy image index |
| `barrier_kcal` | float | Forward barrier (kcal/mol) |
| `delta_kcal` | float | Reaction energy (kcal/mol) |
| `files` | object | Trajectory + HEI files |

### `dft`

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool | SCF converged? |
| `status` | string | `"converged"` or `"not_converged"`; the latter is committed before exit code 3. |
| `energy_hartree` | float | DFT energy |
| `xc_functional` | string | XC functional |
| `basis_set` | string | Basis set |
| `used_gpu` | bool | GPU acceleration used? |
| `charges` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `spin_densities` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `n_atoms` | int | QM-region atom count |
| `grid_level` | int | DFT grid level |
| `conv_tol` | float | SCF convergence tolerance |
| `max_cycle` | int | Effective maximum SCF iterations after YAML/CLI resolution |
| `engine` | string | Actual runtime engine label (`pyscf(cpu)`, `gpu4pyscf`, or low-memory GPU variant) |
| `files` | object | `{"result_yaml": "result.yaml"}` |

### `trj2fig`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_frames` | int | Number of trajectory frames. |
| `min_energy_hartree` / `max_energy_hartree` | float | Minimum and maximum frame energies. |
| `energy_source` | string | `"trajectory_comment"` or `"mlip_recomputed"`. |
| `mlip_backend` / `mlip_model` / `mlip_precision` | string \| null | Resolved recomputation provenance; all are null in comment mode. |
| `charge` / `multiplicity` | int \| null | Resolved recomputation state; null in comment mode. Omitted recomputation values resolve to 0 and 1. |
| `output_files` | string[] | Canonical ordered paths for every output; preserves files with the same basename in different directories. |
| `files` | object | Legacy basename-to-path map; retained for compatibility and therefore lossy when basenames collide. |

Supplying either `-q/--charge` or `-m/--multiplicity` recomputes every frame with the selected MLIP. This is a direct MLIP frame rescore; the command has no topology or model-region input and does not calculate an ONIOM energy.

### `extract`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_atoms_extracted` | int | Atoms after extraction |
| `total_charge` | float | Computed total charge |
| `protein_charge` | float | Protein charge |
| `ligand_total_charge` | float | Ligand charge sum |
| `ion_total_charge` | float | Ion charge sum |
| `unknown_residue_charges` | object | `{resname: charge}` |
| `center` | string | Substrate specification (raw `-c` value): PDB path, residue-ID list (e.g. `'A:123,B:456'`), or residue-name list (e.g. `'GPP,MMT'`) |
| `radius` | float | Extraction radius (angstrom) |
| `input_files` | string[] | Input PDB paths |
| `n_atoms_raw` | int | Atom count in the raw input before extraction |
| `n_link_hydrogens` | int | Count of link H atoms added at severed bonds |
| `files` | object | Map of emitted file names (per-input pocket PDB, etc.) |
| `exclude_backbone` | bool | Value of `--exclude-backbone` at run time |
| `include_h2o` | bool | Value of `--include-h2o` at run time |
| `ligand_charge_input` | string | Raw `-l/--ligand-charge` argument |
| `ion_charges` | array | List of `[resname, charge]` pairs for ion residues encountered |

### `energy-diagram`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"ok"` |
| `n_points` | int | Number of energy data points |
| `files` | object | Output diagram filename-to-path map |

## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json`:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` / `"failed"` for `all`; `"success"` / `"partial"` for `path-search`. |
| `execution_status` / `scientific_status` | string / string | Execution completeness and scientific usability; evaluate these separately from legacy `status`. |
| `scientific_status_reasons` | string[] | Reasons for incomplete or unusable science; omitted on clean success. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed aggregate leaves. |
| `config` | object | Effective settings. `mep_mode` identifies GSM/DMF; `ts_opt_mode` and `endpoint_opt_mode` identify the configured post-processing presets. Generic `opt_mode*` keys retain the resolved CLI inputs. `path_opt_mode` is the single-structure optimizer used for endpoint preoptimization (see `preopt`), not the MEP path algorithm. |
| `n_segments` | int | Segment count |
| `segments` | object[] | Per-segment barrier, delta, bond changes |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Backend name (`uma`, `orb`, `mace`, `aimnet2`, or `custom`) |
| `mlip_model` | string \| null | Exact model/checkpoint name, recorded separately from the backend |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64`; null for custom calculators |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `environment` | object | Hardware info |
| `references` | object[] | Methods actually used by the resolved workflow, as `{method, citation, doi}` records. The same reference set is grouped at the end of `summary.log` and final stdout immediately before elapsed time. |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `n_segments_reactive` | int | Number of non-bridge (reactive) segments |
| `rate_limiting_step` | object | Legacy key for the highest independently referenced local segment barrier. It is not a microkinetic rate-limiting-step assignment. |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `post_segments` | list | Per-segment TS/IRC/freq/DFT results |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | Raw propagation, pre-optimization orientation, and final optimized-endpoint acceptance, respectively. Normal raw stopping is diagnostic; endpoint convergence and optimized connectivity govern MEP-mode acceptance. |
| `post_segments[].thermo_symmetry` | object | Child-reported point-group and rotational-symmetry provenance by state: R/TS/P for MEP runs and E1/TS/E2 for TS-only runs. States with valid symmetry-number provenance are included; missing states are omitted, and the field is absent only when no state has valid provenance. |
| `key_output_files` | object | Current-run output index: root filename → description; each `seg_NN` entry is `{description, files}` with paths relative to that segment directory. |
| `current_output_paths` | string[] | Sorted paths relative to `--out-dir`, limited to artifacts claimed by the current invocation. |

## Usage examples

### Python

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

if result["status"] == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
else:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
```

### jq

```bash
# Check convergence
jq '.status' result.json

# Get barrier from path-opt
jq '.barrier_kcal' result.json

# List imaginary frequencies from tsopt
jq '.imaginary_frequencies_cm' result.json

# Get thermochemistry from freq
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json
```
