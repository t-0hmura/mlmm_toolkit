# JSON Output Reference

This page lists the fields (keys) of the `result.json` and `summary.json` files written with `--out-json`, split into the fields every command shares and the fields of each command.

## `--out-json` flag

Most subcommands that run the ML/MM calculator or write reports support
`--out-json / --no-out-json` (default: off). When enabled, `result.json` and
`summary.json` are written beside the normal outputs. The two files have the
same content, so read `result.json`.

```bash
mlmm opt -i r_complex_layered.pdb --parm7 real.parm7 -q 0 -m 1 \
  --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

Open `result_opt/result.json` and read `execution_status` and `scientific_status` first. `opt`, `tsopt`, and `path-opt` also record `optimization_status`.

## Common envelope

Every result file carries the fields below; fields marked optional appear only
when the command has the corresponding data:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Schema version of the file; a new version signals a structural change. |
| `command` | string | Single commands record the subcommand name (e.g. `"opt"`); the `all` / `path-search` summaries record the full command line |
| `mlmm_version` / `mlmm_toolkit_version` | string | Package version (`mlmm_version` in single-command results; `mlmm_toolkit_version` in the `all` / `path-search` summaries) |
| `execution_status` | string | Execution completion: `completed` / `failed`. |
| `scientific_status` | string | Result usability: `success` / `partial` / `failed`. |
| `run_id` | string | Optional UUID of the current invocation; written when the MCP server starts the command and in every `all` run, including its stages. |
| `elapsed_seconds` | float | Optional wall-clock time; omitted by commands that do not record timing |
| `environment` | object | Hardware info (see below) |

Commands that evaluate the ML/MM calculator also record:

| Field | Type | Description |
|-------|------|-------------|
| `mlip_backend` | string \| null | Backend identifier (`uma`, `orb`, `mace`, `aimnet2`, `dft`, or `custom`); `dft` results record `dft`, and plot-only commands that did not evaluate a calculator record null |
| `mlip_model` | string \| null | Exact model/checkpoint; `filename:factory` for `--calc-file` |
| `mlip_model_label` | string \| null | Publication-facing model label derived from the exact identifier |
| `mlip_task` | string \| null | Backend task used by a multi-domain model (for UMA, `omol` unless `calc.uma_task_name` sets another); `mlip_model` remains the exact identifier |
| `mlip_precision` | string \| null | Effective public precision token (`fp32` or `fp64`); null for a custom calculator whose dtype is controlled by user code |
| `mm_backend` | string \| null | MM Hessian/energy backend (`hessian_ff` or `openmm`); null when a plot-only command did not evaluate a calculator |
| `link_atom_method` | string \| null | Link-atom placement (`scaled` or `fixed`); null for plot-only output |
| `use_cmap` | bool \| null | Whether CMAP terms were enabled; null for plot-only output |

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

### Execution and requested-stage completion

Every result reports `execution_status` and `scientific_status`; multi-stage and scan results also list the outcome of each stage in the fields below. A required optimization or calculation that is missing leaves the result incomplete.

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` | string | `completed` or `failed`. In `all`, an endpoint optimization after IRC that does not converge leaves it `completed`; one that stops on an error makes it `failed`. |
| `scientific_status` | string | `success` when every requested stage converged; otherwise `partial` or `failed`. In `all`, a TS with n_imag ≥ 2 gives `partial`, and a TS with n_imag = 0 stops the run before IRC, so `success` means n_imag = 1; a valid TS plus one failed endpoint optimization is also `partial`. Standalone `tsopt` looks only at convergence, not at n_imag; read `n_imaginary_modes`. |
| `scientific_status_reasons` | string[] | Reasons for unusable or missing stages; omitted on clean success. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed stage identifiers, used to detect missing work. |
| `stage_outcomes` | object[] | One entry per stage with `stage`, `item_id`, `required`, `executed`, `converged`, `usable`, `reason`, and `artifacts`. |
| `point_outcomes` | object[] | Scan points with `point_id`, `executed`, `converged`, `energy_valid`, `artifact_written`, `seed_eligible`, and `reason`. |

### Error envelope (when `execution_status == "failed"`)

| Field | Type | Description |
|-------|------|-------------|
| `error` | string | The error message |
| `error_type` | string | Exception class name (e.g. `"OptimizationError"`) |
| `error_class_chain` | list[string] | The exception class and its parent classes, most specific first (e.g. `["OptimizationError", "RuntimeError", "Exception", "BaseException"]`), so agents can match the hierarchy without parsing text |
| `error_module` | string | Module the exception class was defined in |
| `error_label` | string | High-level CLI stage label (e.g. `"optimization"` for `opt`, `"TS optimization"` for `tsopt`) |

## Error handling

When `opt`, `tsopt`, `freq`, `irc`, `sp`, `dft`, `scan`, `scan2d`, `scan3d`,
`path-opt`, or `path-search` stops on an exception after the output directory
is set up, `result.json` and `summary.json` are written even without
`--out-json`, with `"execution_status": "failed"` and an `"error_type"`; for a
failure before that point, see [Notes](#notes).

When an optimization ends without converging, `result.json` records
`"optimization_status": "not_converged"`. For what to change before a retry,
see {ref}`Troubleshooting › Calculation / convergence <calculation--convergence>`.

An optimizer may also report `"optimization_status": "stalled"`: the energy stopped decreasing over the configured window (an energy plateau) while the force/step convergence criteria remained unmet. A stall is a kind of non-convergence, never `converged`; `stop_reason` records the energy range, the window, and the failed criteria. With `--flatten`, `opt` and `tsopt` still run the flatten loop after a stall while `--max-cycles` cycles remain. With microiteration, the run counts as converged only when both the macro step and the latest MM relaxation have converged; an MM relaxation that stalls with its forces already under the thresholds counts as converged.

## Subcommand schemas

### `sp`

| Field | Type | Description |
|-------|------|-------------|
| `stage` | string | `"sp"` |
| `input` | string | Input structure path |
| `real_parm7` | string | Full-system Amber topology path |
| `charge` / `spin` | int / int | ML-region charge and multiplicity |
| `energy_au` | float | ONIOM single-point energy (Hartree) |
| `forces_path` | string | Path to `forces.npy` |
| `hessian_path` | string \| null | Path to `hessian.npy`, or null without `--hess` |
| `elapsed` | string | Human-readable elapsed-time text |

### `opt`

| Field | Type | Description |
|-------|------|-------------|
| `optimization_status` | string | `"converged"`, `"not_converged"`, or `"stalled"` (energy plateau; see [Error handling](#error-handling)) |
| `stop_reason` | string | Present only when the optimizer stopped without converging (`stalled` or `not_converged`); records why, e.g. the energy-plateau range/window and the failed criteria |
| `energy_hartree` | float | Final ONIOM energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | `"grad"`, `"hess"`, `"lbfgs"`, or `"rfo"` |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `n_atoms` | int | Total atoms (all layers) |
| `n_freeze_atoms` | int | Frozen atoms |
| `thresh` | string | Convergence threshold preset |
| `max_cycles` | int | Maximum allowed cycles |
| `input_file` | string | Input filename |
| `final_max_force` | float | Last max gradient (Hartree/Bohr) |
| `final_rms_force` | float | Last RMS gradient |
| `final_max_step` | float | Last max displacement (Bohr) |
| `final_rms_step` | float | Last RMS displacement |
| `convergence_thresholds` | object | Numeric thresholds of the named preset |
| `files` | object | Output file map |
| `rigid_projection` | object \| null | Optional; present when `--flatten` runs. See [projection provenance](#rigid-projection-provenance). |

### `tsopt`

| Field | Type | Description |
|-------|------|-------------|
| `flatten_requested` / `flatten_enabled` | bool | Whether flatten iterations were configured |
| `flatten_skip_reason` | string \| null | Reason no further flatten step was taken, when applicable |
| `optimization_status` | string | Numerical optimizer outcome: `"converged"`, `"not_converged"`, or `"stalled"`; independent of saddle order |
| `saddle_validation` | string | `"first_order"`, `"higher_order"`, `"no_imaginary"`, or `"unavailable"` from terminal exact PHVA (partial Hessian vibrational analysis) |
| `saddle_order_verified` | bool | `true` only for `saddle_validation: "first_order"` |
| `hessian_status` | string | `"completed"`, `"failed"`, `"skipped"`, or `"unavailable"`; `hessian_error` gives the failure reason |
| `reaction_mode_index` | int\|null | Negative exact-PHVA root that `all` follows in IRC; when no reference-aligned root is available, root 0 is used and labelled as such, which does not confirm that it is the reaction mode |
| `reaction_mode_frequency_cm` | float\|null | Frequency of the selected negative root |
| `reaction_mode_source` | string\|null | How the root was selected (aligned with the reference direction, or the root-0 fallback) |
| `energy_hartree` | float \| null | TS energy (Hartree); `null` when the final energy evaluation failed (every non-finite number is written as `null`), in which case execution and scientific status are `failed` |
| `n_imaginary_modes` | int\|null | Number of imaginary frequencies; `null` if PHVA was not run |
| `n_negative_modes` | int\|null | Number of negative frequencies of any size, including those within the cutoff; a diagnostic next to `n_imaginary_modes`, `null` if PHVA was not run |
| `imaginary_frequencies_cm` | float[]\|null | Imaginary frequencies (cm⁻¹, negative); without PHVA, `[]` with `--skip-final-freq` and otherwise `null` |
| `frequency_zero_cutoff_cm` / `imaginary_mode_criterion` / `imaginary_frequency_threshold_cm` | float / string / float | The defaults are `5.0`, `"frequency_cutoff_cm"`, and `-5.0`: only ν < −5.00 cm⁻¹ counts as imaginary. |
| `opt_mode` | string | `"grad"`, `"hess"`, `"dimer"`, `"rsprfo"`, `"rsirfo"`, or `"trim"`; `hess` selects RS-P-RFO |
| `opt_mode_requested` | string | Requested CLI preset |
| `optimizer` | string | Effective optimizer algorithm used by the run |
| `n_atoms` | int | Total atoms |
| `n_opt_cycles` | int | Optimization cycles |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `reference_mode_file` | string\|null | Path-derived mode file supplied with `--ref-mode`; Hessian-based optimizers only |
| `safeguards` | object | Hessian-TS diagnostics: rejected steps and recovery, exact saddle checks, and target-mode identity |
| `rigid_projection` | object | Rigid-mode and Hessian provenance for Dimer, flatten, and the final saddle analysis; see [projection provenance](#rigid-projection-provenance) |
| `files` | object | Final geometry and vibrational-mode files; includes `hessian_npy` (absolute path) when `--dump-hess` wrote a file |

A successful TS optimization gives one imaginary mode along the reaction
coordinate: `saddle_validation: "first_order"` with `n_imaginary_modes: 1`.
The final PHVA runs after the optimizer converges or stops on an energy plateau
(`stalled`); otherwise it is recorded as skipped, and a failed PHVA gives
`hessian_status: "failed"` and keeps the final geometry.
`optimization_status` and `saddle_validation` are independent, so a converged
run can end with `saddle_validation: "higher_order"`; it is not a first-order
TS. For when `all` goes on to IRC, see
[tsopt › Reading the TS result](tsopt.md#reading-the-ts-result).

### `freq`

| Field | Type | Description |
|-------|------|-------------|
| `n_modes` | int | Total normal modes |
| `n_imaginary` | int | Imaginary frequency count |
| `n_negative_modes` | int | Number of negative frequencies of any size, including those within the cutoff |
| `frequencies_cm` | float[] | All frequencies (cm⁻¹) |
| `imaginary_frequencies_cm` | float[] | Negative frequencies only |
| `thermochemistry` | object\|null | Thermodynamic data (see below) |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `n_atoms` | int | Total atoms |
| `n_freeze_atoms` | int | Frozen atoms |
| `files` | object | Output map; includes `hessian_npy` (absolute path) when `--dump-hess` wrote a file |
| `rigid_projection` | object | Rigid-mode and Hessian provenance used for frequencies and thermochemistry; also written to `thermoanalysis.yaml` with `--dump` |

**`thermochemistry`**:

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

IRC reports `execution_status` and `scientific_status` and keeps the stop reason and trajectory of each direction; `all` reports the endpoint optimizations that follow under `endpoint_opt`. The stitched path runs from its first frame through the TS to its last frame; a standalone `irc` does not decide which end is the reactant and which is the product. How the IRC stopped does not enter `scientific_status`: a standalone `irc` reports `completed` / `success` when the integration returns normally, and `all` judges the TS and endpoint optimizations.

| Field | Type | Description |
|-------|------|-------------|
| `n_frames_forward` / `n_frames_backward` / `n_frames_total` | int | IRC frames |
| `forward_short_branch` / `backward_short_branch` | bool | Branch produced at most three frames without reaching the cycle cap; diagnostic only |
| `energy_first_hartree` | float | Energy of the first frame of the stitched path |
| `energy_ts_hartree` | float | TS energy |
| `energy_last_hartree` | float | Energy of the last frame of the stitched path |
| `endpoint_energy_orientation` | string | `"finished_first_to_finished_last"` |
| `forward_requested` / `backward_requested` | bool | Whether each direction was requested |
| `forward_integration_converged` / `backward_integration_converged` | bool\|null | Whether the direction stopped because the RMS-gradient stationarity criterion fired; diagnostic only, and always `false` under `--never-stop`, which bypasses that criterion. Combine it with `*_downhill_departure_valid` to check that the branch both left the TS downhill and met that criterion. |
| `forward_downhill_departure_valid` / `backward_downhill_departure_valid` | bool\|null | Whether the branch established a downhill departure from the TS |
| `forward_integration_stop_reason` / `backward_integration_stop_reason` | string\|null | Non-empty only for a numerical propagation failure |
| `never_stop` | bool | Whether `--never-stop` (bypass of the physical endpoint stops) was enabled |
| `never_stop_energy_bypasses` | int | Number of energy-rise or one-step energy-change stops actually bypassed |
| `rigid_projection` | object | Rigid-mode and initial/updated-Hessian provenance; see [projection provenance](#rigid-projection-provenance) |
| `rigid_projection.hessian_source` | string | Initial Hessian source: `"file"` (`--read-hess`), `"cache"` (earlier stage in the same run), or `"fresh"` (newly computed, including the EulerPC initialization when `irc.hessian_init` is not `calc`) |
| `bond_changes` | object | Directed first→last `{formed: [...], broken: [...]}`; omitted if the comparison was unavailable |
| `bond_changes_direction` | string | `"finished_first_to_finished_last"` when `bond_changes` is present |
| `files` | object | Trajectory and endpoint files (XYZ, plus PDB/CIF versions when available) |

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `scan_opt_mode` | string | Requested `grad` or `hess` optimizer preset |
| `scan_optimizer` | string | Effective optimizer: `lbfgs` or `rfo` |
| `n_stages` | int | Number of scan stages |
| `stages` | object[] | Per-stage data |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `files` | object | Output files |

**`stages[]`**: `n_steps`, `converged`, `pairs_1based`, `energies_hartree`, `final_energy_hartree`, `bond_changes`, and `optimizer_status` (`converged`/`not_converged`/`stalled`) plus `stop_reason` when the last optimizer of the stage stopped without converging.

### `scan2d` / `scan3d`

| Field | Type | Description |
|-------|------|-------------|
| `n_grid_points` | int | Total grid points |
| `n_points_attempted` | int | Fresh-run grid points attempted, excluding preoptimization |
| `n_points_usable` | int | Fresh-run points with explicit convergence, finite energy/coordinates, and a written geometry file |
| `grid_points` | object[] | Explicit grid-index, distances, energy, convergence, and `geometry_file` mapping |
| `current_output_paths` | string[] | CSV/HTML/PNG files and grid geometries written by this run; files left from an earlier run are not listed |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | Surface minimum energy |
| `charge` | int \| null | ML-region charge; null for plot-only `scan3d --csv` |
| `spin` | int \| null | ML-region multiplicity; null for plot-only `scan3d --csv` |
| `files` | object | CSV + plot files |

Fresh `scan2d`/`scan3d` results include the common calculator fields. Plot-only
`scan3d --csv` writes these keys as null. Plot-only results omit
`n_points_attempted` and include `n_points_usable` only when the imported CSV
records convergence and geometry files for every point.

### `path-opt`

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool \| null | Convergence flag: `true` / `false` from the engine's own convergence signal, `null` when it exposed none (`optimization_status` is then `"completed"`, never a success claim) |
| `mep_mode` | string | `"dmf"` or `"gsm"` |
| `image_energies_hartree` | float[] | All image energies |
| `n_images` | int | Image count |
| `hei_index` | int | Highest-energy image index |
| `barrier_kcal` | float | Forward barrier (kcal/mol) |
| `delta_kcal` | float | Reaction energy (kcal/mol) |
| `files` | object | Trajectory + HEI files |

### `path-search`

`path-search` has no `--out-json` flag. It writes `summary.json`; its fields
are listed in [`summary.json` (`path-search` / `all`)](#summary-json-path-search-all).

### `dft`

> **Note:** With `--out-json`, `dft` writes `result.json` and `summary.json` for
> both converged and non-converged SCF attempts, recording
> `scientific_status: "failed"` and `converged: false` on non-convergence.
> A non-converged SCF exits with code 1.
> An unhandled exception writes the standard `error` envelope.

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool | SCF converged? |
| `energy_hartree` / `energy_kcal_per_mol` | float | ML-region DFT energy; the same values as `model_dft_energy_*` |
| `model_dft_energy_hartree` / `model_dft_energy_kcal_per_mol` | float | ML-region DFT energy |
| `total_dft_mm_energy_hartree` / `total_dft_mm_energy_kcal_per_mol` | float | Recombined DFT/MM energy |
| `xc_functional` | string | XC functional |
| `basis_set` | string | Basis set |
| `used_gpu` | bool | GPU acceleration used? |
| `used_lowmem` / `lowmem_requested` | bool | Effective and requested low-memory state |
| `dft_settings` / `dft_resources` | object | Canonical scientific settings and effective host resources |
| `effective_ecp` | string/object \| null | Effective ECP passed to PySCF |
| `embedding` | object | Point-charge embedding: whether it is on, the cutoff, the number of charges, and a digest that identifies them |
| `charges` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `spin_densities` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `n_atoms` | int | ML-region atom count |
| `grid_level` | int | DFT grid level |
| `conv_tol` | float | SCF convergence tolerance |
| `max_cycle` | int | Maximum SCF iterations in effect after YAML and CLI are applied |
| `engine` | string | Engine label of the run (`pyscf(cpu)`, `gpu4pyscf`, or the low-memory GPU variant) |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `input_file` | string | Input structure path |
| `files` | object | `{"result_yaml": "result.yaml"}` |

### `extract`

| Field | Type | Description |
|-------|------|-------------|
| `n_atoms_extracted` | int | Atoms after extraction |
| `total_charge` | float | Computed total charge |
| `protein_charge` | float | Protein charge |
| `ligand_total_charge` | float | Ligand charge sum |
| `ion_total_charge` | float | Ion charge sum |
| `unknown_residue_charges` | object | `{resname: charge}` |
| `center` | string | Substrate specification (the `-c` value): PDB path, residue-ID list (e.g. `'A:123,B:456'`), or residue-name list (e.g. `'GPP,MMT'`) |
| `radius` | float | Extraction radius (angstrom) |
| `input_files` | string[] | Input PDB paths |
| `n_atoms_raw` | int | Atoms in the selected residues before backbone/truncation filtering (not the whole input structure) |
| `n_link_hydrogens` | int | Link H atoms added at severed bonds |
| `files` | object | Map of written files (pocket PDB per input, etc.) |
| `exclude_backbone` | bool | Value of `--exclude-backbone` for the run |
| `include_h2o` | bool | Value of `--include-h2o` for the run |
| `ligand_charge_input` | string | The `-l/--ligand-charge` argument as given |
| `ion_charges` | array | List of `[resname, charge]` pairs for ion residues encountered |

### `trj2fig`

| Field | Type | Description |
|-------|------|-------------|
| `n_frames` | int | Number of trajectory frames |
| `min_energy_hartree` / `max_energy_hartree` | float | Minimum and maximum frame energies |
| `energy_source` | string | `"trajectory_comment"` or `"mlip_recomputed"` |
| `mlip_backend` / `mlip_model` / `mlip_precision` | string \| null | Calculator used to recompute the energies; all are null in trajectory-comment mode |
| `charge` / `multiplicity` | int \| null | Charge and multiplicity used for the recomputation; null in comment mode. When only one is given, the other defaults to 0 or 1. |
| `output_files` | string[] | Ordered paths of every output; keeps files with the same basename in different directories |
| `files` | object | Basename-to-path map; when two outputs share a basename only one is kept, so prefer `output_files` |

Supplying either `-q/--charge` or `-m/--multiplicity` recomputes every frame with the selected MLIP. This rescoring uses the MLIP alone: the command takes no topology or ML-region input and does not calculate an ONIOM energy.

### `energy-diagram`

| Field | Type | Description |
|-------|------|-------------|
| `n_points` | int | Number of energy data points |
| `files` | object | Output diagram filename-to-path map |

### `bond-summary`

With `--json`, `bond-summary` prints JSON to **stdout** and writes no `result.json`; redirect stdout to keep it.

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` / `scientific_status` | string / string | `completed` / `success` when every pair was compared. If any pair could not be compared, `execution_status` is `failed` and `scientific_status` is `partial` or `failed`, and the command exits with code 1. |
| `comparisons` | object[] | Per-pair comparison with `structure_a`, `structure_b`, `bonds_formed` / `bonds_broken` (counts), and `formed` / `broken` (one record per bond: `atom_i`, `atom_j`, `element_i`, `element_j`, `distance_a_angstrom`, `distance_b_angstrom`); a pair that could not be compared has `error` instead. |

### Rigid projection provenance

`freq`, `irc`, and `tsopt` results include a `rigid_projection` object; `opt` includes it when `--flatten` runs. `freq --dump` also writes the same object to `thermoanalysis.yaml`.

| Field | Type | Description |
|-------|------|-------------|
| `treatment` | string | Fixed rigid-mode treatment: `"constrained"` |
| `algorithm` | string | Name of the projection method |
| `effective_rank` | int | Number of rigid directions removed from the Hessian of the movable atoms |
| `full_rigid_rank` | int | Rank of the rigid motions of the whole system before the frozen atoms are taken into account |
| `frozen_constraint_rank` | int | Rank removed because the frozen atoms must stay in place |
| `svd_rtol` | float | Relative SVD tolerance used for the rank decision |
| `active_atom_count` / `frozen_atom_count` | int | Movable and frozen atom counts |
| `active_atoms` / `frozen_atoms` | int[] | 0-based indices of the movable and frozen atoms |
| `hessian_space` | string | `"full"` or `"active"` input Hessian space |
| `hessian_source` / `source` | string | Hessian provenance. `freq`/`irc` use `hessian_source`: `"file"` (`--read-hess`), `"cache"` (earlier stage in the same run), or `"fresh"`; `opt`/`tsopt` use `source`. |
| `hessian_shape` / `raw_hessian_shape` | int[2] | Input Hessian shape. `freq`/`irc` use `hessian_shape`; `opt`/`tsopt` use `raw_hessian_shape` (`freq` records both). |
| `near_zero_mode_count` / `near_zero_frequencies_cm` | int / float[] | Number and values of the modes within ±`frequency_zero_cutoff_cm` (5.00 cm⁻¹ by default); these modes are also in the full frequency list |

`constrained` removes only the rigid motions of the whole system that keep the frozen atoms in place; see [freq](freq.md#rigid-modes-with-frozen-boundaries).

(summary-json-path-search-all)=
## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json` with a richer structure:

| Field | Type | Description |
|-------|------|-------------|
| `execution_status` / `scientific_status` | string / string | Execution completeness and completion of requested numerical/calculation stages. |
| `scientific_status_reasons` | string[] | Reasons for missing or unusable requested results; omitted on success. |
| `pipeline_stop` | object \| absent | Present only on an early stop. `stage` is `post` (`reason` `no_segments` / `no_reactive_segment`), `before_irc` (a TSOPT reason, plus `segment` and `tsopt_result`), or `endpoint_opt` (`endpoint_execution_failed` and endpoint-specific `failures`). Rendered in `summary.log` as `Pipeline stop`. |
| `expected_item_ids` / `observed_item_ids` | string[] | Expected and observed stage identifiers. |
| `config` | object | Effective settings. `mep_mode` identifies GSM/DMF; `ts_opt_mode` and `endpoint_opt_mode` identify the configured post-processing presets. Generic `opt_mode*` keys record the CLI values in effect. `path_opt_mode` is the single-structure optimizer used for endpoint preoptimization (see `preopt`), not the MEP path algorithm. |
| `n_segments` | int | Segment count |
| `search_max_depth` | int | Effective recursion cap; `0` means subdivision was disabled |
| `path_optimizers` | string[] | Single-structure optimizers actually used during path preparation/refinement (`lbfgs`, `rfo`); includes scan and alignment work in `all`. Also present in `path-opt` `result.json` |
| `preopt_requested` / `preopt_converged` | bool / bool \| null | Whether endpoint preoptimization ran, and whether every endpoint converged; `null` when any endpoint reported no readable signal. In `all`, `preopt_converged` counts toward `scientific_status` unless the requested final TS and both endpoint optimizations have converged for every reactive segment; the field itself is still reported |
| `segments` | object[] | Per-segment `index` (1-based), `tag` (a label such as `seg_001`; a kink segment, which has no covalent bond change, has `kink` in its tag), `kind` (`"seg"` for a reactive segment, `"bridge"` for a [bridge segment](path-search.md#how-it-works), `"tsopt"` in TS-only mode), `converged` (whether every optimization that built the segment converged), `barrier_kcal`, `delta_kcal`, and `bond_changes` (bridge segments emit `""`). `barrier_kcal` is the barrier on the MEP before TS optimization; in TS-only mode it is TS − R, where R is the higher-energy IRC endpoint (a name for the direction, not the chemical direction of the reaction). |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Backend name (`uma`, `orb`, `mace`, `aimnet2`, `dft`, or `custom`) |
| `mlip_model` | string \| null | Exact model/checkpoint name, or `FUNCTIONAL/BASIS` for `dft` |
| `mlip_model_label` | string \| null | Publication-facing model label |
| `mlip_task` | string \| null | Backend task for a multi-domain model |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64`; null for DFT and custom calculators (the DFT engine is recorded separately) |
| `charge` | int | ML-region charge |
| `spin` | int | ML-region multiplicity |
| `environment` | object | Hardware info |
| `references` | object[] | Methods actually used by the run, as `{method, citation, doi}` records. The same reference set is grouped at the end of `summary.log` and final stdout immediately before elapsed time. |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `n_segments_reactive` | int | Number of non-bridge (reactive) segments |
| `rate_limiting_step` | object | Highest local barrier among the reactive segments, as `{segment, barrier_kcal, method}` plus `mep_barrier_kcal` (the barrier on the MEP; not in TS-only mode), taken at the highest-level method available for every segment (`DFT//MLIP/MM_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`; with `--backend dft`, `DFT/MM_Gibbs` and `DFT/MM` take the place of `MLIP_Gibbs` and `MLIP`). It is not a microkinetic rate-limiting-step assignment. |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `overall_reaction_energy_method` | string | Method of the overall reaction energy, from the same list as `rate_limiting_step.method` and no higher than it |
| `post_segments` | list | Per-segment TS/IRC/freq/DFT results |
| `post_segments[].tsopt.n_imaginary_modes` / `.imaginary_frequencies_cm` | int / float[] | n_imag of the optimized TS and its imaginary frequencies (cm⁻¹, negative) |
| `post_segments[].mlip` / `.gibbs_mlip` / `.dft` / `.gibbs_dft_mlip` | object | R, TS, and P energies at one level, with `barrier_kcal`, `delta_kcal`, and `energies_kcal`: ML/MM electronic energies (`--tsopt`), ML/MM Gibbs energies (`--thermo`), DFT energies (`--dft`), and DFT//ML/MM Gibbs energies (`--thermo` with `--dft`) |
| `post_segments[].tsopt.energy_valid` / `.structure_valid` | bool | `energy_valid`: the final TS energy is a finite number. `structure_valid`: the final TS structure file exists and its coordinates are finite. These checks run no extra Hessian or optimization. |
| `post_segments[].tsopt.n_opt_cycles` / `.max_cycles` | int / int\|null | TS optimization cycles executed and configured limit. These are reported for both converged and normally non-converged runs. |
| `post_segments[].irc` / `.endpoint_assignment` / `.endpoint_opt` | object | IRC stop diagnostics, endpoint orientation, and endpoint-OPT convergence, respectively. In TS-only mode, `endpoint_assignment.policy` is `higher_energy_endpoint_as_reactant` and `chemical_direction_known` is false. `endpoint_opt.reactant` and `.product` report `optimization_status`, `n_opt_cycles`, `max_cycles`, and any `stop_reason`. How the IRC stopped and whether the bond changes match do not enter `scientific_status`; whether the optimized endpoints are the intended R and P is for you to check. |
| `post_segments[].thermo_symmetry` | object | Child-reported point-group and rotational-symmetry provenance by state: R/TS/P, including TS-only runs. States with valid symmetry-number provenance are included; missing states are omitted, and the field is absent only when no state has valid provenance. |
| `key_output_files` | object | Current-run output index: root filename → description; each `seg_NN` entry is `{description, files}` with paths relative to that segment directory. |
| `current_output_paths` | string[] | Sorted paths relative to `--out-dir`, limited to files written by the current run. |

## Usage examples

### Python

The script reads `result_opt/result.json`, which `opt --out-json` writes.

```python
import json

with open("result_opt/result.json") as f:
    result = json.load(f)

status = result.get("optimization_status")
if result["execution_status"] == "failed":
    raise RuntimeError(f"{result['error_type']}: {result['error']}")
elif status == "converged":
    print(f"Energy: {result['energy_hartree']:.6f} Hartree")
elif status in {"not_converged", "stalled"}:
    print(f"Not converged after {result['n_opt_cycles']} cycles")
    print(f"Max force: {result['final_max_force']:.6f}")
else:
    print(f"Status: {status}")
```

### jq

```bash
# Check convergence
jq '{execution_status, scientific_status}' result.json

# Get barrier from path-opt
jq '.barrier_kcal' result.json

# List imaginary frequencies from tsopt
jq '.imaginary_frequencies_cm' result.json

# Get thermochemistry from freq
jq '.thermochemistry.sum_EE_and_thermal_free_energy_ha' result.json

# Get each segment's barrier from all (after --tsopt)
jq '.post_segments[] | {index, barrier_kcal: .mlip.barrier_kcal}' result_all/summary.json
```

## Notes

- A run that fails while the CLI options or the input are being checked (before the output directory is set up) can stop without writing any JSON. Treat a nonzero exit code as a failure, and read stderr or the job log for the message.
- `all` and `path-search` write `summary.json` only once the run reaches its summary step, so an early input error leaves no file.

## See Also

- {ref}`Exit codes <exit-codes>` — what each exit code means
- [Output layout](output-layout.md) — where `result.json` and `summary.json` are written
- [Troubleshooting](troubleshooting.md) — what to change after a failed or non-converged run
- [YAML Reference](yaml-reference.md) — configuration inputs whose values surface in these schemas
- [all](all.md), [path-search](path-search.md) — subcommands that write `summary.json` without `--out-json`
- [opt](opt.md), [sp](sp.md), [tsopt](tsopt.md), [freq](freq.md), [irc](irc.md), [scan](scan.md), [scan2d](scan2d.md), [scan3d](scan3d.md), [path-opt](path-opt.md), [dft](dft.md), [extract](extract.md), [trj2fig](trj2fig.md), [energy-diagram](energy-diagram.md), [bond-summary](bond-summary.md) — subcommands that write JSON only with `--out-json` (`--json` for `bond-summary`)
