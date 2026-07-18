# JSON Output Reference

mlmm provides machine-readable JSON output for programmatic consumption by AI agents, scripts, and downstream tools.

## `--out-json` flag

Most MLIP-based subcommands (`opt`, `sp`, `tsopt`, `freq`, `irc`, `scan`, `scan2d`, `scan3d`, `path-opt`, `dft`, `extract`) support `--out-json / --no-out-json` (default: off).
When enabled, a `result.json` file is written to the output directory alongside the normal outputs.

```bash
mlmm opt -i r_complex_layered.pdb --max-cycles 5 --out-json --out-dir result_opt
cat result_opt/result.json | python -m json.tool
```

The `all` and `path-search` commands always write `summary.json` (no `--out-json` flag needed).

### `summary.json` mirror

`write_result_json` mirrors every per-stage `result.json` payload to `summary.json` alongside it. MCP clients and agent scripts can read a single filename (`summary.json`) across every subcommand; the `result.json` written next to it carries the identical payload.

## Common envelope

Every `result.json` (and the mirrored `summary.json`) automatically includes:

| Field | Type | Description |
|-------|------|-------------|
| `schema_version` | string | Envelope schema version; current value comes from `mlmm.core.utils.RESULT_JSON_SCHEMA_VERSION` — pin against that constant rather than the literal in this doc. Bumps signal a structural change. |
| `command` | string | Subcommand name (e.g. `"opt"`) |
| `mlmm_version` | string | Package version |
| `status` | string | Command-specific: `all`/`path-search` use `success`/`partial`/`failed`; `opt` uses `converged`/`not_converged`/`stalled`; `tsopt` uses `converged`/`not_converged`/`stalled`/`unverified`; completed analysis/integration stages use `completed`; exception envelopes use `error`. |
| `elapsed_seconds` | float | Wall-clock time (seconds) |
| `environment` | object | Hardware info (see below) |
| `run_id` | string | Present when an orchestrator (including MCP) assigns a current invocation identity; conflicting caller values are rejected. |

MLIP/ML/MM calculator stages additionally record:

| Field | Type | Description |
|-------|------|-------------|
| `mlip_backend` | string \| null | Backend identifier (`uma`, `orb`, `mace`, `aimnet2`, or `custom`); null when a plot-only command did not evaluate a calculator |
| `mlip_model` | string \| null | Exact model/checkpoint; `filename:factory` for `--calc-file` |
| `mlip_precision` | string \| null | Effective public precision (`fp32` or `fp64`); null for custom calculators |
| `mm_backend` | string \| null | MM Hessian/energy backend (`hessian_ff` or `openmm`); null when a plot-only command did not evaluate a calculator |
| `link_atom_method` | string \| null | Link-atom placement (`scaled` or `fixed`); null for plot-only output |
| `use_cmap` | bool \| null | Whether CMAP terms were enabled; null for plot-only output |

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

An optimizer may also report `"status": "stalled"`: the energy stopped decreasing over the configured window (an energy plateau) while the configured force/step convergence criteria remained unmet. A stall is a distinct, non-converged outcome — it is never reported as `converged` — and it stops further flatten/retry work rather than repeating a non-progressing optimization. When present, a `stop_reason` string records the energy range, window, and the failed criteria. A stall may be retried (e.g. from a perturbed geometry or with tighter step control); it is not an alias for `max_cycles` exhaustion or a generic failure. In microiteration, a stalled macro step or a stalled latest micro (MM) relaxation is reported truthfully and never masquerades as macro convergence.

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
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"lbfgs"`, `"rfo"` (aliases: `light`/`lbfgs` → `grad`; `heavy`/`rfo` → `hess`) |
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
| `status` | string | `"converged"` only when the optimizer converged and `n_imaginary_modes == 1`; otherwise `"not_converged"`, or `"unverified"` with `--skip-final-freq`. An energy-plateau `"stalled"` outcome (see above) wins over all of these and is never reported as `converged`; the dimer (grad) mode also reports `stalled`. |
| `energy_hartree` | float | TS energy (Hartree) |
| `n_imaginary_modes` | int | Number of imaginary frequencies |
| `imaginary_frequencies_cm` | float[] | Imaginary frequencies (cm$^{-1}$, negative) |
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"dimer"`, `"rsirfo"`, `"trim"`, `"rsprfo"` (aliases: `light`/`dimer` → `grad` (PHG-Dimer); `heavy`/`rsirfo` → `hess` (RS-I-RFO); `trim` → TRIM; `rsprfo` → RS-P-RFO) |
| `n_atoms` | int | Total atoms |
| `n_opt_cycles` | int | Optimization cycles |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `reference_mode_file` | string\|null | Advanced path-derived mode supplied with `--ref-mode` |
| `safeguards` | object | Heavy-mode rejection/recovery, exact-saddle, and target-mode diagnostics |
| `rigid_projection` | object | Frozen-boundary TR provenance for Dimer/flatten/final saddle analysis |
| `files` | object | Final geometry + vib mode files |

### `freq`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_modes` | int | Total normal modes |
| `n_imaginary` | int | Imaginary frequency count |
| `frequencies_cm` | float[] | All frequencies (cm$^{-1}$) |
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
| `forward_converged` / `backward_converged` | bool\|null | Directional convergence flags |
| `never_stop` | bool | Whether opt-in energy-rise/plateau bypass mode was enabled |
| `never_stop_energy_bypasses` | int | Number of energy-rise/plateau stops actually bypassed |
| `rigid_projection` | object | Frozen-boundary TR provenance for the initial/updated Hessian |
| `bond_changes` | object | Directed first→last `{formed: [...], broken: [...]}`; omitted if comparison was unavailable |
| `bond_changes_direction` | string | `"finished_first_to_finished_last"` when bond changes are present |
| `files` | object | Trajectory and endpoint files (XYZ plus available PDB/CIF companions) |

**`rigid_projection` provenance:** the object records the selected treatment
(`treatment`), `effective_rank`, active/frozen atom counts and indices, and the
Hessian source/shape used by that workflow. `constrained` is the default;
`legacy-active` is an isolated-active comparison treatment. A `freq --dump`
run writes the same object to `thermoanalysis.yaml`. Field names for the final
two values follow the producing workflow (`hessian_source` / `hessian_shape`,
or `source` / `raw_hessian_shape`).

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
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
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | Surface minimum energy |
| `charge` | int \| null | Model-region charge; null for plot-only `scan3d --csv` |
| `spin` | int \| null | Model-region multiplicity; null for plot-only `scan3d --csv` |
| `files` | object | CSV + plot files |

Fresh `scan2d`/`scan3d` results include the common MLIP/ML/MM calculator
provenance. Plot-only `scan3d --csv` keeps the same keys but writes null because
the imported energy grid does not identify the calculator that produced it.

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
| `include_h2o` | bool | Value of `--include-H2O` at run time |
| `ligand_charge_input` | string | Raw `-l/--ligand-charge` argument |
| `ion_charges` | array | List of `[resname, charge]` pairs for ion residues encountered |

## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json`:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | Segment count |
| `segments` | object[] | Per-segment barrier, delta, bond changes |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Backend name (`uma`, `orb`, `mace`, `aimnet2`, or `custom`) |
| `mlip_model` | string \| null | Exact model/checkpoint name, recorded separately from the backend |
| `mlip_precision` | string \| null | Effective `fp32` / `fp64`; null for custom calculators |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `environment` | object | Hardware info |

The `all` command additionally includes:

| Field | Type | Description |
|-------|------|-------------|
| `n_segments_reactive` | int | Number of non-bridge (reactive) segments |
| `rate_limiting_step` | object | RLS segment index and barrier |
| `overall_reaction_energy_kcal` | float | Overall reaction energy |
| `post_segments` | list | Per-segment TS/IRC/freq/DFT results |

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
