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
| `status` | string | Supplied by each subcommand (not auto-injected): `all`/`path-search` use `success`/`partial`/`failed`, the error path uses `error`, and per-stage subcommands use command-specific values (e.g. `converged`/`not_converged`, `completed`) |
| `elapsed_seconds` | float | Wall-clock time (seconds) |
| `environment` | object | Hardware info (see below) |

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

## Subcommand schemas

### `opt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"converged"` or `"not_converged"` |
| `energy_hartree` | float | Final ONIOM energy (Hartree) |
| `n_opt_cycles` | int | Optimization cycles completed |
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"lbfgs"`, `"rfo"` (aliases: `light`/`lbfgs` → `grad`; `heavy`/`rfo` → `hess`) |
| `backend` | string | ML backend (`"uma"`, `"orb"`, `"mace"`, `"aimnet2"`) |
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
| `files` | object | Output file map |

### `tsopt`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `energy_hartree` | float | TS energy (Hartree) |
| `n_imaginary_modes` | int | Number of imaginary frequencies |
| `imaginary_frequencies_cm` | float[] | Imaginary frequencies (cm$^{-1}$, negative) |
| `opt_mode` | string | One of `"grad"`, `"hess"`, `"light"`, `"heavy"`, `"dimer"`, `"rsirfo"`, `"trim"`, `"rsprfo"` (aliases: `light`/`dimer` → `grad` (PHG-Dimer); `heavy`/`rsirfo` → `hess` (RS-I-RFO); `trim` → TRIM; `rsprfo` → RS-P-RFO) |
| `n_atoms` | int | Total atoms |
| `n_opt_cycles` | int | Optimization cycles |
| `backend` | string | ML backend |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
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
| `backend` | string | ML backend |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `n_atoms` | int | Total atoms |
| `n_freeze_atoms` | int | Frozen atoms |
| `files` | object | `{"frequencies_txt": "frequencies_cm-1.txt"}` |

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
| `n_frames_forward` / `backward` / `total` | int | IRC frames |
| `energy_reactant_hartree` | float | Reactant energy |
| `energy_ts_hartree` | float | TS energy |
| `energy_product_hartree` | float | Product energy |
| `forward_converged` / `backward_converged` | bool | IRC convergence |
| `backend` | string | ML backend |
| `bond_changes` | object | `{formed: [...], broken: [...]}` |
| `files` | object | Trajectory files (xyz + pdb) |

### `scan`

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"completed"` |
| `n_stages` | int | Number of scan stages |
| `stages` | object[] | Per-stage data |
| `backend` | string | ML backend |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `files` | object | Output files |

**`stages[]`**: `n_steps`, `converged`, `pairs_1based`, `energies_hartree`, `final_energy_hartree`, `bond_changes`

### `scan2d` / `scan3d`

| Field | Type | Description |
|-------|------|-------------|
| `n_grid_points` | int | Total grid points |
| `pair1`, `pair2` (,`pair3`) | object | `{i, j, low, high}` |
| `min_energy_hartree` | float | Surface minimum energy |
| `backend` | string | ML backend |
| `charge` | int | Model-region charge |
| `spin` | int | Model-region multiplicity |
| `files` | object | CSV + plot files |

### `path-opt`

| Field | Type | Description |
|-------|------|-------------|
| `converged` | bool | Convergence flag |
| `mep_mode` | string | `"dmf"` or `"gsm"` |
| `backend` | string | ML backend |
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
| `energy_hartree` | float | DFT energy |
| `xc_functional` | string | XC functional |
| `basis_set` | string | Basis set |
| `used_gpu` | bool | GPU acceleration used? |
| `backend` | string | ML backend for ONIOM high-level region |
| `charges` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `spin_densities` | object | `{mulliken, lowdin, iao}` per-atom arrays |
| `n_atoms` | int | QM-region atom count |
| `grid_level` | int | DFT grid level |
| `conv_tol` | float | SCF convergence tolerance |
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
| `ion_charges` | object | Map `{resname: charge}` for ion residues encountered |

## `summary.json` (`path-search` / `all`)

The `all` and `path-search` commands write `summary.json`:

| Field | Type | Description |
|-------|------|-------------|
| `status` | string | `"success"` / `"partial"` |
| `n_segments` | int | Segment count |
| `segments` | object[] | Per-segment barrier, delta, bond changes |
| `energy_diagrams` | object[] | Energy profiles with labels and kcal/mol values |
| `mlip_backend` | string | Backend name (`uma`, `orb`, `mace`, or `aimnet2`) |
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
