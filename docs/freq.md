# `freq`

`frequencies_cm-1.txt`, JSON `frequencies_cm`, and `n_modes` retain the complete signed physical spectrum after the existing frozen-atom and rigid-mode projection. `--max-write` and `--sort` control only which mode files are written and their order. `n_imaginary` (YAML `num_imag_freq`) is the resolved negative count below the reporting threshold; `n_negative_modes` also includes weak negative modes. With `frequency_representation: complete`, `near_zero_frequencies_cm` is a subset of the complete array; do not append it and count modes twice.

Thermochemistry retains the existing QRRHO policy (100 cm⁻¹ rotor cutoff, no imaginary inversion and no positive-frequency floor), including positive low-frequency modes. Changing `freq.zero_cutoff_cm` does not change thermal values computed from the same complete spectrum.


Compute ML/MM vibrational frequencies and thermochemistry (zero-point energy (ZPE), Gibbs energy, etc.) on a layered enzyme PDB, with partial-Hessian vibrational analysis (PHVA) support.

**When to use `mlmm freq`:**

- Validate stationary-point character of an optimized minimum, transition state, or IRC endpoint (a minimum has no imaginary frequencies; a transition state has exactly one).
- Compute quasi-rigid-rotor-harmonic-oscillator (QRRHO) thermochemistry.

The command runs vibrational analysis with the ML/MM calculator, honoring frozen atoms via PHVA. It exports normal-mode trajectories as `_trj.xyz` and `.pdb` (mapped back onto the enzyme ordering), and prints a Gaussian-style thermochemistry summary when the optional `thermoanalysis` package is installed.

Imaginary frequencies appear as negative values. Runtime and memory depend on
the backend and system; compare `Analytical` and `FiniteDifference` on a
representative pilot.

## Examples

Here, `pocket.pdb` contains the full system matching `real.parm7`; `ml_region.pdb` selects the ML atoms.

Basic frequency analysis:

```bash
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --out-dir ./result_freq
```

Limit the number of exported modes for quick inspection:

```bash
# Limit the number of exported modes for quick inspection
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --max-write 6 --out-dir ./result_freq_quick
```

PHVA with explicit frozen atoms and dump thermo payload:

```bash
# PHVA with explicit frozen atoms and dump thermo payload
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

Analytical Hessian mode on VRAM-rich nodes:

```bash
# Analytical Hessian mode on VRAM-rich nodes
mlmm freq -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

## Workflow

1. **ML/MM calculator setup** — The ML region is supplied via `--model-pdb`; Amber parameters are read from `--parm`. `--hessian-calc-mode` selects analytical or finite-difference Hessians. The calculator may return either the full 3N x 3N Hessian or an active degree-of-freedom (DOF) sub-block.
2. **PHVA & translation/rotation (TR) projection** — With frozen atoms, eigenanalysis occurs inside the active subspace. The default constrained projector removes only full-system rigid motions that leave every frozen anchor fixed; it does not treat the active fragment as an isolated molecule. Both 3N x 3N and active-block Hessians are accepted, and frequencies are reported in cm^-1 (negatives = imaginary).
3. **Active DOF mode** — `--active-dof-mode` selects which atoms enter the analysis (default `partial`); see the CLI options table for the four modes.
4. **Mode export** — `--max-write` limits how many mode trajectories are written. Modes are sorted by value (or absolute value with `--sort abs`). Each exported mode writes `_trj.xyz` and `.pdb` trajectories mapped back onto the enzyme ordering. The sinusoidal trajectory amplitude (`--amplitude-ang`) and frame count (`--n-frames`) match the YAML defaults.
5. **Thermochemistry** — If `thermoanalysis` is installed, a QRRHO-like summary (E, ZPE, E/H/G corrections, heat capacities, entropies) is printed using PHVA frequencies. The structure energy is labeled in Hartree as `E + G_corr = G` (electronic energy + Gibbs free-energy correction = Gibbs free energy). CLI pressure in atm is converted internally to Pa. The molecular point group and external rotational symmetry number are detected independently for each analyzed structure, and the resulting `1/sigma` correction is always included. An expert can override the detected number with `thermo.symmetry_number` in YAML. When `--dump`, a `thermoanalysis.yaml` snapshot is also written. **Frequency-treatment policy**: `freq` applies the **standalone-freq policy** — QRRHO with a 100 cm⁻¹ rotor cutoff, unit frequency/ZPE scaling, **no** imaginary-frequency inversion, and **no** positive-frequency floor. This is deliberately different from the internal `Geometry.get_thermoanalysis` policy used by some bundled-engine paths, which additionally inverts small imaginaries (from −15 cm⁻¹) and floors positive frequencies below 25 cm⁻¹. Neither is a universal scientific default; each is tied to its entry point. The effective policy (`kind`, `rotor_cutoff_cm`, `frequency_scale`, `zpe_scale`, `invert_imag_from_cm`, `positive_frequency_floor_cm`) is serialized under `thermo_policy` in `thermoanalysis.yaml` and in `result.json`.
6. **Device selection** — `ml_device="auto"` triggers CUDA when available, otherwise CPU. The internal TR projection/mode assembly runs on the same device to minimize transfers.
7. **Exit behavior** — Keyboard interrupts exit with code 130; other failures print a traceback and exit with code 1.

### Frozen-boundary TR projection

The fixed constrained treatment is used for PHVA. It starts from
the full system's rigid translations and rotations, then retains only
components that do not move any frozen anchor. The generic effective ranks are:

| Frozen-anchor geometry | Effective rank removed |
| --- | ---: |
| none | 6 |
| one anchor | 3 |
| two distinct anchors | 1 |
| at least three non-collinear anchors | 0 |

Realistic ML/MM boundaries normally have several non-collinear anchors, so the
effective rank is usually zero and no active-space direction is removed. An
all-frozen selection has no active DOF and raises an explicit error.

A stale non-constrained `geom.tr_projection` value fails explicitly.

With `--out-json`, `result.json.rigid_projection` records the treatment,
effective rank, Hessian source, and Hessian shape. `--dump` records the same
provenance in `thermoanalysis.yaml`.

## Outputs

```text
out_dir/ (default: ./result_freq/)
├─ result.json                      # Present with --out-json; includes rigid_projection provenance
├─ mode_XXXX_±freqcm-1_trj.xyz   # Per-mode trajectory
├─ mode_XXXX_±freqcm-1.pdb       # PDB trajectory mapped back onto the enzyme ordering
├─ frequencies_cm-1.txt           # Full frequency list using the selected sort order
└─ thermoanalysis.yaml            # Present when thermoanalysis is importable and --dump is True
```
- Console blocks summarizing resolved `geom`, `calc`, `freq`, and thermochemistry settings.

## CLI options

`mlmm freq --help` shows core options; `mlmm freq --help-advanced` shows the full option list. The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
| `-i, --input PATH` | Full enzyme PDB (no link atoms). | Required |
| `--parm PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Explicit ML-region atom indices (alternative to `--model-pdb`). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Indexing convention for `--model-indices`. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Automatically detect ML/MM layers from B-factors. | Enabled |
| `-q, --charge INT` | ML region charge. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g., `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--ref-pdb FILE` | Reference PDB topology for non-PDB inputs. | _None_ |
| **Backend & compute** | | |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--precision [fp32\|fp64]` | MLIP backend precision; unset uses UMA/AIMNet2 fp32 and ORB/MACE fp64. AIMNet2 rejects fp64. | backend-specific |
| `--workers INT` | UMA predictor workers. Values greater than 1 require `fairchem-core[extras]` and cannot be combined with `Analytical`. | `1` |
| `--workers-per-node INT` | Workers per node for the parallel UMA predictor. | _None_ |
| `--mm-backend [hessian_ff\|openmm]` | MM backend. Hessians use finite differences by default; set `calc.mm_fd: false` for the `hessian_ff` analytical path. | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled ($g$-factor) or fixed 1.09/1.01 Å. | `scaled` |
| `--cmap/--no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
| `--hess-device CHOICE` | Device for post-evaluation Hessian placement and diagonalization: `auto`, `cuda`, `cpu`. It does not move Hessian evaluation/assembly; `cpu` moves the evaluated matrix before diagonalization. | `auto` |
| **Active-region freezing & Hessian** | | |
| `--freeze-atoms TEXT` | 1-based comma-separated frozen atom indices. | _None_ |
| `--active-dof-mode CHOICE` | Active DOF selection: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--hess-cutoff FLOAT` | Cutoff distance for Hessian-target MM atoms. | _None_ |
| `--movable-cutoff FLOAT` | Cutoff distance for movable-MM layer. | _None_ |
| `--hessian-calc-mode CHOICE` | Hessian mode (`Analytical` or `FiniteDifference`). | `FiniteDifference` |
| `--dump-hess PATH` | Save Hessian, atom order, Cartesian geometry, active-DOF basis, PHVA metadata, model charge, and multiplicity to `.npz` for a matching `mlmm irc --read-hess` run. | _None_ |
| **Mode export** | | |
| `--max-write INT` | Number of modes to export. | `10` |
| `--sort CHOICE` | Mode ordering: `value` (cm^-1) or `abs`. | `value` |
| `--amplitude-ang FLOAT` | Mode-trajectory amplitude (angstrom). | `0.8` |
| `--n-frames INT` | Frames per mode trajectory. | `20` |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB companions when a PDB template is available. | `True` |
| **Thermochemistry** | | |
| `--temperature FLOAT` | Thermochemistry temperature (K). | `298.15` |
| `--pressure FLOAT` | Thermochemistry pressure (atm). | `1.0` |
| `--dump/--no-dump` | Write `thermoanalysis.yaml`. | `False` |
| **Output & config** | | |
| `-o, --out-dir TEXT` | Output directory. | `./result_freq/` |
| `--out-json/--no-out-json` | Write machine-readable `result.json` to `out_dir`. | `False` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running frequency analysis. Shown in `--help-advanced`. | `False` |

The handoff is identity-checked: IRC rejects files from a different atom order,
geometry, layer selection, Hessian active basis, model charge, or multiplicity.
Schema-1 files predate electronic-state identity and are rejected unless
`--allow-unverified-hess-state` is explicitly supplied to IRC after independent
state verification.

## YAML configuration

An explicit analytical Hessian with `workers > 1` is rejected. Use one worker
for analytical curvature, or select `FiniteDifference` before enabling the UMA
parallel predictor.

Provide mappings with merge order **defaults < config < explicit CLI**.
Shared sections reuse [YAML Reference](yaml-reference.md).
An additional `thermo` section is supported for thermochemistry controls.

```yaml
geom:
 coord_type: cart                  # coordinate type: cartesian vs dlc internals
 freeze_atoms: []                  # 1-based frozen atoms merged with CLI/link detection
 tr_projection: constrained        # fixed internal PHVA treatment
calc:
 model_charge: 0                   # net charge (CLI override)
 model_mult: 1                     # spin multiplicity 2S+1
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
 backend: uma                      # MLIP backend: uma | orb | mace | aimnet2
 uma_model: uma-s-1p2              # uma-s-1p2 | uma-m-1p1
 uma_task_name: omol                # UMA task name (UMA backend only)
 ml_device: auto                   # ML backend device selection
 hessian_calc_mode: FiniteDifference   # Compare both modes on a representative pilot
 out_hess_torch: true              # request torch-form Hessian
 mm_fd: true                       # MM finite-difference toggle
 return_partial_hessian: true      # allow partial Hessians (PHVA default)
freq:
 zero_cutoff_cm: 5.0               # Near-zero classification window (cm^-1); modes are retained
 amplitude_ang: 0.8                # displacement amplitude for modes (Å)
 n_frames: 20                      # number of frames per mode
 max_write: 10                     # maximum number of modes to write
 sort: value                       # sort order: value vs abs
thermo:
 temperature: 298.15               # thermochemistry temperature (K)
 pressure_atm: 1.0                 # thermochemistry pressure (atm)
 symmetry_number: null             # auto-detect; positive integer overrides
 dump: false                       # write thermoanalysis.yaml when true
```

## See Also

- [tsopt](tsopt.md) — Optimize TS candidates (validate with freq/IRC; expected: one imaginary frequency)
- [opt](opt.md) — Geometry optimization (often precedes freq)
- [dft](dft.md) — Single-point DFT for higher-level energy evaluation
- [all](all.md) — End-to-end workflow with `--thermo`
- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [YAML Reference](yaml-reference.md) — Full `freq` and `thermo` configuration options
- [Glossary](glossary.md) — Definitions of ZPE, Gibbs Energy, Enthalpy, Entropy
