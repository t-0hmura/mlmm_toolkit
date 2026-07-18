# `irc`

Runs EulerPC-based IRC (Intrinsic Reaction Coordinate) integration from a transition state toward reactants and products using the ML/MM calculator. Use it to validate that an optimized TS connects the expected reactant and product, or to generate reactant/product structures for downstream thermochemistry and DFT single-point energy evaluation — typically as `tsopt` -> `freq` (confirm **one** imaginary mode) -> `irc`. By default both forward and backward branches are computed. `mlmm irc` keeps the CLI intentionally narrow; parameters not surfaced on the command line should be provided via YAML so the run remains explicit and reproducible. The common input bridge accepts PDB/mmCIF and `geom_loader` formats. With a PDB/mmCIF topology (direct input or `--ref-pdb`) and conversion enabled, trajectories receive PDB companions; mmCIF and oversized-PDB bridge inputs also receive CIF companions with restored identifiers.

## Examples

```bash
# Minimal run from a TS PDB
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --max-cycles 50 --out-dir ./result_irc
```

Forward branch only:

```bash
# Forward branch only
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --no-backward --out-dir ./result_irc_forward
```

Smaller step size with analytical Hessians:

```bash
# Smaller step size for a shallow surface
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 --no-detect-layer -q 0 -m 1 --step-size 0.05 \
 --hessian-calc-mode Analytical --out-dir ./result_irc_analytical
# keep both branches and raise the step limit with --max-cycles 150
```

If an IRC stops almost immediately, first reduce `--step-size` (for example,
from 0.10 to 0.05 Bohr). If a verified small shoulder still triggers only the
energy-rise/plateau stop, opt in to `--never-stop`:

```bash
mlmm irc -i ts.pdb --parm real.parm7 --model-pdb ml_region.pdb -q 0 \
 --step-size 0.05 --never-stop --max-cycles 250 -o result_irc_continue
```

This is not an unlimited loop: integrator convergence, invalid numerical
values, and the cycle cap still stop the run. Inspect both trajectories and
endpoint connectivity before accepting it.

Command form:

```bash
mlmm irc -i TS_STRUCTURE --parm PARM7 --model-pdb ML_REGION [options]
```

`mlmm irc --help` shows core options; `mlmm irc --help-advanced` shows the full option list.

## Workflow

1. **Input preparation** -- Load the TS structure, Amber topology (`--parm`), and ML-region definition (`--model-pdb` / `--model-indices`); resolve charge and spin. Direct PDB/mmCIF input or `--ref-pdb` supplies the topology used for companion output.
2. **ML/MM calculator setup** -- Build the ML/MM calculator from `--parm` and `--model-pdb`. The `-b/--backend` option selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). The `--hessian-calc-mode` controls ML backend Hessian evaluation. When `--embedcharge` is enabled, xTB point-charge embedding (experimental) is applied to correct for MM environment effects on the ML region.
3. **Frozen-boundary TR treatment** -- `--tr-projection constrained` removes only full-system rigid motions that leave all frozen anchors fixed. Its generic effective rank is 6/3/1/0 for zero/one/two/at least three non-collinear anchors; realistic ML/MM boundaries normally have rank 0. `legacy-active` is an isolated-active comparison treatment, not the physical default.
4. **IRC integration** -- The EulerPC integrator propagates along the IRC in both directions (unless `--no-forward` or `--no-backward` disables a branch). Step size and cycle count control integration length.
5. **Output & conversion** -- Trajectories are written as XYZ. PDB companions are generated when a PDB/mmCIF reference topology is available and `--convert-files` is enabled. Bridge inputs additionally produce CIF companions with original identifiers.

## Outputs

```text
out_dir/ (default: ./result_irc/)
├─ result.json                      # Present with --out-json; includes rigid_projection provenance
├─ <prefix>irc_data.h5              # HDF5 dump written every irc.dump_every steps
├─ <prefix>finished_irc_trj.xyz     # Full IRC trajectory (XYZ/TRJ)
├─ <prefix>forward_irc_trj.xyz      # Forward path segment
├─ <prefix>backward_irc_trj.xyz     # Backward path segment
├─ <prefix>finished_irc.pdb         # PDB companion (reference topology + conversion enabled)
├─ <prefix>finished_irc.cif         # Bridge-input companion with restored IDs
├─ <prefix>forward_irc.pdb          # Forward PDB companion (same gating)
├─ <prefix>forward_irc.cif          # Forward CIF companion (bridge input)
├─ <prefix>backward_irc.pdb         # Backward PDB companion (same gating)
├─ <prefix>backward_irc.cif         # Backward CIF companion (bridge input)
├─ <prefix>forward_last.xyz         # Single-frame forward IRC endpoint (XYZ)
├─ <prefix>forward_last.pdb/.cif    # Forward endpoint companions, when available
├─ <prefix>backward_last.xyz        # Single-frame backward IRC endpoint (XYZ)
└─ <prefix>backward_last.pdb/.cif   # Backward endpoint companions, when available
```

When `irc.prefix` is non-empty, EulerPC inserts one underscore before the
filename; for example, `prefix: trial` produces
`trial_finished_irc_trj.xyz`. `result.json.files` records the normalized names.

Standalone IRC records stitched-path `first` / `last` endpoints and their
directed bond changes; it does not assign chemical reactant/product identity.
Inspect or match the endpoint structures before naming them R/P.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation. Do not hand-duplicate the exhaustive list.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Structure file (`.pdb`/`.xyz`/`_trj.xyz`/...). | Required |
| `--parm PATH` | Amber topology for the full enzyme/MM region. Required unless `calc.real_parm7` is set in YAML. | _None_ |
| `--model-pdb PATH` | PDB defining the ML region. Required when `--no-detect-layer` and no `--model-indices` are given. | _None_ |
| `--model-indices TEXT` | Comma-separated ML-region atom indices (ranges allowed, e.g. `1-10,15`). Used when `--model-pdb` is omitted. | _None_ |
| `--model-indices-one-based/--model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer/--no-detect-layer` | Detect ML/MM layers from input PDB B-factors (`B=0/10/20`). | `True` |
| `--freeze-atoms TEXT` | Comma-separated 1-based frozen-atom indices. | _None_ |
| `--tr-projection [constrained\|legacy-active]` | Rigid-mode treatment for the frozen/partial Hessian. `legacy-active` is an isolated-active comparison treatment. | `constrained` |
| `-q, --charge INT` | Net charge; overrides `calc.charge` from YAML. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g., `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1); overrides `calc.spin`. | `1` |
| `--max-cycles INT` | Max number of IRC steps; overrides `irc.max_cycles`. | `125` |
| `--step-size FLOAT` | Step length in Bohr (unweighted Cartesian); overrides `irc.step_length`. | `0.10` |
| `--root INT` | Imaginary mode index for the initial displacement; overrides `irc.root`. | `0` |
| `--forward/--no-forward` | Run the forward IRC; overrides `irc.forward`. | `True` |
| `--backward/--no-backward` | Run the backward IRC; overrides `irc.backward`. | `True` |
| `--never-stop/--no-never-stop` | Ignore energy-rise and plateau stops only. Convergence, invalid values, and `--max-cycles` remain active. | `False` |
| `-o, --out-dir PATH` | Output directory; overrides `irc.out_dir`. | `./result_irc/` |
| `--ref-pdb FILE` | Reference PDB or mmCIF topology to use when `--input` is XYZ (keeps XYZ coordinates). | _None_ |
| `--convert-files/--no-convert-files` | Toggle XYZ/TRJ to PDB/CIF companions when a reference topology is available. | `True` |
| `--hessian-calc-mode CHOICE` | How the ML backend builds the Hessian (`Analytical` or `FiniteDifference`); overrides `calc.hessian_calc_mode`. | `FiniteDifference` |
| `--workers INT` | UMA predictor workers. Values greater than 1 require `fairchem-core[extras]` and cannot be combined with `Analytical`. | `1` |
| `--workers-per-node INT` | Workers per node for the parallel UMA predictor. | _None_ |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layers/config and continue. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects (experimental). | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap/--no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in model parm7. Default: disabled (consistent with Gaussian ONIOM). | `--no-cmap` |
| `--hess-device CHOICE` | Device for initial Hessian storage and IRC operations: `auto`, `cuda`, `cpu`. Use `cpu` for large unfrozen systems. | `auto` |
| `--read-hess PATH` | Read an identified `.npz` from `mlmm freq --dump-hess`; geometry, atom order, layer selection, and active-DOF basis must match. Takes priority over cache/fresh computation. | _None_ |
| `--mm-backend [hessian_ff\|openmm]` | MM backend (analytical Hessian vs OpenMM finite-difference). | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled ($g$-factor) or fixed 1.09/1.01 Å. | `scaled` |
| `--out-json/--no-out-json` | Write machine-readable `result.json` to `out_dir`. | `False` |
| `--dry-run/--no-dry-run` | Validate and print execution plan without running IRC. Shown in `--help-advanced`. | `False` |

Legacy NPZ files without geometry identity metadata are rejected. If the
structure or Hessian target layer changed, regenerate the file with `freq`;
partial-Hessian/PHVA metadata is preserved across a valid handoff.

## YAML configuration

Provide mappings with merge order **defaults < config < explicit CLI**.
Shared sections reuse [YAML Reference](yaml-reference.md) for geometry/calculator keys. For `irc`, `geom.coord_type` is forced to `cart` after YAML/CLI merging. `calc.return_partial_hessian` is forced to `true` (partial Hessian with active-DOF processing).

```yaml
geom:
 coord_type: cart                  # forced to cart for irc (YAML value ignored)
 freeze_atoms: []                  # 1-based frozen atoms merged with CLI/link detection
 tr_projection: constrained        # constrained (default) | legacy-active comparison
calc:
 charge: 0                         # net charge (CLI override)
 spin: 1                           # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7            # Amber parm7 topology
 model_pdb: ml_region.pdb          # ML-region definition
 backend: uma                      # MLIP backend: uma | orb | mace | aimnet2
 embedcharge: false                # xTB point-charge embedding correction
 uma_model: uma-s-1p2              # uma-s-1p2 | uma-m-1p1
 uma_task_name: omol                # UMA task name (UMA backend only)
 ml_device: auto                   # ML backend device selection
 hessian_calc_mode: Analytical        # override; default is FiniteDifference
 return_partial_hessian: true      # forced true for irc (partial Hessian with active-DOF processing)
irc:
 step_length: 0.1                  # integration step length (CLI: --step-size)
 max_cycles: 125                   # maximum steps along IRC (CLI: --max-cycles)
 forward: true                     # propagate forward branch (CLI: --forward)
 backward: true                    # propagate backward branch (CLI: --backward)
 never_stop: false                 # ignore energy-rise/plateau stops only
```

Full schema (every `irc` key and default): [YAML Reference](yaml-reference.md#irc-section).

## Notes

- Both branches run by default; disable one with `--no-forward` or `--no-backward` when you only need a single direction.
- For early stopping, reduce `--step-size` before enabling `--never-stop`; use the latter only after inspecting the surface.
- An all-frozen selection has no IRC direction and raises an explicit error.
- With `--out-json`, `result.json.rigid_projection` records the selected
  treatment, effective rank, initial-Hessian source, and Hessian shape.
- `legacy-active` uses the current common projection kernel and numerical rank
  handling; bitwise identity is not guaranteed for rank-degenerate geometries.

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [tsopt](tsopt.md) — Optimize the TS before running IRC
- [freq](freq.md) — Verify the TS candidate has one imaginary frequency; analyze IRC endpoints
- [opt](opt.md) — Optimize IRC endpoints to true minima
- [all](all.md) — End-to-end workflow that runs IRC after tsopt
- [YAML Reference](yaml-reference.md) — Full `irc` configuration options
- [Glossary](glossary.md) — Definition of IRC (Intrinsic Reaction Coordinate)
