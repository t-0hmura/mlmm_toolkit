# `opt`

Optimizes a single layered enzyme PDB (or XYZ + `--ref-pdb`) to a local minimum with the ML/MM calculator (MLIP region + movable MM shell + frozen outer environment). Use it to relax a full-system layered structure: `--opt-mode grad` (default) runs L-BFGS, `--opt-mode hess` runs RFOptimizer (RFO), `--flatten` flattens imaginary modes after optimization, and `--mm-only` does a cheap MM pre-relaxation before ML/MM ONIOM optimization. Microiteration (`--microiter`, default on) relaxes the movable-MM shell in `hess` mode.

## Examples

Command form:

```bash
mlmm opt -i INPUT --parm PARM7 --model-pdb ML_REGION -q CHARGE [options]
```

`mlmm opt --help` shows core options; `mlmm opt --help-advanced` shows the full option list.

Minimal L-BFGS optimization (grad mode, default):

```bash
# Minimal L-BFGS optimization (grad mode, default)
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --out-dir ./result_opt
```

Tighten convergence and keep an optimization trajectory:

```bash
# Tighten convergence and keep an optimization trajectory
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --thresh gau_tight --dump --out-dir ./result_opt_tight
# add one harmonic distance restraint: --dist-freeze "[(12,45,2.20)]" --bias-k 20.0
```

Switch to heavy mode (RFO):

```bash
# Switch to heavy mode (RFO)
mlmm opt -i pocket.pdb --parm real.parm7 --model-pdb ml_region.pdb \
 -q 0 --opt-mode hess --out-dir ./result_opt_rfo
# use the ORB backend instead of the default: --backend orb
```

## Workflow

1. **Input handling** -- The tool accepts `-i/--input` as a PDB or XYZ file (use `--ref-pdb` with XYZ inputs). The optimizer reads coordinates from this PDB via `pysisyphus.helpers.geom_loader`. ML/MM layer definitions come from `--model-pdb`, `--model-indices`, or `--detect-layer` (B-factor encoding: B=0 ML, B=10 Movable-MM, B=20 Frozen).
2. **ML/MM calculator setup** -- Build the ML/MM calculator (MLIP backend + hessian_ff). The `-b/--backend` option selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). `--parm` provides Amber MM topology; `--model-pdb` defines the ML region. When `--embedcharge` is enabled, xTB point-charge embedding is applied to correct for MM-to-ML environmental electrostatic effects.
3. **Optimization** -- The optimizer runs in the selected `--opt-mode` (`grad` = L-BFGS, `hess` = RFOptimizer); see the CLI options table for the accepted aliases.
   - `--flatten` enables post-optimization flattening of imaginary modes. All detected imaginary modes are flattened each iteration until none remain or the internal loop cap is reached.
4. **Restraints** -- Optional harmonic distance restraints via `--dist-freeze` / `--bias-k` (see CLI options).
5. **Dumping & conversion** -- `--dump` writes `optimization_trj.xyz`; when conversion is enabled, trajectories are mirrored to `.pdb` for PDB inputs (with B-factor annotations). `opt.dump_restart` can emit restart YAML snapshots.
6. **Exit codes** -- `0` success, `2` zero step (step norm < `min_step_norm`), `3` optimizer failure, `130` keyboard interrupt, `1` unexpected error.

## Outputs

```text
out_dir/ (default: ./result_opt/)
├─ final_geometry.xyz          # Always written
├─ final_geometry.pdb          # Only when the input was a PDB and conversion is enabled (B-factors annotated)
├─ optimization_trj.xyz        # Only if dumping is enabled
├─ optimization.pdb            # PDB conversion of the trajectory (PDB inputs, conversion enabled)
├─ optimization_all_trj.xyz    # Concatenated full trajectory (when --dump)
├─ optimization_all.pdb        # PDB companion of the full trajectory (PDB inputs, when --dump)
└─ restart*.yml                # Optional restarts when opt.dump_restart is set
```

Console output prints resolved configuration blocks (`geom`, `calc`, `opt`, `lbfgs`), progress every `print_every` cycles, and a final wall-clock time summary.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation. Default values shown are used when the option is not specified.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input structure accepted by `geom_loader` (`.pdb`, `.xyz`, `_trj.xyz`). | Required |
| `--ref-pdb PATH` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the full enzyme. | Required |
| `--model-pdb PATH` | PDB defining the ML region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed, e.g. `1-5`). Alternative to `--model-pdb`. | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Index convention for `--model-indices`. | 1-based |
| `--detect-layer / --no-detect-layer` | Auto-detect ML/MM layers from B-factors (B=0 ML, B=10 Movable-MM, B=20 Frozen). | Enabled |
| `-q, --charge INT` | Charge of the ML region. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g., `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1). | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze. | _None_ |
| `--radius-freeze FLOAT` | Distance cutoff (Å) from ML region for movable MM atoms. Atoms beyond this are frozen. Providing this disables `--detect-layer`. Alias: `--movable-cutoff`. | _None_ |
| `--radius-partial-hessian, --hess-cutoff FLOAT` | Distance cutoff (Å) from ML region for MM atoms included in Hessian calculation. Combinable with `--detect-layer`. | _None_ |
| `--mm-backend [hessian_ff\|openmm]` | MM backend (analytical Hessian vs OpenMM finite-difference). | `hessian_ff` |
| `--mm-only / --no-mm-only` | Skip the MLIP component and minimize on the MM force field only. Layers are still honored via B-factor / `--radius-freeze`; only `--opt-mode grad` is supported in this mode and microiteration is disabled automatically. Useful for a cheap MM pre-relaxation before ML/MM ONIOM optimization. | `False` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled ($g$-factor) or fixed 1.09/1.01 Å. | `scaled` |
| `--out-json/--no-out-json` | Write machine-readable `result.json` to `out_dir`. | `False` |
| `--dist-freeze TEXT` | Python-literal `(i, j, target_A)` tuples for harmonic restraints (inline literal or YAML/JSON file path); omit `target_A` to restrain the starting distance. | _None_ |
| `--one-based / --zero-based` | Index convention for `--dist-freeze`. | 1-based |
| `--bias-k FLOAT` | Harmonic bias strength (eV/Å²). | `300.0` |
| `--max-cycles INT` | Hard limit on optimization iterations. | `10000` |
| `--opt-mode [grad\|hess\|light\|heavy\|lbfgs\|rfo]` | Optimizer mode: `grad` (LBFGS) or `hess` (RFO). Aliases `light`/`heavy` and `lbfgs`/`rfo` accepted. | `grad` |
| `--microiter/--no-microiter` | Microiteration: alternate ML 1-step (RFO) + MM relaxation (LBFGS). Only effective in `hess` mode (no-op in `--opt-mode grad`). | `True` |
| `--flatten/--no-flatten` | Enable/disable the post-optimization imaginary-mode flatten loop. | `False` |
| `--dump/--no-dump` | Emit trajectory dumps (`optimization_trj.xyz`, `optimization_all_trj.xyz`). | `False` |
| `--convert-files/--no-convert-files` | Enable or disable XYZ/TRJ to PDB companions for PDB inputs. | `True` |
| `-o, --out-dir TEXT` | Output directory for all files. | `./result_opt/` |
| `--thresh TEXT` | Override convergence preset (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`). | _None_ (`gau` applied internally) |
| `--config FILE` | Base YAML configuration file. | _None_ |
| `--show-config/--no-show-config` | Print resolved YAML layer information before execution. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma`, `orb`, `mace`, `aimnet2`. | `uma` |
| `--embedcharge/--no-embedcharge` | Enable xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap/--no-cmap` | Enable CMAP (backbone cross-map dihedral correction) in model parm7. Default: disabled (consistent with Gaussian ONIOM). | `--no-cmap` |
| `--dry-run/--no-dry-run` | Validate options and print execution plan without running optimization. Shown in `--help-advanced`. | `False` |

### Convergence threshold presets

Forces in Hartree/bohr, steps in bohr.

| Preset | Purpose | max\|F\| | RMS(F) | max\|step\| | RMS(step) |
| --- | --- | --- | --- | --- | --- |
| `gau_loose` | Loose/quick pre-optimization; rough path searches | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | Standard Gaussian-like tightness for routine work | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | Tighter; better structures / freq / TS refinement | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | Very tight; benchmarking/high-precision final structures | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | Baker-style rule (converged if `max\|F\| < 3e-4` **and** `\|dE\| < 1e-6 or max\|step\| < 3e-4`) | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

## YAML configuration

Settings are applied with **defaults < config < explicit CLI < override**. The accepted sections are `geom` (`coord_type`, `freeze_atoms`), `calc` / `mlmm` (ML/MM calculator: backends, devices, Hessian mode, embedding), `opt` (shared optimizer controls), and the optimizer-specific `lbfgs` / `rfo` sections.

```yaml
geom:
 coord_type: cart               # cartesian vs dlc internals (dlc needs --opt-mode hess; grad/L-BFGS falls back to cart)
 freeze_atoms: []               # 1-based frozen atoms
calc:
 charge: 0                      # net charge
 spin: 1                        # spin multiplicity 2S+1
mlmm:
 real_parm7: real.parm7         # Amber parm7 topology
 model_pdb: ml_region.pdb       # ML region definition
 backend: uma                   # uma | orb | mace | aimnet2
 hessian_calc_mode: Analytical  # or FiniteDifference
opt:
 thresh: gau                    # convergence preset
 max_cycles: 10000              # optimizer cycle cap
 out_dir: ./result_opt/         # output directory
```

### `microiter`

Used only when `--microiter` is active with `--opt-mode hess`. `micro_thresh` sets the L-BFGS convergence preset for the MM relaxation step. When `null` or omitted, the micro step uses the same preset as the macro optimizer (`--thresh` / `opt.thresh`). There is no `--micro-thresh` CLI flag; set this in YAML.

Full schema (every section, key, and default): [YAML Reference](yaml-reference.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [tsopt](tsopt.md) — Optimize transition states (saddle points) instead of minima
- [freq](freq.md) — Vibrational analysis to confirm a minimum was reached
- [all](all.md) — End-to-end workflow that pre-optimizes endpoints
- [YAML Reference](yaml-reference.md) — Full `opt`, `lbfgs`, `rfo` configuration options
- [Glossary](glossary.md) — Definitions of L-BFGS, RFO
