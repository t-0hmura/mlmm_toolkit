# `tsopt`

`mlmm tsopt` refines a transition-state candidate on a layered enzyme PDB into a first-order saddle point — run it on a standalone TS guess or on the highest-energy image (HEI) extracted by [`path-search`](path-search.md). The default optimiser is **RS-I-RFO** (`--opt-mode hess`) with microiteration (`--microiter`, default on) that alternates ML 1-step RS-I-RFO and MM L-BFGS relaxation; use it as the conservative default when you can afford the Hessian work. The lighter alternative is the **Hessian-Guided Dimer** (`--opt-mode grad`), suited to a lower-cost search or quick iteration from several TS guesses, with `--ml-only-hessian-dimer` using only the ML-region Hessian for dimer orientation (faster). A surplus-imaginary-mode flatten loop (`--flatten`) sanitises extra negative modes via mass-scaled displacements after convergence. A validated TS should show **exactly one** imaginary frequency — always confirm the mode and connectivity with [`freq`](freq.md) / [`irc`](irc.md).

## Examples

The command form is `mlmm tsopt -i TS_GUESS --parm PARM7 --model-pdb ML_REGION -q CHARGE -m MULT [options]`. `mlmm tsopt --help` shows core options; `mlmm tsopt --help-advanced` shows the full option list.

Default run:

```bash
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --out-dir ./result_tsopt
```

Light mode (Dimer) with analytical Hessian:

```bash
# Light mode (Dimer) with analytical Hessian when VRAM allows
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode grad --hessian-calc-mode Analytical --out-dir ./result_tsopt_grad
```

Heavy mode (RS-I-RFO) with YAML overrides:

```bash
# Heavy mode (RS-I-RFO) with YAML overrides
mlmm tsopt -i ts_guess.pdb --parm real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode hess --config tsopt.yaml --out-dir ./result_tsopt_hess
# --dump keeps the full optimisation trajectory; --backend mace uses the MACE backend
```

## Workflow

1. **Input handling** — load the enzyme PDB, Amber topology, and ML-region definition. Resolve charge / spin. Frozen atoms from CLI and YAML are merged.
2. **ML/MM calculator setup** — build the ML/MM calculator (MLIP backend + `hessian_ff`). `-b/--backend` selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). `--hessian-calc-mode` controls whether the ML backend evaluates Hessians analytically or by finite difference. With `--embedcharge`, xTB point-charge embedding provides MM-to-ML environmental corrections.
3. **Light mode (Hessian-Guided Dimer)** — the Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian (active subspace, TR-projected). During the loose / final Dimer loops the `hessian_ff` finite-difference Hessian is disabled (`mm_fd=False`); the ML backend Hessian is embedded into the full 3N × 3N space with MM atoms zero-padded, providing a partial Hessian that still guides the Dimer direction updates. When the flatten loop is enabled (`--flatten`), the stored active Hessian is updated via Bofill using displacements and gradient differences. Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, and runs a Dimer + L-BFGS micro-segment.
4. **Heavy mode (RS-I-RFO)** — runs the RS-I-RFO optimiser with optional Hessian reference files and micro-cycle controls defined in the `rsirfo` YAML section. With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached. Once the search enters the flatten loop, a full ML/MM Hessian (including the MM finite-difference block) is computed exactly once and then updated by Bofill steps in the active subspace between Dimer segments.
5. **Mode export + conversion** — the converged imaginary mode is always written to `vib/imag_*_trj.xyz` and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimisation trajectory and final geometry are also converted to PDB via the input template when `--dump`.

## Outputs

Three artefacts are written to `result_tsopt/`: `final_geometry.pdb` (and `.xyz`) — the optimised first-order saddle point (3-layer B-factor encoding preserved for PDB); `vib/imag_*_trj.xyz` — animation of every detected imaginary mode (expect exactly one for a valid TS); and `vib/imag_*.pdb` — PDB companions of the imaginary modes (PDB inputs only).

```text
out_dir/   (default: ./result_tsopt/)
├── final_geometry.xyz                  # Always written
├── final_geometry.pdb                  # When the input was PDB
├── optimization_all_trj.xyz            # Concatenated Dimer segments (--dump)
├── optimization_all.pdb                # PDB companion (--dump, PDB input)
├── vib/
│   ├── imag_NN_±XXXX.XXcm-1_trj.xyz    # Imaginary-mode trajectory
│   └── imag_NN_±XXXX.XXcm-1.pdb        # PDB companion
└── .dimer_mode.dat                     # Dimer orientation seed (grad mode)
```

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation — do not hand-duplicate the exhaustive list.

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Starting geometry (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the whole enzyme. | Required |
| `--model-pdb PATH` | PDB containing the ML-region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer / --no-detect-layer` | Detect ML/MM layers from input PDB B-factors. | `True` |
| `-q, --charge INT` | Net charge of the ML region. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g. `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--radius-hessian` / `--hess-cutoff FLOAT` | Distance cutoff (Å) from the ML region for MM atoms to include in Hessian calculation. Applied to movable MM atoms. `0.0` means ML-only partial Hessian. | `0.0` |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) for movable MM atoms. | _None_ |
| `--hessian-calc-mode CHOICE` | ML Hessian mode: `Analytical` or `FiniteDifference`. | `FiniteDifference` |
| `--max-cycles INT` | Maximum total optimiser cycles. | `10000` |
| `--opt-mode CHOICE` | TS optimiser mode (Choice: `grad` / `hess` / `light` / `heavy` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` / `light` / `dimer` → Hessian-Guided Dimer; `hess` / `heavy` / `rsirfo` → RS-I-RFO (default, microiter-capable); `trim` → TRIM (Helgaker, non-microiter); `rsprfo` → RS-P-RFO (Banerjee, non-microiter). `trim` / `rsprfo` ignore `--microiter`. | `hess` |
| `--microiter / --no-microiter` | Microiteration: alternate ML 1-step (RS-I-RFO) + MM relaxation (L-BFGS). Only effective in `hess` mode (no-op in `--opt-mode grad`). | `True` |
| `--ml-only-hessian-dimer / --no-ml-only-hessian-dimer` | Use ML-region-only Hessian for dimer orientation in `grad` mode (faster but less accurate). | `False` |
| `--flatten / --no-flatten` | Extra-imaginary-mode flattening loop. `--flatten` uses the default iteration count (50); `--no-flatten` forces it to 0. Applies to both `--opt-mode grad` (Dimer) and `--opt-mode hess` (RS-I-RFO). | _None_ (YAML / defaults; effectively enabled with 50 iterations) |
| `--dump / --no-dump` | Write the concatenated trajectory `optimization_all_trj.xyz`. | `False` |
| `--convert-files / --no-convert-files` | Toggle XYZ / TRJ → PDB companions for PDB inputs. | `True` |
| `-o, --out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--thresh TEXT` | Convergence preset (`gau_loose` / `gau` / `gau_tight` / `gau_vtight` / `baker` / `never`). | _None_ |
| `--partial-hessian-flatten` / `--full-hessian-flatten` | Use partial Hessian (ML only) for imaginary-mode detection in the flatten loop. | `True` (partial) |
| `--active-dof-mode CHOICE` | Active DOF for final frequency analysis: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--skip-final-freq / --no-skip-final-freq` | Skip post-convergence frequency analysis and imaginary-mode flattening. Useful for large unfrozen systems where Hessian diagonalisation is expensive. TS saddle-point order will NOT be verified. | `False` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config / --no-show-config` | Print resolved config layers and continue execution. | `False` |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--precision [fp32\|fp64]` | MLIP backend precision; routed to backend-native kwarg (UMA `precision`, ORB `precision`, MACE `default_dtype`; aimnet2: fp32 no-op, fp64 rejected). | `fp32` |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding correction for MM-to-ML environmental effects. | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap / --no-cmap` | CMAP (backbone cross-map dihedral correction) in the model parm7. Disabled by default, consistent with Gaussian ONIOM. | `--no-cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend (analytical Hessian vs OpenMM finite-difference). | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled (g-factor) or fixed 1.09 / 1.01 Å. | `scaled` |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` to `out_dir`. | `False` |
| `--dry-run / --no-dry-run` | Validate inputs / config and print the execution plan without running TS optimisation (shown in `--help-advanced`). | `False` |

## YAML configuration

Settings are applied with `defaults < config < explicit CLI < override`. Shared sections reuse [YAML Reference](yaml-reference.md).
```yaml
geom:
  coord_type: cart
  freeze_atoms: []
calc:
  charge: 0
  spin: 1
mlmm:
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  backend: uma                  # uma | orb | mace | aimnet2
  hessian_calc_mode: Analytical # or FiniteDifference
opt:
  thresh: baker
  max_cycles: 10000
  out_dir: ./result_tsopt/
rsirfo:                         # --opt-mode hess
  trust_max: 0.10               # bohr; tuned for ML/MM stability near the TS
  hessian_recalc: 500           # lower (50-200) if the TS mode is lost
  track_mode_by_overlap: false  # set true if the TS mode switches root
hessian_dimer:                  # --opt-mode grad
  flatten_max_iter: 50          # 0 with --no-flatten
microiter:
  micro_thresh: null            # MM relaxation preset; null -> same as macro
```

Full schema (every section, key, and default): [YAML Reference](yaml-reference.md).

```{tip}
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimisation (e.g. when multiple imaginary frequencies are present). If TS convergence is slow or the TS mode is lost, lowering `hessian_recalc` (e.g. to 50–200) helps — more frequent exact Hessian recalculations improve robustness at the cost of additional Hessian evaluations.
```

## Notes

Active-DOF projection and mass-weighted translation / rotation removal (PHVA + TR projection) mirror `freq.py`, ensuring consistent imaginary-mode analysis and mode writing.

```{note}
`rsirfo.trust_max` defaults to 0.10 bohr for improved ML/MM stability near the TS.

The shared `opt` block also provides an **energy-plateau fallback** (`energy_plateau: true` by default, `energy_plateau_thresh: 1.0e-4` au over `energy_plateau_window: 50` steps). If the MLIP force noise floor prevents the gradient-based `thresh` preset from being reached, the optimiser still exits cleanly once the energy itself has plateaued. See [yaml-reference](yaml-reference.md#opt) for full details.

For `--microiter`, `rsirfo.thresh` controls the macro RS-I-RFO step. The MM
relaxation threshold is set with `microiter.micro_thresh`; when it is `null` or
omitted, the micro step uses the same preset as the macro step. There is no
`--micro-thresh` CLI flag.
```

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [opt](opt.md) · [freq](freq.md) · [irc](irc.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
