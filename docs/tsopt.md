# `tsopt`

`mlmm tsopt` refines a transition-state candidate on a layered enzyme PDB into a first-order saddle point. Run it on a standalone transition-state (TS) guess, or on the highest-energy image (HEI) extracted by [`path-search`](path-search.md).

Two optimizer families are available. The gradient family provides Hessian-Guided
Dimer (`grad`/`dimer`), while the Hessian family provides RS-I-RFO (`hess`/`rsirfo`,
the default), RS-P-RFO (`rsprfo`), and TRIM (`trim`):

- **Restricted-Step Image-function Rational Function Optimization (RS-I-RFO)** (`--opt-mode hess`) is the default and the conservative choice when you can afford the Hessian work. It runs with microiteration (`--microiter`, default on) that alternates a machine-learning (ML) 1-step RS-I-RFO move with a molecular-mechanics (MM) L-BFGS relaxation.
- **Hessian-Guided Dimer** (`--opt-mode grad`) uses initial and periodic orientation Hessians, which is more robust than a random initial direction for systems with many degrees of freedom. Add `--ml-only-hessian-dimer` to use only the ML-region Hessian for dimer orientation (faster).

`tsopt` always sets `reject_uphill: false` for its saddle-search RFO and
Dimer optimizers, including after YAML overrides. A transition-state search
must be able to raise the physical energy along its reaction mode. The
`--reject-uphill/--no-reject-uphill` toggle belongs only to minimum
optimization (`opt` and post-IRC endpoint re-optimization in `all`). The
inner MM-only relaxation in microiteration remains a minimum subproblem.

When explicitly enabled, the surplus-imaginary-mode flatten loop (`--flatten`) uses mass-scaled displacements to remove extra negative modes. Without `--flatten`, terminal exact PHVA is performed once and the terminal candidate is retained as first-order, higher-order, no-imaginary, or unavailable. First-order TS certification still requires one imaginary mode along the intended reaction coordinate and correct [`irc`](irc.md) connectivity.

### Terminal outcomes and fatal errors

| Condition | `tsopt` artifacts | Composite `all` behavior |
| --- | --- | --- |
| Convergence criteria unmet, explicit cycle limit reached, or opt-in energy plateau | Retain the final geometry and trajectory; skip terminal PHVA | Register the TS result and stop before IRC |
| Terminal PHVA fails or `--skip-final-freq` is explicit | Retain the geometry; record `failed` or `skipped` without inventing frequencies | Stop before IRC after artifact registration |
| Invalid input/geometry or an unrecoverable optimizer exception such as `ZeroStepLength` / `OptimizationError` | Follow the structured error-envelope path; only files already written are retained on a best-effort basis | Abort the stage rather than relabelling it as ordinary non-convergence |


## Building a TS candidate first

`tsopt` refines an existing candidate rather than generating one de novo.
Select the candidate-generation route according to the available structural
information, then continue through `tsopt → irc → freq` (or `mlmm all --tsopt`).

| Route | Subcommand | What it does | Use when |
| --- | --- | --- | --- |
| (a) MEP / path search | [`path-search`](path-search.md) (or [`path-opt`](path-opt.md) for one segment) | Recursive GSM/DMF minimum-energy-path search; brackets the TS between endpoints, bridges gaps between segments, and emits one TS per segment. | You have a reactant (and optionally a product or intermediates) and want the path *discovered*. |
| (b) Distance-restrained build-up | [`scan`](scan.md) | Adds a harmonic restraint `E = ½·k·(r_ij − target)²` to each reacting pair and relaxes everything else with L-BFGS, driving the reacting distance toward the barrier. | You have neither a usable second endpoint nor a TS guess — drive the reacting bond directly. |

```bash
# Route (a): discover the path, then refine its highest-energy image
mlmm path-search -i r.pdb p.pdb --parm enzyme.parm7 -l 'LIG:Q' -o result_mep

# Route (b): drive the reacting distance to build a TS candidate
mlmm scan -i r.pdb --parm enzyme.parm7 -l 'LIG:Q' \
    --scan-lists '[(1,5,1.40)]' -o result_scan
```

```{note}
There is no `opt --restraint` flag. For a restrained minimum optimization, use [`opt`](opt.md) with `--dist-freeze` and set the strength with `--bias-k`. Use [`scan`](scan.md) to drive a distance toward a TS candidate, or [`path-search`](path-search.md) to build a path.
```

## Wrong number of imaginary frequencies

A certified first-order saddle has **exactly one** imaginary mode. Inspect its
displacement and use IRC to establish the connected chemical states. Two or
more imaginary modes fail certification regardless of their magnitudes.

| Symptom | Fix |
| --- | --- |
| `n_imag = 0` (collapsed to a minimum) | Treat the run as failed. Improve the TS guess or MEP; `--flatten` only removes surplus negative modes and cannot create the missing reaction direction. |
| `n_imag > 1` | Recompute at the backend's production precision, try `--coord-type dlc`, and use `--flatten` for residual surplus modes. |
| Exactly one mode, but wrong motion | Improve the path/guess and verify connectivity by IRC; mode count alone does not identify the intended reaction. |

`--flatten` runs the surplus-imaginary-mode flattening loop (`grad`: dimer
loop; `hess`: post-RS-I-RFO); `--no-flatten` forces
`flatten_max_iter=0`. It is opt-in because it adds Hessian evaluations. When
the path itself is too coarse, rerun `all --refine-path` (or refine it with
`path-search`) before TS optimization. Recursive refinement can split a poor
path into multiple segments and substantially increase cost, so it is also off
by default.

```bash
mlmm tsopt -i ts_guess.pdb --parm enzyme.parm7 -l 'LIG:Q' -b uma \
    --precision fp64 --coord-type dlc -o result_ts
```

`--coord-type` selects the optimization coordinate system (`cart` | `redund` |
`dlc` | `tric`; default `cart`). Coordinate-system cost and convergence are
system-dependent; compare alternatives on the same seed.

```{warning}
`--coord-type dlc` needs a **Hessian-based** optimizer. On [`opt`](opt.md) with the default L-BFGS (`--opt-mode grad`) the CLI warns and falls back to `cart`; use it on `tsopt` (RFO / RS-I-RFO) or `opt --opt-mode hess`. `path-opt` / `path-search` accept only `cart` and `dlc`. `DLC + link atom` and `DLC + 3-layer frozen MM` are numerically unverified, so `cart` remains the default.
```

See [Common Error Recipes — Recipe 4](recipes-common-errors.md#recipe-4-convergence-and-post-processing-failures) for symptom-first routing of the same failure.

### Advanced MEP reference mode

`--ref-mode` is an advanced/internal handoff for Hessian-based TS optimizers,
not a normal standalone requirement. It accepts one or more atom-order-matched
Cartesian 3N candidate directions from `.npz`, `.npy`, or whitespace text (a
single vector or a 2-D candidate table). `mlmm all` supplies CPU/file-cached MEP
tangent candidates automatically for RS-I-RFO, RS-P-RFO, and TRIM. Dimer does
not consume `--ref-mode`. With `all --no-tsopt-from-mep-tan`, cache creation/use
is disabled and TSOPT selects its initial root from the initial-structure
Hessian modes.

The reference direction guides negative Hessian-root identity and overlap
tracking; it is not an initial-Hessian replacement. Terminal exact PHVA remains
authoritative for saddle order. A numerically converged higher-order stationary
point remains `optimization_status: "converged"` with
`saddle_validation: "higher_order"`; it is not a certified first-order TS.
When a validated negative root exists, `all` may continue warning-labelled
diagnostic IRC. Numerical non-convergence, zero imaginary modes, failed/skipped
PHVA, or no valid negative root stops `all` after TS artifacts are retained and
before IRC.

## Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

```{important}
Compare activation energies or free energies formed within each system
(`TS - R` or `TS - P`), then compare those barriers. Do not subtract absolute
energies between mutant and wild-type systems with different compositions.
```

Use the same electronic-structure/ML backend, MM force field, convergence
criteria, thermochemistry settings, and temperature. Define chemically
corresponding ML and movable regions, while allowing the atom count to change
where the mutation changes composition. A transferred WT layer assignment can
seed atoms that correspond unambiguously, but assign and inspect every new or
deleted atom. Determine charge and multiplicity independently for each model.
Validate each stationary point independently: a certified TS has exactly one
imaginary mode, its displacement follows the intended coordinate, and IRC
endpoints have the expected chemical identities.

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
# --dump keeps the full optimization trajectory; --backend mace uses the MACE backend
```

## Workflow

1. **Input handling** — load the enzyme PDB, Amber topology, and ML-region definition. Resolve charge / spin. Frozen atoms from CLI and YAML are merged.
2. **ML/MM calculator setup** — build the ML/MM calculator (MLIP backend + `hessian_ff`). `-b/--backend` selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). `--hessian-calc-mode` controls whether the ML backend evaluates Hessians analytically or by finite difference.
3. **Light mode (Hessian-Guided Dimer)** — the Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian in the active subspace. Its fixed constrained treatment removes only full-system rigid motions compatible with the frozen anchors. Every stored, rotated, and trial orientation has frozen Cartesian components set to zero, and every off-center force evaluation retains the central image's frozen coordinates exactly. The mechanics:
   - During the loose/final Dimer loops, the internal
     `mm_hessian_mode: none` policy intentionally uses high-level curvature
     guidance only. Outside those loops, `mm_fd: false` selects the analytical
     subtractive MM Hessian; it is not a high-level-only switch.
   - When the flatten loop is enabled (`--flatten`), the stored active Hessian is updated via Bofill using displacements and gradient differences.
   - Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, and runs a Dimer + L-BFGS micro-segment.
4. **Heavy mode (RS-I-RFO)** — runs the RS-I-RFO optimizer with optional Hessian reference files and micro-cycle controls defined in the `rsirfo` YAML section. The flatten behavior:
   - With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached.
   - Each flatten iteration recomputes a fresh ML/MM Hessian (active-coordinate block by default, or full per `--full-hessian-flatten`) for imaginary-mode detection. There is no Bofill update in this path.
5. **Mode export + conversion** — final frequency analysis writes imaginary modes to `vib/imag_*_trj.xyz` and mirrors them to `.pdb` for PDB input when conversion is enabled. The shared `freq.zero_cutoff_cm` value removes `|frequency| <= cutoff` modes before both saddle classification and trajectory output. With PDB input and conversion enabled, the final geometry is converted to PDB independently; `--dump` additionally writes and converts the optimization trajectory.

## Outputs

`result.json` separates numerical optimization from terminal exact-PHVA
classification. `optimization_status` is `converged`, `not_converged`, or
`stalled`; `saddle_validation` is `first_order`, `higher_order`,
`no_imaginary`, or `unavailable`; and `hessian_status` records whether the
terminal PHVA completed, failed, was skipped, or was unavailable. Terminal PHVA
runs only after numerical convergence; a non-converged or stalled run retains
the geometry and skips PHVA. A PHVA failure is recorded without discarding the
structure or fabricating frequencies.

A numerically converged higher-order stationary point is retained and may be
used only for warning-labelled diagnostic IRC when a validated negative root is
available. It is not first-order certification. Numerical non-convergence,
zero imaginary modes, failed/skipped PHVA, or no valid negative root causes
`all` to stop after registering TS artifacts and before IRC. Explicit
`--skip-final-freq` retains the final structure but leaves saddle order
unverified; in `all`, this therefore stops the pipeline before IRC.

Three artifacts are written to `result_tsopt/`: `final_geometry.pdb` (and `.xyz`) — the final geometry (3-layer B-factor encoding preserved for PDB); `vib/imag_*_trj.xyz` — imaginary-mode trajectories above the configured magnitude threshold; and `vib/imag_*.pdb` — their PDB companions (PDB inputs only).

```text
out_dir/   (default: ./result_tsopt/)
├── result.json                         # With --out-json; includes rigid_projection provenance
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

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
| `-i, --input PATH` | Starting geometry (PDB or XYZ). If XYZ, use `--ref-pdb` for topology. | Required |
| `--ref-pdb FILE` | Reference PDB topology when input is XYZ. | _None_ |
| `--parm PATH` | Amber parm7 topology for the whole enzyme. | Required |
| `--model-pdb PATH` | PDB containing the ML-region atoms. Optional when `--detect-layer` is enabled. | _None_ |
| `--model-indices TEXT` | Comma-separated atom indices for the ML region (ranges allowed). | _None_ |
| `--model-indices-one-based / --model-indices-zero-based` | Interpret `--model-indices` as 1-based or 0-based. | `True` (1-based) |
| `--detect-layer` | Automatically detect ML/MM layers from input PDB B-factors. | Enabled |
| `-q, --charge INT` | Net charge of the ML region. | _None_ (required unless `-l` is given) |
| `-l, --ligand-charge TEXT` | Per-resname charge mapping (e.g. `GPP:-3,SAM:1`). Derives net charge when `-q` is omitted. Requires PDB input or `--ref-pdb`. | _None_ |
| `-m, --multiplicity INT` | Spin multiplicity (2S+1) for the ML region. | `1` |
| **Active-region freezing** | | |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--radius-hessian` / `--hess-cutoff FLOAT` | Distance cutoff (Å) from the ML region for MM atoms to include in Hessian calculation. Unset includes every required movable MM atom. `0.0` requests an ML-only Hessian and must be paired with `--active-dof-mode ml-only` for final frequency validation. | _None_ |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) for movable MM atoms. | _None_ |
| **TS search & optimizer mode** | | |
| `--hessian-calc-mode CHOICE` | ML Hessian mode: `Analytical` or `FiniteDifference`. | `FiniteDifference` |
| `--ref-mode PATH` | Advanced/internal Cartesian reference candidate(s) from `.npz`, `.npy`, or whitespace text (one 3N vector or a 2-D candidate table). Guides negative Hessian-root identity/overlap; does not replace the Hessian and is unsupported by Dimer. `all` supplies it from the MEP for Hessian TS optimizers. | _None_ |
| `--max-cycles INT` | Maximum total optimizer cycles. | `100000` |
| `--opt-mode CHOICE` | TS optimizer mode (Choice: `grad` / `hess` / `light` / `heavy` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` / `light` / `dimer` → Hessian-Guided Dimer; `hess` / `heavy` / `rsirfo` → RS-I-RFO (default); `trim` → TRIM (Helgaker); `rsprfo` → RS-P-RFO (Banerjee). All three Hessian TS optimizers (`rsirfo` / `rsprfo` / `trim`) are microiter-capable. | `hess` |
| `--microiter / --no-microiter` | Microiteration: alternate a 1-step macro TS move (RS-I-RFO / RS-P-RFO / TRIM) + MM relaxation (L-BFGS). Effective in any Hessian mode (`hess` / `rsirfo` / `rsprfo` / `trim`); no-op in `--opt-mode grad` / `dimer`. | `True` |
| `--ml-only-hessian-dimer / --no-ml-only-hessian-dimer` | Use ML-region-only Hessian for dimer orientation in `grad` mode (faster but less accurate). | `False` |
| **Convergence & flatten** | | |
| `--thresh TEXT` | Convergence preset (`gau_loose` / `gau` / `gau_tight` / `gau_vtight` / `baker` / `never`). | _None_ |
| `--flatten / --no-flatten` | Extra-imaginary-mode flattening loop. `--flatten` uses the default iteration count (50); `--no-flatten` forces it to 0. Applies to both `--opt-mode grad` (Dimer) and `--opt-mode hess` (RS-I-RFO). | _None_ → disabled by default (0 iterations); `--flatten` enables it (50), and YAML/config can also enable it |
| `--partial-hessian-flatten` / `--full-hessian-flatten` | Use the active-coordinate Hessian block or the full Hessian for imaginary-mode detection in the flatten loop. | `True` (active block) |
| `--active-dof-mode CHOICE` | Active DOF for final frequency analysis: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--skip-final-freq / --no-skip-final-freq` | Skip terminal frequency/PHVA validation. The final TS candidate is retained, but saddle order and a negative IRC direction are unverified; `all` stops before IRC. | `False` |
| **Backend & compute** | | |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--precision [fp32\|fp64]` | MLIP backend precision; unset uses UMA/AIMNet2 fp32 and ORB/MACE fp64. AIMNet2 rejects fp64. | backend-specific |
| `--workers INT` | UMA predictor workers. Values greater than 1 require `fairchem-core[extras]` and cannot be combined with `Analytical`. | `1` |
| `--workers-per-node INT` | Workers per node for the parallel UMA predictor. | _None_ |
| `--allow-charge-mult-mismatch` | Skip ML-region charge/multiplicity electron-parity validation after emitting a warning. An open-shell ML region needs a matching multiplicity; use this only for an intentional nonstandard input. | off |
| `--cmap / --no-cmap` | Preserve CMAP in both REAL and MODEL MM layers. | `--cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend. Hessians use finite differences by default; set `calc.mm_fd: false` for the `hessian_ff` analytical path. | `hessian_ff` |
| `--link-atom-method [scaled\|fixed]` | Link-atom placement: scaled (g-factor) or fixed 1.09 / 1.01 Å. | `scaled` |
| **Output & config** | | |
| `--dump / --no-dump` | Write the concatenated trajectory `optimization_all_trj.xyz`. | `False` |
| `--convert-files / --no-convert-files` | Toggle XYZ / TRJ → PDB companions for PDB inputs. | `True` |
| `-o, --out-dir TEXT` | Output directory. | `./result_tsopt/` |
| `--config FILE` | Base YAML configuration applied before explicit CLI options. | _None_ |
| `--show-config / --no-show-config` | Print resolved config layers and continue execution. | `False` |
| `--out-json / --no-out-json` | Write a machine-readable `result.json` to `out_dir`. | `False` |
| `--dry-run / --no-dry-run` | Validate inputs / config and print the execution plan without running TS optimization (shown in `--help-advanced`). | `False` |

## YAML configuration

Settings are applied with `defaults < config < explicit CLI`. Shared sections reuse [YAML Reference](yaml-reference.md).
```yaml
geom:
  coord_type: cart
  freeze_atoms: []
  tr_projection: constrained      # fixed internal PHVA treatment
calc:
  model_charge: 0
  model_mult: 1
  real_parm7: real.parm7
  model_pdb: ml_region.pdb
  backend: uma                  # uma | orb | mace | aimnet2
  hessian_calc_mode: Analytical # or FiniteDifference
opt:
  thresh: baker
  max_cycles: 100000
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
Set `rsirfo.track_mode_by_overlap: true` if the TS mode switches root during optimization (e.g. when multiple imaginary frequencies are present). If TS convergence is slow or the TS mode is lost, lowering `hessian_recalc` (e.g. to 50–200) helps — more frequent exact Hessian recalculations improve robustness at the cost of additional Hessian evaluations.
```

## Notes

Frozen-boundary PHVA and mass-weighted TR treatment mirror `freq.py`. With
`constrained` (default), only full-system rigid motions that leave every frozen
anchor fixed are removed. The generic effective rank is 6/3/1/0 for
zero/one/two/at least three non-collinear anchors; realistic ML/MM boundaries
normally have rank 0. An all-frozen selection raises an explicit error.
The Dimer rebuilds this basis whenever its central image changes and applies it
to orientations and rotation forces; it does not subtract active-fragment
translations that are finite-curvature motions against the frozen boundary.

The fixed constrained rigid-mode treatment is unrelated to `--ref-mode`, which
supplies an advanced 3N MEP tangent for TS root selection and overlap tracking. A stale
non-constrained `geom.tr_projection` value fails explicitly. With `--out-json`,
`result.json.rigid_projection` records the treatment, effective rank, Hessian
source, and Hessian shape.

```{note}
`rsirfo.trust_max` defaults to 0.10 bohr for improved ML/MM stability near the TS.

The shared `opt` block also provides an **energy-plateau stop**, off by default and turned on with `--stop-plateau` (`energy_plateau_thresh: 1.0e-4` au over `energy_plateau_window: 50` steps). A plateau stops the search as `stalled` and skips terminal PHVA, as does reaching `max_cycles` without convergence. It never applies to MM micro iterations. See [yaml-reference](yaml-reference.md#opt) for details.

For `--microiter`, `rsirfo.thresh` controls the macro RS-I-RFO step. The MM
relaxation threshold is set with `microiter.micro_thresh`; when it is `null` or
omitted, the micro step uses the same preset as the macro step. There is no
`--micro-thresh` CLI flag.
```

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [opt](opt.md) · [freq](freq.md) · [irc](irc.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
