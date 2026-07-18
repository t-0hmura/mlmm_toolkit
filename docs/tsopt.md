# `tsopt`

`mlmm tsopt` refines a transition-state candidate on a layered enzyme PDB into a first-order saddle point. Run it on a standalone transition-state (TS) guess, or on the highest-energy image (HEI) extracted by [`path-search`](path-search.md).

Two optimizers are available, and you pick between them with `--opt-mode`:

- **Restricted-Step Image-function Rational Function Optimization (RS-I-RFO)** (`--opt-mode hess`) is the default and the conservative choice when you can afford the Hessian work. It runs with microiteration (`--microiter`, default on) that alternates a machine-learning (ML) 1-step RS-I-RFO move with a molecular-mechanics (MM) L-BFGS relaxation.
- **Hessian-Guided Dimer** (`--opt-mode grad`) is the lighter alternative, suited to a lower-cost search or quick iteration from several TS guesses. Add `--ml-only-hessian-dimer` to use only the ML-region Hessian for dimer orientation (faster).

After convergence, a surplus-imaginary-mode flatten loop (`--flatten`) removes extra negative modes via mass-scaled displacements. A validated TS should show **exactly one** imaginary frequency — always confirm the mode and connectivity with [`freq`](freq.md) / [`irc`](irc.md).

## Building a TS candidate first

`tsopt` refines a *candidate* — it does not find one from scratch. Pick the route that matches the information you already have, then feed the result into `tsopt → irc → freq` (or `mlmm all --tsopt`).

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
There is no `opt --restraint` flag. Plain [`opt`](opt.md) is an *un-restrained* optimizer; the restrained build-up of a TS candidate is done with [`scan`](scan.md) (drive the distance) or with [`path-search`](path-search.md) (route a).
```

## Wrong number of imaginary frequencies

A clean first-order saddle has **exactly one** dominant imaginary mode along the reaction coordinate. Two common failures are a spurious second small imaginary mode, or no dominant reaction mode at all.

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

`--coord-type` selects the optimization coordinate system (`cart` | `redund` | `dlc` | `tric`; default `cart`). `dlc` (delocalized internal coordinates) is slower but converges more robustly on torsion-rich systems and is more likely to reach a clean first-order saddle.

```{warning}
`--coord-type dlc` needs a **Hessian-based** optimizer. On [`opt`](opt.md) with the default L-BFGS (`--opt-mode grad`) the CLI warns and falls back to `cart`; use it on `tsopt` (RFO / RS-I-RFO) or `opt --opt-mode hess`. `path-opt` / `path-search` accept only `cart` and `dlc`. `DLC + link atom` and `DLC + 3-layer frozen MM` are numerically unverified, so `cart` remains the default, and is the setting used to produce the published results.
```

See [Common Error Recipes — Recipe 4](recipes-common-errors.md#recipe-4-convergence-and-post-processing-failures) for symptom-first routing of the same failure.

### Advanced MEP reference mode

`--ref-mode` reads a non-zero Cartesian 3N vector (`.npy` or whitespace text)
that identifies the MEP tangent during saddle recovery. It is an internal,
advanced input for the end-to-end workflow: `mlmm all` derives and passes it
from the MEP. Ordinary standalone `mlmm tsopt` users should omit it unless they
have constructed a matching vector in exactly the same atom order.

## Controlled mutant-vs-WT (or mechanism-vs-mechanism) comparison

```{important}
For a mutant-versus-wild-type (or mechanism-versus-mechanism) barrier comparison, **all compared models must use the same atom set — identical atom count and residues.** Otherwise the energy reference differs and the comparison is not a controlled experiment. A geometrically re-derived ML/movable/frozen partition on the mutant also produces spurious soft modes (`tsopt.n_imaginary ≥ 2`, both tiny → the IRC aborts).
```

In mlmm, preserve the wild-type ML/MM layering by **transplanting the WT B-factor layer encoding onto the mutated structure** and running with `--detect-layer`:

| Step | Action |
| --- | --- |
| 1 | Build the mutant; keep the same residue set as WT (only the mutated residue's identity differs). |
| 2 | Copy WT's per-atom B-factor layer codes onto the mutant by `(resid, atom-name)`: **ML = 0.0, MovableMM = 10.0, FrozenMM = 20.0**. |
| 3 | Run with `--detect-layer` (default on) so the *same* layer assignment is reused. |

```bash
mlmm all -i mutant_layered.pdb -l 'LIG:Q' \
    --tsopt --thermo -o result_mutant
```

| Flag | Action | Why |
| --- | --- | --- |
| `--detect-layer` | keep (default `True`) | reads the transplanted B-factors so the ML / movable / frozen layers are byte-identical to WT |
| `-c/--center`, `-r/--radius` | **omit** | geometric extraction would re-derive a *different* pocket on the mutant; omitting `-c` skips extraction and uses the structure as-is |
| `-l/--ligand-charge` | keep | a non-standard ligand's charge is not in the standard amino-acid table; `-l 'RES:Q'` auto-derives the total charge (preferred over a hardcoded `-q`) |
| `--movable-cutoff` | do not pass | it disables `--detect-layer` |

The same-atom-set principle applies equally when comparing two mechanisms on the same enzyme: keep an identical atom set across both models and vary only the reaction coordinate.

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
2. **ML/MM calculator setup** — build the ML/MM calculator (MLIP backend + `hessian_ff`). `-b/--backend` selects the MLIP (`uma`, `orb`, `mace`, or `aimnet2`; default `uma`). `--hessian-calc-mode` controls whether the ML backend evaluates Hessians analytically or by finite difference. With `--embedcharge`, xTB point-charge embedding provides MM-to-ML environmental corrections.
3. **Light mode (Hessian-Guided Dimer)** — the Dimer stage periodically refreshes the dimer direction by evaluating an exact Hessian in the active subspace. Its TR treatment follows `--tr-projection`: the default removes only full-system rigid motions compatible with the frozen anchors. Every stored, rotated, and trial orientation has frozen Cartesian components set to zero, and every off-centre force evaluation retains the central image's frozen coordinates exactly. The mechanics:
   - During the loose / final Dimer loops the `hessian_ff` finite-difference Hessian is disabled (`mm_fd=False`). The ML backend Hessian is then embedded into the full 3N × 3N space with MM atoms zero-padded, giving a partial Hessian that still guides the Dimer direction updates.
   - When the flatten loop is enabled (`--flatten`), the stored active Hessian is updated via Bofill using displacements and gradient differences.
   - Each loop estimates imaginary modes, flattens once, refreshes the dimer direction, and runs a Dimer + L-BFGS micro-segment.
4. **Heavy mode (RS-I-RFO)** — runs the RS-I-RFO optimizer with optional Hessian reference files and micro-cycle controls defined in the `rsirfo` YAML section. The flatten behaviour:
   - With `--flatten`, when more than one imaginary mode remains after convergence the workflow flattens extra modes and reruns RS-I-RFO until only one imaginary mode remains or the flatten-iteration cap is reached.
   - Each flatten iteration recomputes a fresh ML/MM Hessian (partial ML-only by default, or full per `--partial-hessian-flatten`) for imaginary-mode detection. There is no Bofill update in this path.
5. **Mode export + conversion** — the converged imaginary mode is always written to `vib/imag_*_trj.xyz` and mirrored to `.pdb` when the input was PDB and conversion is enabled. The optimization trajectory and final geometry are also converted to PDB via the input template when `--dump`.

## Outputs

With final frequency validation enabled, `result.json` reports `status:
"converged"` only when the optimizer converged and the final Hessian has exactly
one imaginary mode. It reports `not_converged` for zero or multiple modes, and
`unverified` when `--skip-final-freq` suppresses saddle-order validation.

Three artifacts are written to `result_tsopt/`: `final_geometry.pdb` (and `.xyz`) — the optimized first-order saddle point (3-layer B-factor encoding preserved for PDB); `vib/imag_*_trj.xyz` — animation of every detected imaginary mode (expect exactly one for a valid TS); and `vib/imag_*.pdb` — PDB companions of the imaginary modes (PDB inputs only).

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

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation — do not hand-duplicate the exhaustive list.

| Option | Description | Default |
| --- | --- | --- |
| **Input & charge** | | |
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
| **Active-region freezing** | | |
| `--freeze-atoms TEXT` | Comma-separated 1-based indices to freeze (merged with YAML `geom.freeze_atoms`). | _None_ |
| `--tr-projection [constrained\|legacy-active]` | TR treatment for Cartesian PHVA, Dimer refresh, flattening, and final saddle validation. `legacy-active` is an isolated-active comparison treatment. | `constrained` |
| `--radius-hessian` / `--hess-cutoff FLOAT` | Distance cutoff (Å) from the ML region for MM atoms to include in Hessian calculation. Applied to movable MM atoms. `0.0` means ML-only partial Hessian. | `0.0` |
| `--movable-cutoff FLOAT` | Distance cutoff (Å) for movable MM atoms. | _None_ |
| **TS search & optimizer mode** | | |
| `--hessian-calc-mode CHOICE` | ML Hessian mode: `Analytical` or `FiniteDifference`. | `FiniteDifference` |
| `--ref-mode PATH` | Advanced Cartesian 3N path-direction hint. `all` supplies it from the MEP; ordinary standalone `tsopt` runs omit it. | _None_ |
| `--max-cycles INT` | Maximum total optimizer cycles. | `10000` |
| `--opt-mode CHOICE` | TS optimizer mode (Choice: `grad` / `hess` / `light` / `heavy` / `dimer` / `rsirfo` / `trim` / `rsprfo`). `grad` / `light` / `dimer` → Hessian-Guided Dimer; `hess` / `heavy` / `rsirfo` → RS-I-RFO (default); `trim` → TRIM (Helgaker); `rsprfo` → RS-P-RFO (Banerjee). All three Hessian TS optimizers (`rsirfo` / `rsprfo` / `trim`) are microiter-capable. | `hess` |
| `--microiter / --no-microiter` | Microiteration: alternate a 1-step macro TS move (RS-I-RFO / RS-P-RFO / TRIM) + MM relaxation (L-BFGS). Effective in any Hessian mode (`hess` / `rsirfo` / `rsprfo` / `trim`); no-op in `--opt-mode grad` / `dimer`. | `True` |
| `--ml-only-hessian-dimer / --no-ml-only-hessian-dimer` | Use ML-region-only Hessian for dimer orientation in `grad` mode (faster but less accurate). | `False` |
| **Convergence & flatten** | | |
| `--thresh TEXT` | Convergence preset (`gau_loose` / `gau` / `gau_tight` / `gau_vtight` / `baker` / `never`). | _None_ |
| `--flatten / --no-flatten` | Extra-imaginary-mode flattening loop. `--flatten` uses the default iteration count (50); `--no-flatten` forces it to 0. Applies to both `--opt-mode grad` (Dimer) and `--opt-mode hess` (RS-I-RFO). | _None_ → disabled by default (0 iterations); `--flatten` enables it (50), and YAML/config can also enable it |
| `--partial-hessian-flatten` / `--full-hessian-flatten` | Use partial Hessian (ML only) for imaginary-mode detection in the flatten loop. | `True` (partial) |
| `--active-dof-mode CHOICE` | Active DOF for final frequency analysis: `all`, `ml-only`, `partial`, `unfrozen`. | `partial` |
| `--skip-final-freq / --no-skip-final-freq` | Skip post-convergence frequency analysis and imaginary-mode flattening. Useful for large unfrozen systems where Hessian diagonalization is expensive. TS saddle-point order will NOT be verified. | `False` |
| **Backend & compute** | | |
| `-b, --backend CHOICE` | MLIP backend for the ML region: `uma` (default), `orb`, `mace`, `aimnet2`. | `uma` |
| `--precision [fp32\|fp64]` | MLIP backend precision; unset uses UMA/AIMNet2 fp32 and ORB/MACE fp64. AIMNet2 rejects fp64. | backend-specific |
| `--workers INT` | UMA predictor workers. Values greater than 1 require `fairchem-core[extras]` and cannot be combined with `Analytical`. | `1` |
| `--workers-per-node INT` | Workers per node for the parallel UMA predictor. | _None_ |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding correction for MM-to-ML environmental effects (experimental). | `False` |
| `--embedcharge-cutoff FLOAT` | Cutoff radius (Å) for embed-charge MM atoms. | `12.0` |
| `--cmap / --no-cmap` | CMAP (backbone cross-map dihedral correction) in the model parm7. Disabled by default, consistent with Gaussian ONIOM. | `--no-cmap` |
| `--mm-backend [hessian_ff\|openmm]` | MM backend (analytical Hessian vs OpenMM finite-difference). | `hessian_ff` |
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
  tr_projection: constrained      # constrained (default) | legacy-active comparison
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

`--tr-projection` is unrelated to `--ref-mode`: the former controls
frozen-boundary rigid-mode projection, while the latter supplies an advanced 3N
MEP tangent for saddle recovery. `legacy-active` is an isolated-active
comparison treatment using the current common kernel and numerical rank
handling; bitwise identity is not guaranteed for rank-degenerate cases. With
`--out-json`, `result.json.rigid_projection` records the treatment, effective
rank, Hessian source, and Hessian shape.

```{note}
`rsirfo.trust_max` defaults to 0.10 bohr for improved ML/MM stability near the TS.

The shared `opt` block also provides an **energy-plateau fallback** (`energy_plateau: true` by default, `energy_plateau_thresh: 1.0e-4` au over `energy_plateau_window: 50` steps). If the MLIP force noise floor prevents the gradient-based `thresh` preset from being reached, a plateau stops the search instead of spending the remaining cycles, and reports `status: "stalled"` — a distinct non-converged outcome, never `converged`. The terminal exact Hessian still runs, so the saddle diagnosis is reported either way; the flatten/retry loop does not, because a stalled root is not a validated TS mode. See [yaml-reference](yaml-reference.md#opt) for full details.

For `--microiter`, `rsirfo.thresh` controls the macro RS-I-RFO step. The MM
relaxation threshold is set with `microiter.micro_thresh`; when it is `null` or
omitted, the micro step uses the same preset as the macro step. There is no
`--micro-thresh` CLI flag.
```

## See Also

[Common Error Recipes](recipes-common-errors.md) · [Troubleshooting](troubleshooting.md) · [path-search](path-search.md) · [opt](opt.md) · [freq](freq.md) · [irc](irc.md) · [all](all.md) · [YAML Reference](yaml-reference.md) · [Glossary](glossary.md).
