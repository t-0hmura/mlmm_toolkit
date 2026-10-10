# `opt` (geometry optimization)

`opt` optimizes one structure of a layered ML/MM enzyme model to a local minimum. It uses L-BFGS (`--opt-mode grad`, the default) or RFO (`--opt-mode hess`).

---

## What it is for

* **Preparing R, P, and intermediates**: relax the reactant, product, and intermediate structures before a path search or a frequency calculation, and confirm each minimum (n_imag = 0) with [`freq`](freq.md).
* **Relaxing with fixed distances**: keep chosen atom pairs at a set distance while everything else relaxes.
* **Turning IRC endpoints into R and P**: optimize the endpoints of an [`irc`](irc.md) run to the minima they lead to.
* **MM pre-relaxation**: relax the whole system on the MM force field alone (`--mm-only`) before the ML/MM optimization.

The ML region uses **UMA** (Meta) by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or **DFT**. The MM atoms use the Amber force field of `--parm7`.

---

## Examples

### 1. Basic minimization

Optimize the full system `system_layered.pdb` with its Amber topology `real.parm7` and the ML region `ml_region.pdb`, and write a summary with `--out-json`.

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --out-json --out-dir ./result_opt
```

The run converged when the console prints `[opt] Converged!` and `result_opt/result.json` has `"optimization_status": "converged"`.

### 2. Tighter threshold with trajectory

Use the `gau_tight` criteria and keep the optimization trajectory.

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --thresh gau_tight --dump --out-dir ./result_opt_tight
```

### 3. Distance restraint

Pull atoms 12 and 45 toward 2.20 Å with a weak harmonic restraint (20 eV·Å⁻²).

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --distance-restraint '[(12,45,2.20)]' --restraint-k 20.0 \
    --out-dir ./result_opt_rest
```

### 4. RFO with microiteration

Switch to RFO, which starts from an exact Hessian, with `--opt-mode hess`; microiteration is on by default.

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --opt-mode hess --out-dir ./result_opt_rfo
```

---

## How it works

1. **Building the ML/MM system**: `opt` reads the full system from `-i`, the Amber topology from `--parm7`, and the {ref}`ML region <mlmm-options>` from `--model-pdb`; the other atoms are movable or frozen MM atoms. `-q` and `-m` are the charge and the spin multiplicity of the ML region, and `--freeze-atoms` freezes more atoms.
2. **Choosing the optimizer** (`--opt-mode`): `grad` (alias `lbfgs`) runs **L-BFGS**, which uses gradients only; `hess` (alias `rfo`) runs **RFO**, which starts from an exact Hessian, updates it with TS-BFGS (the default of YAML [`rfo.hessian_update`](yaml-reference.md#rfo)), and recomputes it every 500 cycles. With `hess`, microiteration alternates one RFO step of the ML atoms and the MM parent atoms of the link atoms with an L-BFGS relaxation of the other movable MM atoms on MM forces only, as in Gaussian's microiteration.
3. **Adding distance restraints** (`--distance-restraint`): each `(i, j, target)` adds a harmonic term with force constant `--restraint-k` (eV·Å⁻²) that pulls atoms i and j toward `target` in Å; `(i, j)` keeps their starting distance. Indices are 1-based unless `--zero-based` is given.
4. **Minimizing**: the optimizer moves the structure until the convergence criteria are met or `--max-cycles` is reached. The default `--thresh gau` asks for a max force below 4.5 × 10⁻⁴ and an RMS force below 3.0 × 10⁻⁴ hartree/bohr, and a max step below 1.8 × 10⁻³ and an RMS step below 1.2 × 10⁻³ bohr, the same as Gaussian's default.
5. **Removing imaginary modes (only with `--flatten`)**: after the optimization, `opt` computes the Hessian, displaces the structure by 0.10 Å along every imaginary mode (ν < −5.00 cm⁻¹), and optimizes again, for up to 50 rounds or until no imaginary mode is left. With `--flatten`, the console prints n_imag in the line `[Imaginary modes] n=…` after each round, and `[flatten] WARNING: Remaining imaginary modes after the flatten loop: N` when modes are left after the last round.

---

## Checking convergence

How the run ended is printed on the console and recorded in `result.json` (`--out-json`):

| How it ended | `optimization_status` | Console line | `scientific_status` / exit code |
| --- | --- | --- | --- |
| Converged | `converged` | `[opt] Converged!` | `success` / 0 |
| Reached `--max-cycles` without converging | `not_converged` | `[opt] Reached max cycles (N/M).` | `failed` / 1 |
| Stopped on an energy plateau (`--stop-plateau`) | `stalled` | `[opt] Stalled (energy plateau; not converged)` | `failed` / 1 |

Each of these lines is followed by `[opt] Total cycles: N`. A `stalled` run is not converged: the energy stopped changing while the force criteria were still unmet.

Convergence gives a stationary point, not necessarily a minimum. `opt` computes no final Hessian unless `--flatten` is on, so run [`freq`](freq.md) on the final geometry and check that n_imag = 0.

---

## Output files

When the run finishes, `--out-dir` (default `./result_opt/`) contains:

```text
result_opt/
├─ final_geometry.xyz        # Final geometry (always written)
├─ final_geometry.pdb        # Same, for PDB/mmCIF input or --ref-pdb
├─ optimization_trj.xyz      # Optimization trajectory (--dump)
├─ optimization.pdb          # Same trajectory as PDB (--dump)
├─ optimization_all_trj.xyz  # All optimization steps in one trajectory (--dump)
├─ optimization_all.pdb      # Same trajectory as PDB (--dump)
├─ restart_NNN.yaml          # Optimizer state (--dump with YAML opt.dump_restart)
└─ result.json               # Summary (--out-json)
```

{ref}`mmCIF input <mmcif-input>`, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers.

* **Final geometry**: `final_geometry.*` is the optimized structure to pass to [`freq`](freq.md) or to a path search.
* **Summary**: with `--out-json`, [`result.json`](json-output.md) records `optimization_status`, the final energy `energy_hartree` (without the restraint energy), and the number of cycles `n_opt_cycles`; with microiteration, also the number of MM relaxation cycles `n_micro_cycles`.
* **Console**: the cycle table and the elapsed time.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `opt`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-l, --ligand-charge` | text | `None` | Total charge of unknown ligand residues (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-b, --backend` | text | `uma` | ML-region backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--opt-mode` | `grad` / `hess` | `grad` | Optimizer: L-BFGS / RFO (`lbfgs` and `rfo` are aliases) |
| `--microiter/--no-microiter` | flag | `True` | With `hess`, alternate RFO steps of the ML region with L-BFGS relaxations of the movable MM atoms |
| `--mm-only/--no-mm-only` | flag | `False` | Minimize the whole system on the MM force field only (`grad` only) |
| `--thresh` | preset | `gau` | Convergence criteria (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`; see below) |
| `--max-cycles` | integer | `100000` | Maximum number of optimization cycles, shared with the `--flatten` rounds |
| `--coord-type` | `cart` / `redund` / `dlc` / `tric` | `cart` | Optimization coordinates; keep `cart` for ML/MM |
| `--dump/--no-dump` | flag | `False` | Write the trajectories `optimization_trj.xyz` and `optimization_all_trj.xyz` |
| `--distance-restraint` | text | `None` | Harmonic distance restraints, inline (`'[(i,j,target_Å),...]'`) or as a YAML/JSON file; `(i,j)` keeps the starting distance |
| `--restraint-k` | float | `300` | Force constant of the distance restraints (eV·Å⁻²) |
| `--one-based/--zero-based` | flag | `--one-based` | Count `--distance-restraint` indices from 1 or from 0 |
| `--freeze-atoms` | text | `None` | Atoms to freeze (1-based, comma-separated, for example `'1,3,5'`) |
| `--hessian-cutoff` | float | `None` | Put only the movable MM atoms within this distance (Å) of the ML region into the Hessian; by default all movable MM atoms |
| `--flatten/--no-flatten` | flag | `False` | Remove imaginary modes after the optimization |
| `--reject-uphill/--no-reject-uphill` | flag | `False` | With `hess`, reject RFO steps that raise the energy by more than 1e-4 hartree and shrink the trust radius |
| `--stop-plateau/--no-stop-plateau` | flag | `False` | Stop when the energy stops changing (range below 1e-4 hartree over 50 cycles) and report `stalled` |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |

See the [generated CLI reference](reference/commands/opt.md) for every option.

The `--thresh` presets set these limits (forces in hartree/bohr, steps in bohr):

| Preset | Max force | RMS force | Max step | RMS step |
| --- | --- | --- | --- | --- |
| `gau_loose` | 2.5e-3 | 1.7e-3 | 1.0e-2 | 6.7e-3 |
| `gau` | 4.5e-4 | 3.0e-4 | 1.8e-3 | 1.2e-3 |
| `gau_tight` | 1.5e-5 | 1.0e-5 | 6.0e-5 | 4.0e-5 |
| `gau_vtight` | 2.0e-6 | 1.0e-6 | 6.0e-6 | 4.0e-6 |
| `baker` | 3.0e-4 | 2.0e-4 | 3.0e-4 | 2.0e-4 |

`baker` also requires an energy change below 1e-6 hartree between cycles. `never` never reports convergence, so the run continues to `--max-cycles`.

> **Note:** in YAML (`--config`), every key is listed under [`geom`](yaml-reference.md#geom), [`opt`](yaml-reference.md#opt), [`lbfgs`](yaml-reference.md#lbfgs), [`rfo`](yaml-reference.md#rfo), and [`microiter`](yaml-reference.md#microiter) in the YAML Reference.

---

## Notes

* **Microiteration**: with `--distance-restraint` or `--embedcharge`, `opt` runs RFO without microiteration (with `--embedcharge` because the MM-only steps leave out the embedding forces). The MM relaxation converges to the same preset as `--thresh`; YAML `microiter.micro_thresh` sets another preset.
* **`--mm-only` works only with `grad`**: with `--opt-mode hess` it stops with an error (exit code 2). The Movable-MM and Frozen-MM layers still apply.
* **Plateau stop**: `--stop-plateau` saves cycles when force noise keeps the force criteria out of reach, but a flat energy is no evidence of a stationary point. `--max-cycles` remains the real limit, and the MM relaxation of microiteration is never stopped this way. `--stop-plateau-thresh` and `--stop-plateau-window` set the energy range and the number of cycles.
* **Restraint strength**: the default force constant, 300 eV·Å⁻², holds the distance firmly; the 20 eV·Å⁻² of the example guides the structure gently toward the target.
* **`--reject-uphill` works only with `hess`**: with `grad` (L-BFGS) it is ignored.
* **One structure per run**: `-i` takes a single structure, and `.xyz` input needs `--ref-pdb` for the atom order and the layers. Extract the frame you need from a trajectory to `.xyz` first.
* **`--flatten` needs a nearly converged structure**: with more than 25 imaginary modes, `opt` skips the flatten loop and prints a warning; optimize the structure first and rerun with `--flatten`.
* **Rigid motions with frozen atoms**: `--flatten` treats rigid motions as [`freq`](freq.md#rigid-modes-with-frozen-boundaries) does, and `result.json` records them under `rigid_projection`.
* **Frozen atoms and restraints in general**: how to choose frozen atoms and restraints is described in {ref}`Frozen atoms and distance restraints <freeze-atoms-and-restraints>`.
* **Optimizer state dumps**: with `--dump`, set YAML `opt.dump_restart` to a positive integer N to write `restart_NNN.yaml` every N cycles. mlmm-toolkit does not read this file back, so rerun `opt` from the final geometry to continue a stopped calculation.
* **Model and precision**: `--backend-model` selects the model of the backend and `--precision` its precision; see the generated reference.

---

## See also

* [freq](freq.md) — check that the optimized structure is a minimum (n_imag = 0)
* [tsopt](tsopt.md) — optimize a TS (saddle point) instead of a minimum
* [irc](irc.md) — trace the reaction path from a TS to the endpoints to optimize
* [define-layer](define-layer.md) — write the ML and MM layers into the B-factors before optimizing
* [all](all.md) — the full workflow, which also optimizes the IRC endpoints
* [Troubleshooting](troubleshooting.md) — when a run fails
* [YAML Reference](yaml-reference.md) — every `opt`, `lbfgs`, `rfo`, and `microiter` setting
* [Glossary](glossary.md) — L-BFGS, RFO, and other terms
* {ref}`Exit codes <exit-codes>` — what each exit status means
