# `sp` (single point)

`sp` computes the **ML/MM ONIOM energy and atomic forces** of one structure, and with `--hess` also the **Hessian** of the atoms that move. The ML region is computed with the selected backend and the rest of the enzyme with the Amber force field of `--parm7`. It runs no optimization: the geometry stays as given.

---

## What it is for

* **Check before an optimization**: confirm that the ML region, charge, and multiplicity are accepted and that the backend returns a finite energy and forces.
* **Compare backends**: evaluate the same structure and ML region with UMA, ORB, MACE, AIMNet2, or DFT (`-b dft`).
* **Reference values**: forces and Hessians as `.npy` files, and the energy in the console or `result.json`, for your own analysis.

---

## Examples

### 1. Energy and forces

Evaluate a neutral singlet ML region with the default backend (UMA). `enzyme.pdb` is the full system, `real.parm7` its Amber topology, and `ml_region.pdb` the atoms of the ML region. `-q` and `-m` are the charge and multiplicity of the ML region.

```bash
mlmm sp -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --out-json
```

The console prints `[sp] energy = … a.u.  |force|_max = … a.u./bohr`, and `result_sp/` has `forces.npy` and `result.json` with `energy_au`.

### 2. Add the Hessian

`--hess` also computes the Hessian of the moving atoms.

```bash
mlmm sp -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --hess
```

---

## How it works

1. **Building the ML/MM system**:
`sp` reads the full system from `-i`, the Amber topology from `--parm7`, and the ML region from `--model-pdb`, `--model-indices`, or the B-factors of the input. The charge comes from `-q`, or from `-l` with PDB/mmCIF input. The Frozen-MM layer and the atoms given with `--freeze-atoms` are frozen.
2. **Energy and forces**:
The backend computes the ML region and the force field the MM atoms, once at the input geometry, and the two are combined into the ONIOM energy and forces. `sp` prints the energy and the largest force component and saves the forces to `forces.npy`; frozen atoms get zero force.
3. **Hessian (with `--hess`)**:
The Hessian covers the ML region and the movable MM atoms, without frozen atoms. `--hessian-calc-mode FiniteDifference` (the default) differentiates the forces numerically; `Analytical` uses the analytical Hessian of UMA, ORB, MACE, or AIMNet2 for the ML region and cannot run with `--uma-workers` (parallel MLIP predictor workers) above 1. The MM part uses finite differences by default; YAML `calc.mm_fd: false` selects the analytical MM Hessian of `hessian_ff`.

---

## Output files

`sp` writes these files to `--out-dir`:

| File | Contents | Written |
| --- | --- | --- |
| `forces.npy` | ONIOM forces as an `(N, 3)` array over all atoms of the full system, in Hartree/bohr | Always |
| `hessian.npy` | ONIOM Hessian without mass weighting (Hartree/bohr²): `(3M, 3M)` for the M atoms of the Hessian, in input order | With `--hess` |
| `result.json` | Energy (`energy_au`), backend, model, charge, multiplicity, ML region (source and atom count), paths to the `.npy` files, elapsed time | With `--out-json` |
| `summary.json` | Copy of `result.json`; read `result.json` | With `--out-json` |

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `sp`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure (`.pdb`, `.cif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Per-residue formal charges (e.g. `'SAM:1,GPP:-3'`) or one total ligand charge, used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-b, --backend` | text | `uma` | ML-region backend (`uma`, `orb`, `mace`, `aimnet2`, `dft`); for the `-b dft` settings see [Refine an MLIP TS with DFT](dft-backend.md) |
| `--hess/--no-hess` | flag | `False` | Also compute the Hessian and write `hessian.npy` |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | Hessian method (finite difference / analytical); used with `--hess` |
| `--hessian-cutoff` | float | `None` | Put only the movable MM atoms within this distance (Å) of the ML region into the Hessian; by default all movable MM atoms |
| `--freeze-atoms` | text | `None` | Atoms to freeze (1-based, comma-separated, e.g. `'1,3,5'`) |
| `--embedcharge/--no-embedcharge` | flag | `False` | Electrostatic embedding of the MM point charges (xTB correction for an MLIP, PySCF point charges for `-b dft`) |
| `-o, --out-dir` | path | `./result_sp/` | Output directory |
| `--out-json/--no-out-json` | flag | `False` | Write `result.json` and `summary.json` |

See the [generated CLI reference](reference/commands/sp.md) for every option.

> **Note:** In YAML (`--config`), `calc` sets the backend and `geom.freeze_atoms` adds frozen atoms (1-based), merged with `--freeze-atoms`.

---

## Notes

* **Failed run**: a failed run prints a one-line `Error: …`, such as `ML region electron count inconsistent`, or `Unhandled error during single-point:` with a traceback, and exits with a nonzero code.
* **Energy looks wrong**: if the energy is finite but looks wrong, re-check the ML region and its charge and multiplicity ({ref}`Charge / spin <charge--spin>`).
* **Frozen atoms**: indices are 1-based, and frozen atoms get zero force. The Frozen-MM layer is frozen as well.
* **Atomic charges**: `sp -b dft` gives the ML(DFT)/MM energy and forces only. For Mulliken, meta-Löwdin, and IAO charges of the ML region, use [`dft`](dft.md).
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [opt](opt.md) — optimize the structure
* [tsopt](tsopt.md) — optimize a transition-state (TS) candidate
* [freq](freq.md) — vibrational analysis and thermochemistry
* [dft](dft.md) — DFT single point of the ML region with atomic charges
* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` settings (`--func-basis`, `--dft-engine`) and GPU memory
* [MLIP Backends](backends.md) — choosing a backend, precision, and workers
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
