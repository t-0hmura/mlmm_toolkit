# `freq` (vibrational analysis and thermochemistry)

## Overview

`freq` computes **harmonic vibrational frequencies** and **thermochemical corrections** such as the zero-point energy (ZPE), enthalpy, and Gibbs free energy for a layered ML/MM enzyme model.

### What it is for

* **Checking a stationary point**: confirm that an optimized structure is a minimum (no imaginary frequency, n_imag = 0) or a transition state (TS: exactly one, n_imag = 1).
* **Thermochemistry**: free energies and other thermodynamic quantities from the QRRHO (quasi-rigid-rotor harmonic oscillator) model.
* **Seeing the modes**: atomic-displacement animations of the imaginary mode or any other mode, as `.xyz` / `.pdb` / `.cif` files of the whole enzyme.

The ML region uses **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md), by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or DFT. The MM region uses the Amber parameters in `--parm7`.

---

## Examples

In these examples, `pocket.pdb` is the full system that matches `real.parm7`, and `ml_region.pdb` defines the ML region (see [Building the ML region and layers](model-setup.md)).

### 1. Minimal run (explicit charge and multiplicity)

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_freq
```

### 2. Extra frozen atoms and the detailed thermochemistry file

`--freeze-atoms` freezes more atoms on top of the frozen MM layer, and `--dump` also writes the detailed thermochemistry file `thermoanalysis.yaml`.

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --freeze-atoms "1,3,5,7" --dump --out-dir ./result_freq_phva
```

### 3. Analytical Hessian

Use this to avoid the step-size error of finite differences. It uses more GPU memory, so test it on your system first.

```bash
mlmm freq -i pocket.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --hessian-calc-mode Analytical --out-dir ./result_freq_analytical
```

---

## How it works

1. **Reading the layers and freezing (PHVA)**:
`freq` takes the ML region from `--model-pdb` and the movable and frozen MM layers from the B-factors of the input PDB. `--active-dof-mode` chooses the atoms in the vibrational analysis; the default `partial` takes the ML region and the movable MM atoms. The frozen MM layer and the atoms given with `--freeze-atoms` stay fixed.
2. **Hessian**:
`--hessian-calc-mode` selects `FiniteDifference` (finite differences, the default) or `Analytical` for the ML region. `--hess-device` chooses where the computed Hessian is held and diagonalized.
3. **Thermochemistry (QRRHO)**:
From the frequencies, the QRRHO model with a rotor cutoff of 100 cm⁻¹ gives the Gibbs free-energy correction `G_corr`. The thermochemistry summary on the console prints `G_corr` and the Gibbs energy `E + G_corr = G` in Hartree (E: electronic energy); `--dump` also writes them to `thermoanalysis.yaml`, and `--out-json` to `thermochemistry` in `result.json`. With frozen atoms, the vibrational terms come from the PHVA frequencies. Besides the vibrational terms from the positive frequencies, G always includes the translational and rotational terms of the whole structure; the point group and rotational symmetry number are detected from the structure, and YAML `thermo.symmetry_number` overrides the detected number.
4. **Writing the modes**:
Mode animations (trajectory files) are written starting from the imaginary or lowest modes; `--max-write` sets how many and `--sort` the order.

### Rigid modes with frozen boundaries

Without frozen atoms, `freq` removes the six rigid motions (three translations and three rotations), which are not vibrations, before reporting frequencies. With frozen atoms, it removes only the rigid motions that keep every frozen atom in place. With three or more frozen atoms that do not lie on one line, the normal case for an ML/MM model with a frozen MM layer, nothing is removed and every vibrational mode of the movable atoms is kept. With one frozen atom, three motions are removed (rotations about that atom); with two, one is removed (rotation about the axis through them).

`irc`, the TS frequency check and the Dimer direction in `tsopt`, and `--flatten` (removing extra imaginary modes) in `opt` and `tsopt` treat rigid motions the same way. With `--out-json`, `result.json` records the number of removed motions and the Hessian used under `rigid_projection`; see [JSON Output Reference](json-output.md#rigid-projection-provenance).

---

## Reading the frequencies

How `frequencies_cm-1.txt` and the JSON record treat each case:

| Item | Value / behavior | Meaning |
| --- | --- | --- |
| **Imaginary modes** | Negative values (ν < 0 cm⁻¹) | An imaginary mode is listed as a negative frequency. |
| **Imaginary threshold** | ν < −5.00 cm⁻¹ | Such a mode counts as imaginary (n_imag; the JSON field is `n_imaginary`). The cutoff is YAML `freq.zero_cutoff_cm` (default `5.0`). |
| **Tiny negative modes** | −5.00 ≤ ν < 0 cm⁻¹ | Small negative modes from numerical noise do not count toward n_imag. `n_negative_modes` counts every negative frequency, including these. |
| **Thermochemistry** | No inversion, no floor | Imaginary modes are not flipped and small positive modes are not raised. QRRHO uses only the positive modes, so imaginary modes are left out of ZPE and G. |

---

## Output files

`freq` writes these files to `--out-dir` (default `./result_freq/`):

```text
result_freq/
├─ frequencies_cm-1.txt          # All frequencies (cm⁻¹)
├─ mode_0001_-385.20cm-1_trj.xyz # Animation of each mode (XYZ)
├─ mode_0001_-385.20cm-1.pdb     # Same animation as PDB
├─ mode_0001_-385.20cm-1.cif     # Same animation as mmCIF (mmCIF or very large PDB input)
├─ thermoanalysis.yaml           # Detailed thermochemistry (with --dump)
└─ result.json                   # Summary (with --out-json)
```

* **Minimum or TS?** The thermochemistry summary on the console prints n_imag as `Number of Imaginary Freq = N`; `freq` does not judge it, so `scientific_status` in `result.json` is `success` whatever n_imag is. A minimum has n_imag = 0. A successful TS optimization gives one imaginary mode along the reaction coordinate: the top of `frequencies_cm-1.txt` should hold **exactly one** clear imaginary frequency (a negative value), and every value after it should be positive or within the tolerance. A TS then goes to [`irc`](irc.md). If a structure meant to be a minimum has imaginary modes, optimize it again with [`opt`](opt.md) `--flatten`; if a TS has none or several, see {ref}`When the TS search fails <ts-search-fails>`.
* **Watching the motion**: open `mode_*_trj.xyz` or `.pdb` in PyMOL, VMD, or another viewer to animate the vibration.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `freq`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure matching `--parm7` (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Total charge of unknown ligands or a charge per residue name (e.g. `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB input or `--ref-pdb`) |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `-b, --backend` | text | `uma` | Backend for the ML region (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `--hessian-calc-mode` | `FiniteDifference` / `Analytical` | `FiniteDifference` | How the ML-region Hessian is computed (finite differences / analytical) |
| `--hess-device` | `auto` / `cuda` / `cpu` | `auto` | Device that holds and diagonalizes the computed Hessian; `cpu` moves it off the GPU first |
| `--active-dof-mode` | `all` / `ml-only` / `partial` / `unfrozen` | `partial` | Atoms in the analysis: all atoms / ML region only / ML region and movable MM atoms / every atom outside the frozen layer |
| `--freeze-atoms` | text | `None` | Extra atoms to freeze (1-based, comma-separated, e.g. `'1,3,5'`) |
| `--read-hess` | path | `None` | Read the Hessian from a `.npy` file (for example one saved by `freq` or `tsopt --dump-hess`) instead of computing it |
| `--dump-hess` | path | `None` | Save the Hessian as a `.npy` file for `--read-hess` in `freq`, `tsopt`, or `irc` |
| `--max-write` | integer | `10` | Maximum number of mode animations to write |
| `--sort` | `value` / `abs` | `value` | Order of the modes (by value / by absolute value) |
| `--temperature` | float | `298.15` | Temperature for thermochemistry (K) |
| `--pressure` | float | `1.0` | Pressure for thermochemistry (atm) |
| `--dump/--no-dump` | flag | `False` | Write the detailed thermochemistry file `thermoanalysis.yaml` |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` |

See the [generated CLI reference](reference/commands/freq.md) for every option.

> **Note:** In YAML (`--config`), the {ref}`freq <freq-section>` section sets the imaginary threshold `zero_cutoff_cm` and the number and amplitude of the written modes, and the [`thermo`](yaml-reference.md#thermo) section sets temperature and pressure.

---

## Notes

* **`freq` or `tsopt`?** `tsopt` already checks the imaginary frequencies. Run `freq` on its own when you need detailed thermochemistry (ZPE, Gibbs energy) or mode animations.
* **At least one atom must move.** If every atom is frozen there is no vibration to analyze, and `freq` stops with an error.
* **Hessian mode priority**: `--hessian-calc-mode` follows the priority default < YAML config < command line.
* **Analytical Hessian and `--uma-workers`**: with UMA, `--hessian-calc-mode Analytical` cannot run with `--uma-workers` (parallel MLIP predictor workers) above 1 and stops with an error. Use `--uma-workers 1` for an analytical Hessian ([details](backends.md#workers-and-hessian-mode)).
* **`all --thermo` keeps the thermochemistry file**: `all` builds its Gibbs energy diagram from `thermoanalysis.yaml`, so with `--thermo` its `freq` step writes this file even under `--no-dump`.
* **The `--read-hess` / `--dump-hess` file** is one NumPy array (`numpy.save`): the Cartesian Hessian in Hartree/bohr², not mass-weighted, with atoms in input order. It covers all atoms (3N × 3N) or only the atoms selected by `--active-dof-mode`. `--read-hess` checks only that the matrix is square, finite, symmetric, and one of these two sizes, so pass a Hessian computed for the same geometry, charge, multiplicity, layers, and calculator settings.

---

## See also

* [opt](opt.md) — geometry optimization to a minimum
* [tsopt](tsopt.md) — TS optimization
* [irc](irc.md) — IRC from a TS
* [dft](dft.md) — DFT single-point energies
* [all](all.md) — the full workflow: model building, path search, TS optimization, and vibrational analysis
* [YAML Reference](yaml-reference.md) — configuration file format
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* {ref}`Exit codes <exit-codes>` — what each exit status means
