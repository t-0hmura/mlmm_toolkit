# `dft` (DFT single point)

`dft` runs a **DFT (density functional theory) single point on the ML region** of one ML/MM structure with GPU4PySCF (GPU) or PySCF (CPU), and combines it with the MM energies into the **ML(DFT)/MM total energy**, `E_total = E_REAL_low + E_ML(DFT) - E_MODEL_low`. It also reports the **atomic charges** of the ML region. It computes energies only, no forces.

---

## What it is for

* **DFT energies on ML/MM geometries**: single points on the reactant (R), transition state (TS), and product (P) optimized with an MLIP for the ML region.
* **Charge distribution**: per-atom charges of the ML region, and spin densities for open shells.
* **Electrostatics of the protein**: `--embedcharge` puts the MM point charges into the DFT Hamiltonian.

---

## Examples

### 1. GPU single point

Compute the energy and charges of a neutral singlet ML region on the GPU. `enzyme.pdb` is the full system, `real.parm7` its Amber topology, and `ml_region.pdb` the atoms computed with DFT. `-q` and `-m` are the charge and multiplicity of the ML region.

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --out-dir ./result_dft
```

The console prints `E_DFT (Hartree): …` and `E_total ML(dft)/MM (Hartree): …`, and `result_dft/result.yaml` has `energy.converged: true`.

### 2. Tighter SCF and a larger basis

Tighten the SCF (self-consistent field) and use a larger basis.

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-tzvpd' --scf-tol 1e-10 --scf-max-cycles 200 \
    --out-dir ./result_dft_tight
```

### 3. CPU only

Run with CPU PySCF on a machine without a GPU.

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 \
    --dft-engine cpu --out-dir ./result_dft_cpu
```

### 4. ML-region charge from ligand charges

Without `-q`, `-l` gives the formal charges of the ligands, and `dft` adds the charges of the amino-acid residues and ions in the ML region to get its charge; the console prints the breakdown.

```bash
mlmm dft -i enzyme.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -l 'SAM:1,GPP:-3' -m 1 --out-dir ./result_dft_ligand
```

---

## How it works

1. **Building the ML region**:
`dft` reads the full system from `-i`, the Amber topology from `--parm7`, and the ML region from `--model-pdb`, `--model-indices`, or the B-factors of the input. For an XYZ input, `--ref-pdb` gives the PDB/mmCIF topology. Link hydrogens cap the bonds of `--parm7` that the ML/MM boundary cuts, and the ML region is saved without and with them.
2. **SCF**:
`--func-basis` sets the functional and basis; a basis whose name begins with `def2` gets the matching def2 effective core potential (ECP). `--dft-engine` selects GPU4PySCF (`gpu`, the default) or PySCF (`cpu`). A closed shell runs RKS and an open shell UKS. Low-memory mode, on by default, builds J and K directly without density fitting; on the GPU, a closed shell then uses GPU4PySCF's low-memory RKS. `--no-dft-low-memory` uses density fitting instead. With `--embedcharge`, the MM point charges of `--parm7` within `--embedcharge-cutoff` of the ML region enter the DFT Hamiltonian.
3. **ML(DFT)/MM energy**:
The DFT energy of the ML region with link hydrogens takes the place of the ML energy in the ONIOM sum: `E_REAL_low` and `E_MODEL_low` are the MM energies of the full system and of the ML region.
4. **Charges and the result file**:
After the SCF, `dft` computes Mulliken, meta-Löwdin, and IAO (intrinsic atomic orbital) charges and spin densities of the ML region and writes them with the energies (Hartree and kcal/mol) to `result.yaml`. An analysis that fails gives `null` in its column.

---

## Output files

`dft` writes these files to `--out-dir`:

```text
result_dft/
├─ ml_region_without_linkH.xyz   # ML region as selected, without link hydrogens
├─ ml_region_with_linkH.xyz      # ML region with link hydrogens, as passed to PySCF
├─ ml_region_without_linkH.pdb   # Same as PDB (PDB input with --convert-files)
├─ ml_region_with_linkH.pdb      # Same as PDB (PDB input with --convert-files)
├─ result.yaml                   # Energies, convergence, engine, per-atom charges and spin densities
├─ result.json                   # Machine-readable summary (with --out-json)
└─ summary.json                  # Copy of result.json; read result.json (with --out-json)
```

* **`energy`** in `result.yaml`: the DFT energy of the ML region (`hartree`, `kcal_per_mol`), `converged`, and the engine used (`engine`: `gpu4pyscf(rks_lowmem)`, `gpu4pyscf`, or `pyscf(cpu)`; `used_gpu`; `used_lowmem`).
* **`mlmm_energy`** in `result.yaml`: the MM energies `E_real_low_hartree` and `E_model_low_hartree` and the total `E_total_ml_dft_mm_hartree` (also in kcal/mol).
* **`charges [index, element, mulliken, lowdin, iao]`**: one row per atom of the ML region with link hydrogens; `index` starts at 0. The console prints the same table.
* **`spin_densities [index, element, mulliken, lowdin, iao]`**: the same layout. It is always written to `result.yaml`, and the console prints it only for an open shell.
* **`result.json`** holds the energies, the charges and spin densities as `mulliken`, `lowdin`, and `iao` arrays, the charge, multiplicity, functional, basis, and SCF settings; see [JSON Output Reference](json-output.md#dft).

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `dft`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system structure (`.pdb`, `.cif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` or YAML `calc.model_charge` gives it |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-l, --ligand-charge` | text | `None` | Per-residue formal charges (e.g. `'SAM:1,GPP:-3'`) or one total ligand charge, used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `--func-basis` | text | `wb97m-v/def2-svp` | Functional and basis as `FUNCTIONAL/BASIS` |
| `--scf-tol` | float | `1e-9` | SCF convergence threshold (Hartree) |
| `--scf-max-cycles` | integer | `100` | Maximum number of SCF iterations |
| `--dft-grid-level` | integer | `3` | Integration grid level (PySCF `grids.level`) |
| `--dft-engine` | `gpu` / `cpu` | `gpu` | GPU4PySCF or CPU PySCF |
| `--dft-low-memory/--no-dft-low-memory` | flag | `True` | Build J and K directly; `--no-dft-low-memory` uses density fitting |
| `--scf-stepwise-grid/--no-scf-stepwise-grid` | flag | `False` | Converge the SCF on a [coarse grid](dft-backend.md#notes) first, then on the final grid |
| `--dft-nprocs` | integer | auto | PySCF CPU threads (detected from the scheduler and the host) |
| `--dft-memory` | text | auto | PySCF host RAM limit (e.g. `64GB`); this is not GPU memory |
| `--embedcharge/--no-embedcharge` | flag | `False` | Put the MM point charges into the DFT Hamiltonian |
| `--embedcharge-cutoff` | float | `12.0` | Distance (Å) from the ML region within which MM point charges are embedded |
| `--convert-files/--no-convert-files` | flag | `True` | Also write the ML region as PDB (PDB input only) |
| `-o, --out-dir` | path | `./result_dft/` | Output directory |

See the [generated CLI reference](reference/commands/dft.md) for every option.

> **Note:** In YAML (`--config`), the [`dft`](yaml-reference.md#dft-section) section holds the same settings. `dft.pyscf` passes attributes to PySCF objects by name, for example `pyscf: {mf: {level_shift: 0.2}}` for a hard-to-converge SCF. The charge and multiplicity of the ML region go in `calc.model_charge` and `calc.model_mult`; `-q`, `-l`, and `-m` come first, then YAML.

---

## Notes

* **Requirements**: `dft` needs the DFT extra: `pip install "mlmm-toolkit[dft]"` with the `cu130` or `cu132` PyTorch wheel, or `pip install "mlmm-toolkit[dft-cuda12]"` with `cu126`.
* **Basis cost**: `def2-tzvpd` costs much more than `def2-svp`. There is no fixed limit on atoms or GPU memory; the cost depends on the number of basis functions, the elements, the functional, the grid (`--dft-grid-level`), and the GPU. Run one representative structure first and watch the peak memory; if it runs out, use a smaller basis or a GPU with more memory.
* **GPU**: if GPU4PySCF cannot run, `dft` stops with an error that suggests the CPU engine; it does not switch to the CPU by itself. On a new GPU generation such as Blackwell (RTX 50xx), an out-of-memory or unsupported-kernel error can come from the GPU4PySCF and CuPy versions rather than from memory, so check the versions and the traceback first.
* **CPU**: `--dft-engine cpu` needs no GPU. How large an ML region is practical depends on the method and the machine, so time one representative single point.
* **Machines other than x86**: prebuilt GPU4PySCF wheels may not support them; build GPU4PySCF from source (https://github.com/pyscf/gpu4pyscf).
* **ECP**: the def2 ECP is attached for any basis whose name begins with `def2`, whatever the elements; the console prints `[dft] Using ECP: …`.
* **IAO analysis** can fail on difficult systems; its column in `result.yaml` is then `null`.
* **SCF not converged**: `dft` prints `WARNING: SCF did not converge.`, still writes `result.yaml` (and `result.json` with `--out-json`) with `converged: false`, and exits with code 1. In low-memory mode it suggests retrying with density fitting, `--no-dft-low-memory` (alias `--no-lowmem`), when memory allows.
* **Multiplicity** below 1 is rejected.
* **Earlier results**: a new run first removes `result.yaml`, `result.json`, `summary.json`, and the four `ml_region_*` files left in the output directory.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` and `--dft` in a workflow, DFT settings, and GPU memory
* [sp](sp.md) — single-point ML/MM energy and forces with any backend, including `-b dft`
* [all](all.md) — the full workflow; `--dft` adds DFT single points on R, TS, and P
* [MLIP Backends](backends.md) — choosing a backend
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
