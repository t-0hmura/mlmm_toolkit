# `oniom-export` (Gaussian ONIOM / ORCA QM/MM input)

`oniom-export` **writes an mlmm ML/MM system as an input file for Gaussian ONIOM (`--mode g16`) or ORCA QM/MM (`--mode orca`)**. It reads the Amber topology (`--parm7`) and a PDB whose B-factors hold the layers. It uses the ML region as the QM region and writes the coordinates, the QM and movable atoms, and the MM parameters into one input file. The topology must be free of CMAP terms; [mm-parm](mm-parm.md#cmap-free-topology-for-oniom-export) shows how to build one.

---

## What it is for

* **Gaussian ONIOM**: take a structure from mlmm, such as a TS candidate, into a Gaussian ONIOM calculation with a DFT high layer
* **ORCA QM/MM**: run the same system with the QM/MM module of ORCA
* **Round trip**: edit the exported input outside mlmm and bring it back with its atom and residue names through [`oniom-import`](oniom-import.md) `--ref-pdb`

---

## Examples

### 1. Gaussian ONIOM (--mode g16)

Write the TS candidate from `mlmm tsopt` as a Gaussian ONIOM input. Here `result_tsopt/final_geometry.pdb` is the full system with the layers in its B-factors, `real.parm7` is the topology of the same system, and `ml_region.pdb` selects the QM atoms.

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1
g16 < ts_refine.com > ts_refine.log
```

The console prints `[oniom-gaussian] Wrote 'ts_refine.com'` followed by the numbers of QM atoms, movable atoms, and link boundaries.

When you already trust the atom order, `--no-element-check` skips the atom-by-atom element comparison with the topology.

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.gjf -q 0 -m 1 --no-element-check
```

### 2. ORCA QM/MM (--mode orca)

Write the same structure as an ORCA QM/MM input. The `.inp` suffix selects ORCA mode, so `--mode orca` can be left out.

```bash
mlmm oniom-export --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1
```

The console prints `[oniom-orca] Wrote 'ts_refine.inp'` and then `[oniom-orca] ORCAFF.prms: <path>` when the force-field file for ORCA is ready.

Set the charge and multiplicity of the whole QM+MM system (`Charge_Total`, `Mult_Total`) yourself:

```bash
mlmm oniom-export --mode orca --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1 --total-charge -1 --total-mult 1
```

Reuse an `ORCAFF.prms` you already have and skip the conversion step:

```bash
mlmm oniom-export --mode orca --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.inp -q 0 -m 1 \
    --orcaff ./ORCAFF.prms --no-convert-orcaff
```

### 3. Method, processors, and memory

Change the QM method and the resources written into the Gaussian input (`%nprocshared`, `%mem`).

```bash
mlmm oniom-export --mode g16 --parm7 real.parm7 -i result_tsopt/final_geometry.pdb \
    --model-pdb ml_region.pdb -o ts_refine.com -q 0 -m 1 \
    --method 'wb97xd/def2-svp' --nproc 16 --mem 32GB
```

---

## How it works

1. **Topology and layers**:
`oniom-export` reads the atoms, bonds, charges, and Amber parameters from the parm7, and the coordinates and B-factors from the PDB at `-i`. The PDB must list the same atoms in the parm7 order; `--element-check` compares their elements one by one. B-factors of 0, 10, and 20 (within ±1.0) mark the ML, movable MM, and frozen MM atoms.
2. **QM region**:
With `--model-pdb`, its atoms form the QM region; they are matched to `-i` by atom name, residue name, chain, residue number, and insertion code. Without it, the atoms with B-factor 0 form the QM region. Every atom except the frozen MM atoms is movable, and the QM atoms are always movable.
3. **QM/MM boundary**:
For Gaussian, a link H replaces the MM atom of each cut QM–MM bond: `--link-atom-method scaled` (the default) places it with the Morokuma/Dapprich g-factor, and `fixed` places it 1.09 Å (QM carbon) or 1.01 Å (QM nitrogen) from the QM atom. ORCA builds the caps itself from `QMAtoms` and `ORCAFF.prms`, and the input lists the estimated cap positions only as comments.
4. **Writing the input**:
The Gaussian input has the route `#p oniom(<method>:amber=softonly)`, the coordinates with the movable flag (`0` movable, `-1` frozen) and the layer (`H` or `L`), the connectivity, and the Amber parameters. Its charge and multiplicity line has three pairs: the whole system (the topology total charge and `-m`), then the QM region twice (`-q` and `-m`). The ORCA input has `! <method>` and `! QMMM`, a `%qmmm` block with `ORCAFFFilename`, `QMAtoms`, `ActiveAtoms`, `Charge_Total`, and `Mult_Total`, and the coordinates under `* xyz` with the QM charge and multiplicity (`-q`, `-m`). In ORCA mode, `oniom-export` also finds or creates `ORCAFF.prms`. The `%qmmm` keywords are described in the [ORCA 6.0 manual (QM/MM)](https://www.faccts.de/docs/orca/6.0/manual/contents/typical/qmmm.html).

---

## Output files

* **Gaussian input** (`--mode g16`): the `-o` file (`.com` or `.gjf`). The console prints `[oniom-gaussian] Wrote '<file>'`, `QM atoms: N, Movable atoms: M`, and `Link boundaries: K`.
* **ORCA input** (`--mode orca`): the `-o` file (`.inp`). The console prints `[oniom-orca] Wrote '<file>'`, `QM atoms: N, Active atoms: M`, and `Link boundaries (auto-capped by ORCA): K`.
* **`ORCAFF.prms`** (ORCA): the `--orcaff` file, or `<parm7 stem>.ORCAFF.prms` in the directory of `-o`. An existing file is reused; a missing one is created with `orca_mm -convff -AMBER <parm7>` when `--convert-orcaff` is on and `orca_mm` is on `PATH`. If the file still does not exist, the console prints `[oniom-orca] NOTE: ORCAFF.prms not found at '<path>'. Run manually: cd <dir> && orca_mm -convff -AMBER <parm7>`, and the `.inp` is complete only after you run that command.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--parm7` | path | (required) | Amber topology of the full system, free of CMAP terms |
| `-i, --input` | path | (required) | PDB of the full system in the parm7 atom order, with the layers in the B-factors (0, 10, 20) |
| `--model-pdb` | path | `None` | PDB of the QM atoms; without it, the atoms with B-factor 0 in `-i` |
| `-o, --output` | path | (required) | Input file to write: `.com` / `.gjf` (g16) or `.inp` (ORCA) |
| `--mode` | `g16` / `orca` | from the `-o` suffix | Program to write the input for |
| `--method` | text | `wB97XD/def2-TZVPD` (g16), `B3LYP D3BJ def2-SVP` (ORCA) | QM method and basis set, written into `oniom(<method>:amber=softonly)` (g16) or the `!` line (ORCA) |
| `-q, --charge` | integer | (required) | Charge of the QM region |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity of the QM region |
| `--nproc` | integer | `8` | Number of processors (`%nprocshared` for g16, `%pal nprocs` for ORCA) |
| `--mem` | text | `16GB` | g16: memory (`%mem`) |
| `--total-charge`, `--total-mult` | integer | topology total charge, `-m` | ORCA: charge and multiplicity of the whole QM+MM system (`Charge_Total`, `Mult_Total`) |
| `--orcaff` | path | `<parm7 stem>.ORCAFF.prms` in the directory of `-o` | ORCA: an existing `ORCAFF.prms` to use |
| `--convert-orcaff/--no-convert-orcaff` | flag | `True` | ORCA: when `--orcaff` is not given and the default file is missing, create it with `orca_mm -convff -AMBER` |
| `--element-check/--no-element-check` | flag | `True` | Compare the elements of `-i` with the topology atom by atom |
| `--link-atom-method` | `scaled` / `fixed` | `scaled` | g16: place link H atoms with the g-factor (`scaled`) or at a fixed bond length (`fixed`) |

See the [generated CLI reference](reference/commands/oniom_export.md) for every option.

---

## Notes

* **CMAP**: Gaussian ONIOM cannot represent the CMAP terms of a parm7, and the MM engine of ORCA does not apply them, so a topology with CMAP stops the export before any file is written. Build a CMAP-free topology for the export; the ML/MM calculations in mlmm can keep CMAP, which they apply in both MM layers.
* **Choosing the mode**: `--mode` takes precedence over the suffix; without it, a suffix other than `.gjf`, `.com`, or `.inp` is an error.
* **Atom order**: `-i` must be a PDB (`.pdb` or `.ent`) with the same number of atoms as the parm7. A different atom count stops the export even with `--no-element-check`; with the check on, the first different element stops it with `Element sequence mismatch at atom index …` (counted from 0).
* **Whole-system charge**: the charge of the Gaussian real system, and the default ORCA `Charge_Total`, is the sum of the parm7 partial charges rounded to an integer. If the sum is more than 0.05 from an integer, the export stops; in ORCA mode, give `--total-charge` instead.
* **Gaussian boundaries**: each cut QM–MM bond needs its own MM atom. If two QM atoms are bonded to the same MM atom, the Gaussian export stops.
* **Job type**: the exported input has no job keyword (such as `opt` or `freq`) on the Gaussian route line or the ORCA `!` line, so as written it is a single-point calculation. Add the keywords for the job you want before running it.
* **`ORCAFF.prms` before running ORCA**: the `.inp` refers to `ORCAFF.prms` by its absolute path. Check that this file exists before you run the `.inp`, also after moving the input to another machine.
* **Atom-order marker**: the exported file carries `MLMM_REF_PDB_ORDER_V1_SHA256=<digest>`, a hash of the names, numbers, and other identity fields of every atom in `-i`; coordinates, occupancy, and B-factors are left out. [`oniom-import`](oniom-import.md) `--ref-pdb` checks this marker before it copies the names back.
* **Multiplicity**: a value of `-m` below 1 is rejected on the command line.
* **Requirements**: Gaussian and ORCA are not part of mlmm-toolkit; install and license them separately.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [oniom-import](oniom-import.md) — bring an edited ONIOM input back as XYZ and a layered PDB
* [mm-parm](mm-parm.md) — build the Amber topology, including a CMAP-free one
* [define-layer](define-layer.md) — write the layer B-factors into the full-system PDB
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
