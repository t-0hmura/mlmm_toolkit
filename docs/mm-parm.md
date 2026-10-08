# `mm-parm` (build the Amber topology)

`mm-parm` builds an Amber topology (`parm7`), coordinates (`rst7`), and a matching PDB from a PDB of the whole enzyme–substrate complex with AmberTools tleap. Residues that the force field does not know, such as a substrate or a cofactor, get GAFF2 parameters with AM1-BCC charges. Every ML/MM command reads the resulting `parm7` through `--parm7`.

---

## What it is for

* **Building the MM topology**: write `parm7`, `rst7`, and a PDB for the whole system.
* **Parameterizing ligands**: give unknown residues GAFF2 parameters with the formal charge and multiplicity from `-l` and `--ligand-mult`.
* **Preparing a model by hand**: pass the PDB that `mm-parm` writes to `extract` and `define-layer`; its atoms are in the same order as the `parm7`.

---

## Examples

### 1. Build with ligand charges and multiplicities

Give each ligand's formal charge and spin multiplicity by residue name.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3,MMT:-1' --ligand-mult 'GPP:1,MMT:1'
```

The console prints `[mm-parm] Wrote:` lines for `complex.pdb`, `complex.parm7`, and `complex.rst7`.

### 2. Add hydrogens at pH 7.0

Let PDBFixer add hydrogens before the build.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3,MMT:-1' --ligand-mult 'GPP:1,MMT:1' \
    --add-ter --ff-set ff19SB --add-h --ph 7.0
```

### 3. Input that already has hydrogens

Leave the input as it is.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    -l 'GPP:-3' --no-add-h
```

### 4. Prepare the model by hand

Build the topology, cut the ML region out of the PDB that `mm-parm` writes, and assign the layers on the same PDB.

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

The calculation commands then take `system_layered.pdb` with `--parm7 system.parm7`.

### 5. Topology for mmCIF inputs

`mm-parm` reads PDB only. For mmCIF inputs, write the full system to a PDB with the same atom order and elements, build the topology from it, and pass the `parm7` to `all` with the mmCIF structures.

```bash
mlmm mm-parm -i reactant_topology.pdb -l 'SAM:1,GPP:-3' \
    --out-prefix full_system
mlmm all -i reactant.cif product.cif --parm7 full_system.parm7 \
    -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result
```

---

## How it works

1. **Input**: the PDB is used as it is. With `--add-h`, PDBFixer adds hydrogens at `--ph`; it adds no missing heavy atoms or residues.
2. **TER records**: with `--add-ter` (default), a `TER` record is inserted before and after each block of residues named in `-l`, waters, and ions, without splitting a block of consecutive such residues. A `TER` is also inserted between neighboring amino acids that are not joined by a peptide C–N bond (different chains, or C–N > 1.9 Å).
3. **Disulfides**: CYS/CYX pairs whose SG atoms are within 2.5 Å are bonded, and a bonded CYS is renamed CYX so that tleap removes its HG. With `--no-auto-disulfide`, only residues already named CYX are bonded.
4. **Unknown residues**: tleap first runs with the force field alone. Each residue name that tleap reports as unknown is parameterized with antechamber (GAFF2, AM1-BCC) and parmchk2 from the first residue of that name in the file, with the charge from `-l` (0 if not given) and the multiplicity from `--ligand-mult` (1 if not given). Before antechamber runs, the electron count of the residue is checked against that charge and multiplicity.
5. **Topology**: tleap runs again with the new parameters and writes the topology, the coordinates, and the PDB. `mm-parm` fills blank element columns of the PDB from the `parm7` and keeps every record and the atom order unchanged.

---

## Output files

```text
./
├─ <prefix>.parm7   # Amber topology
├─ <prefix>.rst7    # Amber coordinates (ASCII)
└─ <prefix>.pdb     # tleap's PDB with element columns filled; same atoms and order as the parm7
```

`<prefix>` defaults to the input file name without its extension, in the current directory. Choose a prefix different from the input name, so that `<prefix>.pdb` does not replace the input. The PDB is written when `--out-prefix` is given, and as `<input name>_parm.pdb` when `--add-h` is given without `--out-prefix`; otherwise only `parm7` and `rst7` are written. If the build fails after `--add-h`, the hydrogen-added structure is written to that PDB path, unless a file already exists there. `--keep-temp` keeps the working directory `parm7build_*`, with the tleap logs, in the current directory.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input PDB, used as it is unless `--add-h` is given |
| `-o, --out-prefix` | text | input name without extension | Prefix of the output files |
| `-l, --ligand-charge` | text | `None` | Formal charge per residue name (e.g. `'GPP:-3,MMT:-1'`) |
| `--ligand-mult` | text | `1` | Spin multiplicity per residue name (e.g. `'HEM:1,NO:2'`) |
| `--keep-temp/--no-keep-temp` | flag | `False` | Keep the working directory and the tleap logs |
| `--add-ter/--no-add-ter` | flag | `True` | Insert `TER` before and after ligand, water, and ion blocks, and between amino acids that are not joined by a peptide bond |
| `--auto-disulfide/--no-auto-disulfide` | flag | `True` | Bond CYS/CYX pairs with SG–SG ≤ 2.5 Å and rename a bonded CYS to CYX. Off: bond only residues already named CYX |
| `--add-h/--no-add-h` | flag | `False` | Add hydrogens with PDBFixer at `--ph` |
| `--ph` | float | `7.0` | pH for `--add-h` |
| `--ff-set` | `ff19SB` or `ff14SB` | `ff19SB` | Force-field set (see [Notes](#notes)) |

See the [generated CLI reference](reference/commands/mm_parm.md) for every option.

---

## CMAP-free topology for `oniom-export`

Use ff14SB and confirm that the resulting topology has no CMAP terms:

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --ff-set ff14SB --out-prefix system
python -c "import parmed as pmd; p=pmd.load_file('system.parm7'); assert not p.cmaps"
```

---

## Notes

* **Systems that need a topology built elsewhere**: `mm-parm` works best when the substrate is a typical organic molecule. For the following systems, build the topology yourself and pass it with `--parm7`.
  * **Metalloenzymes**: metal centers need dedicated bonded and non-bonded parameters (MCPB.py, the bonded model, or ZAFF); GAFF2 cannot describe metal–ligand coordination.
  * **Glycans**: the force-field set loads GLYCAM_06j-1, but `mm-parm` creates no bonds other than disulfides, so glycosidic and other covalent links between residues need your own tleap `bond` commands.
  * **Non-standard amino acids and post-translational modifications**: modified residues may need their own `frcmod`/`lib` files.
  * **Structures from MD**: reuse the `parm7` of the MD run, so that the ML/MM calculation uses the same MM energy surface as the MD and does not change partial charges or atom types.

  ```bash
  # Use a topology built for MD
  mlmm opt -i snapshot_layered.pdb --parm7 md_system.parm7 -q -1 -m 1 \
    --opt-mode grad --out-dir result
  ```
* **The `parm7` follows atom order**: every coordinate input must contain the atoms of the whole system in the same order as the `parm7`. Before a calculation, the ML/MM calculator compares the atom count and, atom by atom, the element, the atom name (`1HB` matches `HB1`), the residue name, and the residue order, and stops at the first mismatch. Keep atom names, residue names, and the atom order when you prepare reactant, intermediate, and product structures. A `--model-pdb` file is only an unchanged subset of the full system that selects the ML atoms.
* **Residue numbers in the output PDB**: tleap numbers residues 1, 2, … in the order they appear, so the residue numbers of the PDB that `mm-parm` writes can differ from the input. In `examples/beza/1.R.pdb`, ARG 38 is the first residue and SAM 320 the 283rd. When you cut that PDB by hand, select residues by name or read their numbers from the file. `all` cuts the original input, so its selectors use the original numbers.
* **Amino acids that the force field does not know**: a residue that `extract` lists as an amino acid (see its appendix) but tleap does not know stops the build. The message offers three ways out: list the residue in `-l` to give it GAFF2 parameters, change the residue in the input, or build the topology yourself with tleap.
* **Ligand charge and hydrogens**: the electron-count check stops the build when the hydrogens of a residue do not match its charge and multiplicity (SAM: 22 H for charge 0, 23 H for +1). A `-l` or `--ligand-mult` entry for a residue that tleap already knows is not used, and a warning says so.
* **Force-field sets**: `ff19SB` loads ff19SB with phosaa19SB and ff19SB_modAA, OPC3 water, and its ion parameters. `ff14SB` loads ff14SB with phosaa14SB and ff14SB_modAA, TIP3P water, and its ion parameters. Both also load lipid21, RNA.OL3, DNA.OL21, GLYCAM_06j-1, and GAFF2.
* **Water with virtual sites**: the default MM backend, `hessian_ff`, rejects topologies whose water has massless virtual sites (OPC, TIP4P/-Ew, TIP5P) and reports their number and atom numbers. Use 3-point water, or run the calculation with `--mm-backend openmm`.
* **Requirements**: tleap, antechamber, and parmchk2 from AmberTools must be on `PATH`, and `--add-h` also needs PDBFixer; see [Installation](installation.md).

---

## See also

* [Building the ML region and layers](model-setup.md) — choose the ML region and the MM layers on the PDB that `mm-parm` writes
* [all](all.md) — the full workflow; runs `mm-parm` when `--parm7` is not given
* [extract](extract.md) — cut the ML region out of the topology-matched PDB
* [define-layer](define-layer.md) — assign the ML, Movable-MM, and Frozen-MM layers on the topology-matched PDB
* [oniom-export](oniom-export.md) — write Gaussian ONIOM or ORCA QM/MM input; needs the CMAP-free topology above
* [Troubleshooting](troubleshooting.md) — topology and atom-order errors
