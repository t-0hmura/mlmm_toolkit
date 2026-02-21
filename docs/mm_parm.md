# `mm-parm`

## Overview

> **Summary:** Build Amber prmtop/rst7 topology files from a PDB with automatic GAFF2 ligand parameterization, disulfide detection, and optional hydrogen addition via PDBFixer.

### Quick reference
- **Use when:** You need AMBER topology and coordinate files for an ML/MM calculation.
- **Input:** One PDB file (used as-is unless `--add-h`).
- **Defaults:** ff19SB (proteins) + OPC3 (water) + GAFF2 (organics) + Lipid21 + GLYCAM_06j-1 + OL3/OL21.
- **Outputs:** `<prefix>.parm7`, `<prefix>.rst7`, and optionally `<prefix>.pdb` (LEaP output).
- **Requirements:** AmberTools (tleap, antechamber, parmchk2) on PATH; PDBFixer + OpenMM only when `--add-h`.

`mlmm mm-parm` generates Amber topology/coordinate files from a PDB. The input PDB is used without structural fixing by default. Unknown residues are automatically parameterized with antechamber (GAFF2, AM1-BCC charges) and parmchk2. Disulfide bonds are inferred geometrically from SG-SG (or S-S) distances within 2.3 A. Nonstandard amino acids (residues containing N/CA/C that are not recognized by the selected force field) are **not** handled automatically -- the build aborts with a message asking you to parameterize them manually.

When `--add-h`, hydrogens are added at the specified `--pH` using PDBFixer before tleap processing. No other structural fixing is performed. With `--ff-set ff14SB`, the force field switches to ff14SB (proteins) + TIP3P (water) (+ phosaa14SB).

## Usage
```bash
mlmm mm-parm -i INPUT.pdb [--out-prefix PREFIX] \
             [--ligand-charge "RES=-3,ABC:+1"] \
             [--ligand-mult "RES=2,ABC:1"] \
             [--keep-temp] \
             [--add-ter/--no-add-ter] \
             [--add-h/--no-add-h] [--pH 7.0] \
             [--ff-set ff19SB|ff14SB]
```

### Example
```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
    --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
    --add-ter --ff-set ff19SB --add-h --pH 7.0
```

## Description

### Disulfide detection
CYS/CYM/CYX disulfide bonds are inferred purely geometrically from PDB coordinates: residue pairs with SG-SG (or S-S) distance <= 2.3 A are bonded.

### Force fields
- **ff19SB** (default): ff19SB (proteins) + OPC3 (water) + GAFF2 (general organics) + Lipid21 + GLYCAM_06j-1 + OL3/OL21 (+ phosaa19SB).
- **ff14SB**: ff14SB (proteins) + TIP3P (water) (+ phosaa14SB).

### Unknown residues
If LEaP reports unknown residues, they are parameterized automatically with antechamber (GAFF2, AM1-BCC charges) and parmchk2. Formal charge and spin multiplicity can be controlled with `--ligand-charge` and `--ligand-mult`.

### Automatic TER insertion
When `--add-ter` (default), TER records are inserted **before and after** contiguous blocks of residues whose names appear in `--ligand-charge`, water residues, or ions. Consecutive such residues do not receive TER records between them.

### Failure behavior
Even if the build fails, when `--add-h` and hydrogen addition succeeded, the hydrogen-added PDB is still exported to disk using the same naming rule as the LEaP PDB (`<out_prefix>.pdb` or `<input_stem>_parm.pdb` when `--out-prefix` is omitted).

### Output naming rules
- `<prefix>.parm7` -- prmtop topology
- `<prefix>.rst7` -- ASCII inpcrd (the LEaP-produced `complex.inpcrd` is copied)
- `<prefix>.pdb` -- LEaP `savepdb` output:
  - If `--out-prefix` is given: `<out_prefix>.pdb`
  - If `--out-prefix` is omitted and `--add-h`: `<input_stem>_parm.pdb`
  - If `--out-prefix` is omitted and `--no-add-h`: no PDB is written

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB (used as-is unless `--add-h`). | Required |
| `--out-prefix TEXT` | Output prefix for parm7/rst7/pdb files. | Stem of input PDB |
| `--ligand-charge TEXT` | Map residue name to formal charge, e.g. `"GPP=-3,MMT=-1"`. | _None_ |
| `--ligand-mult TEXT` | Map residue name to spin multiplicity (`-m`), e.g. `"HEM=1,NO:2"`. | `1` |
| `--keep-temp` | Keep intermediate files/logs in a working directory (for debugging). | `False` |
| `--add-ter/--no-add-ter` | Insert TER before/after ligand/water/ion blocks. | `True` |
| `--add-h/--no-add-h` | Add hydrogens at `--pH` using PDBFixer. | `False` |
| `--pH FLOAT` | pH for PDBFixer hydrogen addition (used only with `--add-h`). | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | Force field set: ff19SB (default) or ff14SB. | `ff19SB` |

## Outputs
```
<prefix>.parm7         # Amber prmtop topology
<prefix>.rst7          # Amber ASCII inpcrd coordinates
<prefix>.pdb           # LEaP savepdb output (see naming rules above)
```

## Requirements
- **AmberTools** (`tleap`, `antechamber`, `parmchk2`) must be available on PATH.
- **PDBFixer** + **OpenMM** are required **only** when using `--add-h`.

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [all](all.md) -- End-to-end workflow (calls mm-parm internally)
- [extract](extract.md) -- Extract active-site pocket before parameterization
- [define-layer](define_layer.md) -- Define ML/MM layers after building topology
