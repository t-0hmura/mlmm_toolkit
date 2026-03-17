# `mm-parm`

## Overview

> **Summary:** Build Amber prmtop/rst7 topology files from a PDB with automatic GAFF2 ligand parameterization, disulfide detection, and optional hydrogen addition via PDBFixer.

`mlmm mm-parm` generates Amber topology/coordinate files from a PDB. The input PDB is used without structural fixing by default. Unknown residues are automatically parameterized with antechamber (GAFF2, AM1-BCC charges) and parmchk2. Residues explicitly listed in `--ligand-charge` are treated as ligand/cofactor definitions and are prioritized for GAFF2 parameterization. Disulfide bonds are inferred geometrically from SG-SG (or S-S) distances within 2.3 Å. Amino-acid residues listed in `AMINO_ACIDS` but still unrecognized by the selected force field are **not** handled automatically -- the build aborts with a message asking you to parameterize them manually.

When `--add-h`, hydrogens are added at the specified `--ph` using PDBFixer before tleap processing. No other structural fixing is performed. With `--ff-set ff14SB`, the force field switches to ff14SB (proteins) + TIP3P (water) (+ phosaa14SB).

## Minimal example

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

## Output checklist

- `<prefix>.parm7` -- Amber prmtop topology
- `<prefix>.rst7` -- Amber ASCII inpcrd coordinates
- `<prefix>.pdb` -- LEaP savepdb output (see naming rules below)

## Common examples

1. Basic topology build with hydrogen addition.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

2. Build without hydrogen addition.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3" --no-add-h
```

## Workflow

1. **Input preparation** -- The input PDB is read as-is (no structural fixing). If `--add-h` is set, hydrogens are added via PDBFixer at the specified `--ph`.
2. **TER insertion** -- When `--add-ter` (default), TER records are inserted before and after contiguous blocks of ligand/water/ion residues.
3. **Unknown residue parameterization** -- Residues not recognized by the force field are parameterized with antechamber (GAFF2, AM1-BCC) and parmchk2. Residues named in `--ligand-charge` are prioritized for this route. Formal charge and spin multiplicity are controlled via `--ligand-charge` and `--ligand-mult`.
4. **Disulfide detection** -- CYS/CYM/CYX pairs with SG-SG (or S-S) distance <= 2.3 Å are bonded automatically.
5. **Topology build** -- tleap generates parm7/rst7/pdb files using the selected force field set.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB (used as-is unless `--add-h`). | Required |
| `--out-prefix TEXT` | Output prefix for parm7/rst7/pdb files. | Stem of input PDB |
| `-l, --ligand-charge TEXT` | Map residue name to formal charge, e.g. `"GPP=-3,MMT=-1"`. | _None_ |
| `--ligand-mult TEXT` | Map residue name to spin multiplicity, e.g. `"HEM=1,NO=2"`. Unspecified residues default to singlet (1). | _None_ |
| `--keep-temp/--no-keep-temp` | Keep intermediate files/logs in a working directory (for debugging). | `False` |
| `--add-ter/--no-add-ter` | Insert TER before/after ligand/water/ion blocks. | `True` |
| `--add-h/--no-add-h` | Add hydrogens at `--ph` using PDBFixer. | `False` |
| `--ph FLOAT` | pH for PDBFixer hydrogen addition (used only with `--add-h`). | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | Force field set: ff19SB (default) or ff14SB. | `ff19SB` |

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [all](all.md) -- End-to-end workflow (calls mm-parm internally)
- [extract](extract.md) -- Extract active-site pocket before parameterization
- [define-layer](define_layer.md) -- Define ML/MM layers after building topology
