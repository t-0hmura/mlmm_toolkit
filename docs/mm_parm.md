# `mm-parm`

## Overview

> **Summary:** Build Amber prmtop/rst7 topology files from a PDB with automatic GAFF2 ligand parameterization, disulfide detection, and optional hydrogen addition via PDBFixer.

### At a glance
- **Use when:** You need AMBER topology and coordinate files for an ML/MM calculation.
- **Method:** tleap + antechamber (GAFF2, AM1-BCC) + parmchk2; optional PDBFixer hydrogen addition.
- **Outputs:** `<prefix>.parm7`, `<prefix>.rst7`, and optionally `<prefix>.pdb` (LEaP output).
- **Defaults:** ff19SB (proteins) + OPC3 (water) + GAFF2 (organics) + Lipid21 + GLYCAM_06j-1 + OL3/OL21.
- **Next step:** [define-layer](define_layer.md) to assign ML/MM layers.

`mlmm mm-parm` generates Amber topology/coordinate files from a PDB. The input PDB is used without structural fixing by default. Unknown residues are automatically parameterized with antechamber (GAFF2, AM1-BCC charges) and parmchk2. Disulfide bonds are inferred geometrically from SG-SG (or S-S) distances within 2.3 A. Nonstandard amino acids (residues containing N/CA/C that are not recognized by the selected force field) are **not** handled automatically -- the build aborts with a message asking you to parameterize them manually.

When `--add-h`, hydrogens are added at the specified `--ph` using PDBFixer before tleap processing. No other structural fixing is performed. With `--ff-set ff14SB`, the force field switches to ff14SB (proteins) + TIP3P (water) (+ phosaa14SB).

## Minimal example

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

## Output checklist

- `<prefix>.parm7` -- Amber prmtop topology
- `<prefix>.rst7` -- Amber ASCII inpcrd coordinates
- `<prefix>.pdb` -- LEaP savepdb output (see naming rules below)

## Common examples

1. Basic topology build with hydrogen addition.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

2. Build without hydrogen addition.

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 --ligand-charge "GPP=-3" --no-add-h
```

## Workflow

1. **Input preparation** -- The input PDB is read as-is (no structural fixing). If `--add-h` is set, hydrogens are added via PDBFixer at the specified `--ph`.
2. **TER insertion** -- When `--add-ter` (default), TER records are inserted before and after contiguous blocks of ligand/water/ion residues.
3. **Unknown residue parameterization** -- Residues not recognized by the force field are parameterized with antechamber (GAFF2, AM1-BCC) and parmchk2. Formal charge and spin multiplicity are controlled via `--ligand-charge` and `--ligand-mult`.
4. **Disulfide detection** -- CYS/CYM/CYX pairs with SG-SG (or S-S) distance <= 2.3 A are bonded automatically.
5. **Topology build** -- tleap generates parm7/rst7/pdb files using the selected force field set.

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB (used as-is unless `--add-h`). | Required |
| `--out-prefix TEXT` | Output prefix for parm7/rst7/pdb files. | Stem of input PDB |
| `--ligand-charge TEXT` | Map residue name to formal charge, e.g. `"GPP=-3,MMT=-1"`. | _None_ |
| `--ligand-mult TEXT` | Map residue name to spin multiplicity (`-m`), e.g. `"HEM=1,NO:2"`. | `1` |
| `--allow-nonstandard-aa` | Allow antechamber parameterization for amino-acid-like modified residues (N/CA/C present). | `False` |
| `--keep-temp` | Keep intermediate files/logs in a working directory (for debugging). | `False` |
| `--add-ter/--no-add-ter` | Insert TER before/after ligand/water/ion blocks. | `True` |
| `--add-h/--no-add-h` | Add hydrogens at `--ph` using PDBFixer. | `False` |
| `--ph FLOAT` | pH for PDBFixer hydrogen addition (used only with `--add-h`). | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | Force field set: ff19SB (default) or ff14SB. | `ff19SB` |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- **AmberTools** (`tleap`, `antechamber`, `parmchk2`) must be available on PATH.
- **PDBFixer** + **OpenMM** are required **only** when using `--add-h`.
- Even if the build fails, when `--add-h` and hydrogen addition succeeded, the hydrogen-added PDB is still exported to disk.

### Force fields
- **ff19SB** (default): ff19SB (proteins) + OPC3 (water) + GAFF2 (general organics) + Lipid21 + GLYCAM_06j-1 + OL3/OL21 (+ phosaa19SB).
- **ff14SB**: ff14SB (proteins) + TIP3P (water) (+ phosaa14SB).

### Output naming rules
- `<prefix>.parm7` -- prmtop topology
- `<prefix>.rst7` -- ASCII inpcrd (the LEaP-produced `complex.inpcrd` is copied)
- `<prefix>.pdb` -- LEaP `savepdb` output:
 - If `--out-prefix` is given: `<out_prefix>.pdb`
 - If `--out-prefix` is omitted and `--add-h`: `<input_stem>_parm.pdb`
 - If `--out-prefix` is omitted and `--no-add-h`: no PDB is written

### Quick failure recovery
1. AmberTools commands are missing.

```bash
which tleap antechamber parmchk2
mlmm mm-parm -i input.pdb --keep-temp
```

2. Build stops on nonstandard residues.

```bash
mlmm mm-parm -i input.pdb --keep-temp
# Check tleap log in the kept temp directory, then provide custom residue params.
```

3. `--add-h` fails in the active environment.

```bash
python -c "import pdbfixer, openmm; print('ok')"
mlmm mm-parm -i input.pdb --no-add-h
```

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [all](all.md) -- End-to-end workflow (calls mm-parm internally)
- [extract](extract.md) -- Extract active-site pocket before parameterization
- [define-layer](define_layer.md) -- Define ML/MM layers after building topology
