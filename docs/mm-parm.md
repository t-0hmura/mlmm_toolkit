# `mm-parm`

`mlmm mm-parm` generates Amber topology/coordinate files (parm7/rst7/pdb) from a PDB using AmberTools tleap. Unknown residues are auto-parameterized with GAFF2 (AM1-BCC charges); see Workflow for the full pipeline, and CLI options for the force-field and hydrogen-addition flags.

## Examples

Basic build (ligand charges + multiplicities):

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1"
```

Add TER records, ff19SB, and hydrogens at pH 7:

```bash
mlmm mm-parm -i input.pdb --out-prefix complex \
 -l "GPP=-3,MMT=-1" --ligand-mult "GPP=1,MMT=1" \
 --add-ter --ff-set ff19SB --add-h --ph 7.0
```

Skip hydrogen addition (input already protonated):

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

## Outputs

- `<prefix>.parm7` -- Amber prmtop topology
- `<prefix>.rst7` -- Amber ASCII inpcrd coordinates
- `<prefix>.pdb` -- LEaP savepdb output

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB (used as-is unless `--add-h`). | Required |
| `-o, --out-prefix TEXT` | Output prefix for parm7/rst7/pdb files. | Stem of input PDB |
| `-l, --ligand-charge TEXT` | Map residue name to formal charge, e.g. `"GPP=-3,MMT=-1"`. | _None_ |
| `--ligand-mult TEXT` | Map residue name to spin multiplicity, e.g. `"HEM=1,NO=2"`. Unspecified residues default to singlet (1). | _None_ |
| `--keep-temp/--no-keep-temp` | Keep intermediate files/logs in a working directory (for debugging). | `False` |
| `--add-ter/--no-add-ter` | Insert TER before/after ligand/water/ion blocks. | `True` |
| `--add-h/--no-add-h` | Add hydrogens at `--ph` using PDBFixer. | `False` |
| `--ph FLOAT` | pH for PDBFixer hydrogen addition (used only with `--add-h`). | `7.0` |
| `--ff-set {ff19SB\|ff14SB}` | Force field set: ff19SB (default) or ff14SB. | `ff19SB` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## Notes

`mm-parm` relies on AmberTools tleap with GAFF2 automatic parameterization and works well when the substrate is a **typical organic molecule**. For the following cases, it is strongly recommended to prepare your own topology externally (e.g. with tleap, MCPB.py, or glycam.org tools) and supply it via the `--parm` flag of each subcommand:

- **Metalloenzymes** -- Metal centers require specialized bonded/non-bonded parameters (e.g. MCPB.py, the bonded model, or ZAFF). Automatic GAFF2 parameterization cannot handle metal-ligand coordination.
- **Glycans and carbohydrate-containing systems** -- Glycan linkages need GLYCAM force field parameters that are not included in the standard GAFF2/ff19SB setup.
- **Non-standard amino acids or post-translational modifications** -- Phosphorylated, methylated, or other modified residues may require custom `frcmod`/`lib` files.
- **MD snapshot initial structures** -- When starting from an MD trajectory snapshot, reusing the same `.parm7` file from the MD simulation is the most appropriate approach. This ensures consistency between the MM energy surface used for ML/MM and the one used in the preceding MD, avoiding artifacts from re-parameterization (e.g. different partial charges or atom-type assignments).
- Amino-acid residues listed in `AMINO_ACIDS` but still unrecognized by the selected force field are not handled automatically -- the build aborts with a message asking you to parameterize them manually.
- `--ff-set ff14SB` switches the force field to ff14SB (proteins) + TIP3P (water) (+ phosaa14SB); the default `ff19SB` set is used otherwise.
- **4-point water with a virtual site (OPC, TIP4P/-Ew, TIP5P) is _not yet_ supported by mlmm-toolkit's default MM backend** — the massless extra point (Amber `EPW`, element `EP`) is treated as a free atom, so the geometry optimizer stalls. Use a **3-point water model**: the default `ff19SB` set already builds **OPC3** (the recommended 3-point model), and `--ff-set ff14SB` builds TIP3P. To keep 4-point water, run with `--mm-backend openmm`, which places the virtual sites correctly but is slower (finite-difference MM Hessian).

```bash
# Example: supply a pre-built topology from MD
mlmm opt -i snapshot_layered.pdb --parm md_system.parm7 -q -1 -m 1 \
  --opt-mode grad --out-dir result
```

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide
- [all](all.md) — End-to-end workflow (calls mm-parm internally)
- [extract](extract.md) — Extract active-site pocket before parameterization
- [define-layer](define-layer.md) — Define ML/MM layers after building topology
