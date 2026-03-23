# mlmm mm-parm

```
mlmm-toolkit ver. 0.2.5.dev18

Usage: cli mm-parm [OPTIONS]

  Generate Amber parm7/rst7 (and a LEaP-exported PDB) from a PDB using
  AmberTools only.

Options:
  --help-advanced           Show all options (including advanced settings) and
                            exit.
  -i, --input FILE          Input PDB file (used as-is; optional hydrogens via
                            --add-h/--ph).  [required]
  --out-prefix TEXT         Output prefix (default: input PDB stem). For LEaP
                            PDB: if omitted and --add-h True,
                            <input_stem>_parm.pdb is used.
  -l, --ligand-charge TEXT  Comma-separated mapping of residue=charge or
                            residue:charge (e.g., "GPP=-3,MMT=-1" or
                            "GPP:-3,MMT:-1")
  --ligand-mult TEXT        Comma-separated mapping of residue=multiplicity or
                            residue:multiplicity (e.g., "HEM=1,NO:2")
  --ff-set [ff19SB|ff14SB]  Force-field set for proteins/backbone typing and
                            water/ion parameters (default: ff19SB).
  -h, --help                Show this message and exit.
```
