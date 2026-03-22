# `mlmm extract`

```text
mlmm-toolkit ver. 0.2.5.dev18

Usage: mlmm extract [OPTIONS]

  Extract a binding pocket around substrate residues (from a PDB or residue
  IDs/names), with biochemically aware truncation and optional link-H; supports
  multi-structure input and multi-MODEL output.

Options:
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input TEXT                Protein-substrate complex PDB(s). If multiple,
                                  they must have identical atom counts and
                                  ordering.  [required]
  -c, --center TEXT               Substrate specification: a PDB path, a
                                  comma/space-separated residue-ID list like
                                  '123,124' or 'A:123,B:456' (insertion codes
                                  supported), or a residue-name list like
                                  'GPP,SAM'.  [required]
  -o, --output TEXT               Output PDB path(s). One path for multi-MODEL
                                  PDB, or N paths for per-file output. If
                                  omitted: single input -> pocket.pdb; multiple
                                  inputs -> pocket_{filename}.pdb.
  -r, --radius FLOAT              Cutoff (angstrom) around substrate atoms for
                                  pocket inclusion.  [default: 2.6]
  --radius-het2het FLOAT          Cutoff (angstrom) for substrate-protein
                                  hetero-atom proximity (non-C/H). 0 disables.
                                  [default: 0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL).
                                  [default: include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Delete main-chain atoms from non-substrate
                                  amino acids.  [default: no-exclude-backbone]
  --add-linkh / --no-add-linkh    Add carbon-only link-H at 1.09 angstrom along
                                  cut-bond directions.  [default: no-add-linkh]
  --selected-resn TEXT            Comma/space-separated residue IDs to force-
                                  include.
  -l, --ligand-charge TEXT        Total charge number or per-resname mapping
                                  like 'GPP:-3,SAM:1'.
  -v, --verbose / --no-verbose    Enable INFO-level logging.  [default: v]
  -h, --help                      Show this message and exit.
```
