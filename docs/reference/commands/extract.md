# `mlmm extract`

```text
Usage: mlmm extract [OPTIONS]

  Extract an active site model around substrate residues (from PDB/mmCIF or
  residue IDs/names), with biochemically aware truncation and optional cap-H;
  mmCIF inputs also produce mmCIF outputs.

Options:
  -v, --verbose LEVEL             Console verbosity 0-3 (default 2). 0=silent;
                                  1=milestones only; 2=+detailed step logging
                                  and deliverable paths; 3=everything (full
                                  config blocks, per-file paths, DEBUG logging).
                                  [0<=x<=3]
  --help-advanced                 Show all options (including advanced settings)
                                  and exit.
  -i, --input TEXT                Protein-substrate complex PDB/mmCIF file(s).
                                  Multiple files may be given space-separated
                                  after one -i or by repeating -i. PDBs beyond
                                  fixed-column residue/atom limits are handled
                                  through an internal safe bridge. If multiple,
                                  they must have identical atom counts and
                                  ordering.  [required]
  -c, --center TEXT               Substrate specification: a PDB/mmCIF path, a
                                  comma/space-separated residue-ID list like
                                  '123,124' or 'A:123,B:456' (insertion codes
                                  supported), a residue-name list like
                                  'GPP,SAM', or a chain-qualified name like
                                  'A:SAM' (all matches in chain A) / 'A:SAM:123'
                                  (one residue).  [required]
  -o, --output TEXT               Internal/output PDB path(s). For mmCIF or
                                  oversized-PDB input, a .cif companion with the
                                  original chain/residue IDs is written
                                  automatically. One path creates multi-MODEL
                                  output; N paths create one output per input.
  -r, --radius FLOAT              Cutoff (angstrom) around substrate atoms for
                                  active site model inclusion.  [default: 2.6]
  --radius-het2het FLOAT          Cutoff (angstrom) for substrate hetero-atom
                                  (non-C/H) to neighbor hetero-atom proximity. 0
                                  is treated as 0.001 angstrom (effectively
                                  off).  [default: 0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/TIP3/SOL).  [default:
                                  include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Delete main-chain atoms from non-substrate
                                  amino acids.  [default: no-exclude-backbone]
  --add-linkh / --no-add-linkh    Add cap hydrogens (carbon boundaries only) at
                                  1.09 angstrom along cut-bond directions.
                                  [default: no-add-linkh]
  --selected-resn TEXT            Comma/space-separated residue IDs/names to
                                  force-include; chain-qualified A:SAM is
                                  supported.
  --modified-residue TEXT         Comma-separated residue names (with optional
                                  charge) to treat as amino acids for backbone
                                  truncation and charge assignment. Examples:
                                  'HD1,HD2,HD3' (charge defaults to 0) or
                                  'HD1:0,SEP:-2'.
  -l, --ligand-charge TEXT        Total charge number or per-resname mapping
                                  like 'GPP:-3,SAM:1'.
  --out-json / --no-out-json      Write machine-readable result.json next to the
                                  output PDB.  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
