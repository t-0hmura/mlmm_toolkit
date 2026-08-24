# `mlmm extract`

```text
Usage: mlmm extract [OPTIONS]

  Extract an active site model around substrate residues (from PDB/mmCIF or
  residue IDs/names), with biochemically aware truncation and optional link-H;
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
  -r, --radius FLOAT RANGE        Cutoff (Å) around substrate atoms for active-
                                  site inclusion. Zero is accepted and evaluated
                                  internally as 0.001 Å (effectively off for
                                  ordinary radius-based neighbors).  [default:
                                  2.6; x>=0.0]
  --radius-het2het FLOAT RANGE    Cutoff (Å) for substrate hetero-atom (non-C/H)
                                  to neighbor hetero-atom proximity. 0 is
                                  treated as 0.001 Å (effectively off).
                                  [default: 0; x>=0.0]
  --include-h2o / --no-include-h2o
                                  Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL).
                                  [default: include-h2o]
  --exclude-backbone / --no-exclude-backbone
                                  Delete main-chain atoms from non-substrate
                                  amino acids.  [default: no-exclude-backbone]
  --add-linkh / --no-add-linkh    Add link hydrogens (carbon boundaries only) at
                                  1.09 Å along cut-bond directions.  [default:
                                  no-add-linkh]
  --selected-resn TEXT            Force-include residues using IDs ('123',
                                  'A:123A'), names ('SAM'), or chain-qualified
                                  names ('A:SAM', 'A:SAM:123'); comma/space
                                  separated.  [default: ""]
  --modified-residue TEXT         Comma-separated residue names with charges to
                                  treat as amino acids for backbone truncation
                                  and charge assignment. A known catalog residue
                                  may omit its charge. Example: 'HD1:0,SEP'.
  -l, --ligand-charge TEXT        Total charge number or per-resname mapping
                                  like 'GPP:-3,SAM:1'.
  --out-json / --no-out-json      Write machine-readable result.json next to the
                                  output PDB.  [default: no-out-json]
  -h, --help                      Show this message and exit.
```
