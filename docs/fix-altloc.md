# `fix-altloc`

## Overview

> **Summary:** Remove alternate location (altLoc) indicators from PDB files by selecting the best conformer for each atom based on occupancy and dropping duplicates.

### What it does
1. Blank the PDB altLoc column (column 17, 1-based) with a single space.
 - This is a 1-character replacement (no shifting / no reformatting).
2. If the same atom appears multiple times due to alternate locations
    (altLoc like A/B/... or custom labels like H/L), keep the "best" one:
 - Highest occupancy first
 - If tied (or occupancy missing), keep the earliest one in the file

### Handled records
- `ATOM` / `HETATM`: altLoc selection and blanking
- `ANISOU`: kept only if the corresponding ATOM/HETATM line (same serial) is kept

### Automatic skip behavior
By default, if a file contains **no altLoc characters** (all column 17 positions
are blank), the file is **skipped** and no output is written. Use `--force` to
process files regardless of altLoc presence.

## Minimal example

```bash
mlmm fix-altloc -i 1abc.pdb
```

## Output checklist

- `<input>_clean.pdb` -- PDB with altLoc column blanked and duplicates removed
- Original file unchanged (unless `--inplace` is set)

## Common examples

1. Specify output file.

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

2. Process a directory recursively.

```bash
mlmm fix-altloc -i ./structures -o ./cleaned --recursive
```

3. Overwrite input files in-place (creates .bak backups).

```bash
mlmm fix-altloc -i ./structures --inplace --recursive
```

4. Force processing even if no altLoc is detected.

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb --force
```

## Workflow
1. Check if the input file contains any non-blank altLoc characters (column 17).
 - If no altLoc is found and `--force` is not set, skip the file.
2. For each ATOM/HETATM record, build an identity key ignoring the altLoc field:
 - record name, atom name, residue name, chain ID, residue sequence, insertion code, segID
3. Among atoms with the same identity key, select the one with:
 - Highest occupancy (columns 55–60)
 - If tied, the earliest appearance in the file
4. Write output with:
 - Only the selected atoms retained
 - altLoc column (17) blanked to a single space
 - ANISOU records filtered to match retained atoms

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file or directory. | Required |
| `-o, --out PATH` | Output file (if input is a file) or directory (if input is a directory). | File input: `<input>_clean.pdb`; directory input: `<input>_clean/` |
| `--recursive/--no-recursive` | Process `*.pdb` files recursively when input is a directory. | `False` |
| `--inplace/--no-inplace` | Overwrite input file(s) in-place (creates `.bak` backup). | `False` |
| `--overwrite/--no-overwrite` | Allow overwriting existing output files. | `False` |
| `--force/--no-force` | Process files even if no altLoc is detected. | `False` |

## Outputs
- A PDB file with alternate locations removed:
 - File input: `<input>_clean.pdb` by default (when `-o/--out` is omitted)
 - Directory input: `<input>_clean/` directory by default (mirrors subpaths)
 - `OUTPUT.pdb` if `-o/--out` is provided
 - Original file overwritten if `--inplace` is set (backup saved as `<input>.pdb.bak`)

## Integration with `all` workflow
When running the `mlmm all` workflow, `fix-altloc` can be automatically
invoked as a preflight step **after** `add-elem-info` (if element fields were
missing) and **before** pocket extraction. This ensures that:
1. Element symbols are populated first
2. Alternate locations are resolved to a single conformer
3. The extraction step receives clean, unambiguous coordinates

Files are only processed if altLoc characters are detected; otherwise, the
original file is passed through unchanged.

## Handling different atom counts between altLoc states

When different altLoc states contain different atoms (e.g., altLoc A has atoms
N, CA, CB, CG while altLoc B has N, CA, CB, CD), `fix-altloc` handles this correctly:

- **Duplicate atoms** (same residue + atom name in multiple altLocs, e.g., N, CA, CB):
  The best one is selected based on occupancy (highest first, then earliest in file).
- **Unique atoms** (only present in one altLoc, e.g., CG in A, CD in B):
  ALL unique atoms are preserved in the output.

This ensures the output structure contains all atoms from all altLoc states,
with only true duplicates resolved to a single conformer.

**Example:**
```
Input:
 ATOM 1 N AALA A 1... 0.50 # altLoc A
 ATOM 2 CA AALA A 1... 0.50 # altLoc A
 ATOM 3 CG AALA A 1... 0.50 # altLoc A only
 ATOM 4 N BALA A 1... 0.40 # altLoc B
 ATOM 5 CA BALA A 1... 0.40 # altLoc B
 ATOM 6 CD BALA A 1... 0.40 # altLoc B only

Output:
 ATOM 1 N ALA A 1... 0.50 # from A (higher occ)
 ATOM 2 CA ALA A 1... 0.50 # from A (higher occ)
 ATOM 3 CG ALA A 1... 0.50 # kept (A only)
 ATOM 6 CD ALA A 1... 0.40 # kept (B only)
```

## API usage
For programmatic use, the module exports:
```python
from mlmm.fix_altloc import has_altloc, fix_altloc_file

# Check if a file has altLoc
if has_altloc(Path("input.pdb")):
 # Fix altLoc
 was_processed = fix_altloc_file("input.pdb", "output.pdb", overwrite=True)
```

---

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [add-elem-info](add-elem-info.md) -- Repair PDB element columns before altLoc fixing
- [extract](extract.md) -- Extract active-site pocket after altLoc resolution
- [all](all.md) -- End-to-end workflow that includes automatic altLoc fixing
