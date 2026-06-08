# `fix-altloc`

Remove alternate location (altLoc) indicators from PDB files by selecting the best conformer for each atom based on occupancy and dropping duplicates. The altLoc column (column 17, 1-based) is blanked with a single space (a 1-character replacement; no shifting or reformatting), and when the same atom appears in multiple altLoc states the highest-occupancy copy is retained (earliest in file on ties, or when occupancy is missing). `ATOM` / `HETATM` records undergo altLoc selection and blanking; `ANISOU` records are kept only if the corresponding ATOM/HETATM line (same serial) is kept.

## When to use

- Clean a PDB that carries altLoc characters before running downstream ML/MM preparation.
- Repair element columns first with [add-elem-info](add-elem-info.md), then resolve altLocs.
- Files with no altLoc characters are skipped unless `--force` is set.

## Quick examples

```bash
mlmm fix-altloc -i 1abc.pdb
```

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

```bash
mlmm fix-altloc -i ./structures -o ./cleaned --recursive
```

```bash
mlmm fix-altloc -i ./structures --inplace --recursive
```

## Inputs

Command form:

```bash
mlmm fix-altloc -i INPUT [-o OUTPUT] [options]
```

| Input | Required | Notes |
| --- | --- | --- |
| `-i, --input PATH` | yes | Input PDB file or directory. |
| `-o, --out PATH` | no | Output file (if input is a file) or directory (if input is a directory). Defaults to `<input>_clean.pdb` (file) or `<input>_clean/` (directory). |

## Workflow

1. Check if the input file contains any non-blank altLoc characters (column 17).
 - If no altLoc is found and `--force` is not set, skip the file (left unchanged).
2. For each ATOM/HETATM record, build an identity key ignoring the altLoc field:
 - record name, atom name, residue name, chain ID, residue sequence, insertion code, segID
3. Among atoms with the same identity key, select the best one using the occupancy / earliest-appearance rule (highest occupancy; ties or missing occupancy keep the earliest in file; occupancy is read from columns 55–60).
4. Write output with:
 - Only the selected atoms retained
 - altLoc column (17) blanked to a single space
 - ANISOU records filtered to match retained atoms

### Handling different atom counts between altLoc states

When different altLoc states contain different atoms (e.g., altLoc A has atoms
N, CA, CB, CG while altLoc B has N, CA, CB, CD), `fix-altloc` handles this correctly:

- **Duplicate atoms** (same residue + atom name in multiple altLocs, e.g., N, CA, CB):
  The best one is selected using the same occupancy / earliest-appearance rule.
- **Unique atoms** (only present in one altLoc, e.g., CG in A, CD in B):
  ALL unique atoms are preserved in the output.

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

## Outputs

- A PDB file with alternate locations removed:
 - File input: `<input>_clean.pdb` by default (when `-o/--out` is omitted)
 - Directory input: `<input>_clean/` directory by default (mirrors subpaths)
 - `OUTPUT.pdb` if `-o/--out` is provided
 - Original file overwritten if `--inplace` is set (backup saved as `<input>.pdb.bak`)

## Python API

For programmatic use, the module exports:
```python
from pathlib import Path
from mlmm.io.pdb_fix import has_altloc, clean_pdb_file

# Check if a file has altLoc
if has_altloc(Path("input.pdb")):
    # Resolve altLoc into a cleaned PDB (always overwrites output)
    clean_pdb_file(Path("input.pdb"), Path("output.pdb"))
```

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file or directory. | Required |
| `-o, --out PATH` | Output file (if input is a file) or directory (if input is a directory). | File input: `<input>_clean.pdb`; directory input: `<input>_clean/` |
| `--recursive/--no-recursive` | Process `*.pdb` files recursively when input is a directory. | `False` |
| `--inplace/--no-inplace` | Overwrite input file(s) in-place (creates `.bak` backup). | `False` |
| `--overwrite/--no-overwrite` | Allow overwriting existing output files. | `False` |
| `--force/--no-force` | Process files even if no altLoc is detected. | `False` |

The full flag list is in the generated [command reference](reference/commands/index.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [add-elem-info](add-elem-info.md) — Repair PDB element columns before altLoc fixing
- [extract](extract.md) — Extract active-site pocket after altLoc resolution
- [all](all.md) — End-to-end ML/MM workflow (run `fix-altloc` beforehand if your inputs carry altLocs)
