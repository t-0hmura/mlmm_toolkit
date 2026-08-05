# `fix-altloc`

Remove alternate locations by selecting one coherent non-blank altLoc label per
residue. The label with the highest mean occupancy across that residue's
labeled atoms is selected; ties are broken by first appearance. Blank/shared
atoms are retained, atoms from other labels are dropped, and column 17 is
blanked on surviving records. This prevents a per-atom selection from creating
an A/B hybrid that corresponds to no deposited conformer.

## Examples

Command form:

```bash
mlmm fix-altloc -i INPUT [-o OUTPUT] [options]
```

Resolve altLocs in a single file (writes `<input>_clean.pdb`):

```bash
mlmm fix-altloc -i 1abc.pdb
```

Resolve altLocs in a single file with an explicit output name:

```bash
mlmm fix-altloc -i 1abc.pdb -o 1abc_fixed.pdb
```

Process a directory recursively into a new output directory:

```bash
mlmm fix-altloc -i ./structures -o ./cleaned --recursive
```

Process a directory recursively, overwriting files in place:

```bash
mlmm fix-altloc -i ./structures --inplace --recursive
```

## Workflow

1. Check if the input file contains any non-blank altLoc characters (column 17).
 - If no altLoc is found and `--force` is not set, skip the file (left unchanged).
2. Group labeled ATOM/HETATM records by site (chain ID,
   residue sequence, insertion code, and segID).
3. Select one non-blank label per residue using the highest mean parsed
   occupancy (columns 55–60). A label with no parsed occupancy ranks below any
   label with a parsed mean; earliest appearance breaks equal scores, including
   the case where every label lacks parsed occupancy.
4. Keep blank/shared atoms and atoms from the selected label. Resolve malformed
   duplicates that remain by occupancy and file order.
5. Write output with:
 - Only blank/shared atoms and the selected residue conformer retained
 - altLoc column (17) blanked to a single space
 - ANISOU records filtered to match retained atoms

### Handling different atom counts between altLoc states

When different altLoc states contain different atoms (e.g., altLoc A has atoms
N, CA, CB, CG while altLoc B has N, CA, CB, CD), `fix-altloc` processes them as follows:

Only atoms belonging to the selected residue label are retained. An atom unique
to an unselected label is dropped.

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

## Notes

- Files with no altLoc characters are skipped unless `--force` is set.

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [add-elem-info](add-elem-info.md) — Repair PDB element columns before altLoc fixing
- [extract](extract.md) — Extract active-site pocket after altLoc resolution
- [all](all.md) — End-to-end ML/MM workflow (run `fix-altloc` beforehand if your inputs carry altLocs)
