# `mlmm add-elem-info`

## Purpose

Repair the element column (PDB cols 77–78) when it is blank or
inconsistent with the atom name. Many PDBs from PyMOL / Maestro come
out with empty element columns; running `extract` on them then fails
because the element-aware truncation logic can't classify atoms.

Always run `add-elem-info` (and `fix-altloc.md`) on a freshly-downloaded
RCSB PDB before any other subcommand.

## Synopsis

```bash
mlmm add-elem-info -i in.pdb -o out.pdb
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB |
| `-o, --output` | path | required | Output PDB with element column populated |

## Examples

```bash
mlmm add-elem-info -i raw.pdb -o cleaned.pdb
```

## Algorithm

The element is inferred from the **atom name** (cols 13–16) using
the standard PDB convention: the leading 1–2 alphabetic characters
(after stripping leading whitespace) are the element symbol. For
4-character names where the first character is a digit (e.g. `1HG2`),
the second character is taken instead.

Special cases (Mg, Mn, Fe, Zn, Ca, …) are handled by an internal
case-sensitivity table.

## Caveats

- Existing element columns are **overwritten** by default. Use a
  diff to confirm the changes are sensible.
- Atom names that don't follow the standard convention (e.g.
  exotic ligand names) may be misclassified; verify by spot-check.

## See also

- `fix-altloc.md` — typically run together as a "PDB cleanup" step.
- `extract.md` — depends on a well-formed element column.
- `mlmm-structure-io/pdb.md` — PDB column layout reference.
