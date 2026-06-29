# `mlmm add-elem-info`

## Purpose

Repair the element column (PDB cols 77–78) when it is blank or
inconsistent with the atom name. Many PDBs from PyMOL / Maestro come
out with empty element columns; running `extract` on them fails
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
| `-o, --out` | path | optional | Output PDB with element column populated. Omit to overwrite the input in place. |
| `--overwrite / --no-overwrite` | flag | `no-overwrite` | Re-infer the element column even when it is already populated. |

## Examples

```bash
mlmm add-elem-info -i raw.pdb -o cleaned.pdb
```

## Algorithm

The element is inferred from the **atom name and residue name** (with
the HETATM flag), via `guess_element()`. The atom name is read with
Biopython (`atom.get_name()`); the existing element field occupies
cols 77–78. Inference follows a fixed priority:

1. **Ion residues** (`ION` table): polyatomic ions (NH4, H3O+, …) are
   decided per atom name (H/D→H, N→N, O→O); monatomic metals/halogens
   are taken from the residue name (with CL/BR/I/F atom-name fallback).
2. **Protein / nucleic-acid / water** polymers: by convention —
   H/D→H, water→O (non-H), Se→Se, first letter P/N/O/S, then C*→C.
3. **Other ligands/cofactors**: by atom-name prefix — H/D→H,
   C*→C (excluding CL), P*→P.
4. **Fallback normalization**: 2-letter-then-1-letter match against
   the known-element set (deuterium D→H); returns nothing if still
   ambiguous.

## Caveats

- Existing element columns are **preserved** by default; pass
  `--overwrite` to re-infer them. Diff before committing.
- Atom names that don't follow the standard convention (e.g.
  exotic ligand names) may be misclassified; verify by spot-check.

## See also

- `fix-altloc.md` — typically run together as a "PDB cleanup" step.
- `extract.md` — depends on a well-formed element column.
- `mlmm-structure-io/pdb.md` — PDB column layout reference.
