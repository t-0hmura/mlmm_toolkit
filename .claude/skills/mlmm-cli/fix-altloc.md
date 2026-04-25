# `mlmm fix-altloc`

## Purpose

Resolve PDB alternate locations (`altloc` field, col 17). When a
crystal structure has multiple conformations for the same residue, the
PDB stores them with `altloc='A'`, `altloc='B'`, etc. Most downstream
tools (including `mlmm extract`) expect a single conformation
per atom — `fix-altloc` keeps one and discards the rest.

Run on raw RCSB PDBs before `extract`.

## Synopsis

```bash
mlmm fix-altloc -i in.pdb -o out.pdb [--keep A]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB (with altloc field populated) |
| `-o, --output` | path | required | Output PDB (altloc cleared) |
| `--keep` | str | `A` | Which alt-loc letter to retain (`A`, `B`, …) |
| `--by-occupancy` | flag | off | Keep the conformation with highest occupancy instead of letter-based |

## Examples

```bash
# Keep alt-loc A (most common PDB convention)
mlmm fix-altloc -i raw.pdb -o cleaned.pdb --keep A

# Keep highest-occupancy alt-loc per residue
mlmm fix-altloc -i raw.pdb -o cleaned.pdb --by-occupancy
```

## Caveats

- Choosing the wrong alt-loc can place a side chain in the wrong
  conformation for the chemistry you care about. **Inspect the
  active-site residues** of both alt-loc copies before picking.
- After fix-altloc, col 17 is blanked out. If a downstream tool needs
  it, this won't help.

## See also

- `add-elem-info.md` — usually run together with fix-altloc.
- `extract.md` — depends on cleaned PDB input.
- `mlmm-structure-io/pdb.md` — col 17 (altloc) reference.
