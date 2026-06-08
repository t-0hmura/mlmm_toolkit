# `mlmm fix-altloc`

## Purpose

Blank the PDB altLoc column (col 17) without shifting any other field,
and keep **one** altLoc per atom. The default rule is **highest
occupancy first, then earliest appearance** — there is no letter-based
selection. Run on raw RCSB PDBs before `extract`; most downstream tools
(including `mlmm extract`) expect one conformation per atom.

## Synopsis

```bash
mlmm fix-altloc -i in.pdb [-o out.pdb] [--help-advanced]
```

`-i` can be a single PDB or a directory; `-o` matches accordingly
(file → file, directory → directory).

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB file or directory |
| `-o, --out` | path | derived | Output file (if input is file) or directory (if directory). Omitted: writes `<input>_clean.pdb` (file) or `<input>_clean/` (directory). To overwrite the original use `--inplace` (creates a `.bak`). |
| `--recursive / --no-recursive` | flag | — | Recurse into subdirectories (directory input) |
| `--help-advanced` | flag | — | Reveal advanced flags (`--inplace`, `--overwrite`, `--force`) |

The selection rule is **fixed**: highest occupancy, then earliest
appearance.

## Examples

```bash
# Single PDB
mlmm fix-altloc -i raw.pdb -o cleaned.pdb

# Whole directory
mlmm fix-altloc -i raw_pdbs/ -o cleaned_pdbs/

# In place (advanced)
mlmm fix-altloc -i raw.pdb --help-advanced     # see --inplace
mlmm fix-altloc -i raw.pdb --inplace --force
```

## Output

- Single-file mode: one PDB at `-o` (or in place).
- Directory mode: one PDB per input file, mirrored under `-o/`.
- altLoc column (col 17) is blanked for all surviving atoms.

## Caveats

- If the active site has a high-occupancy alt-loc that is **not**
  the chemistry you want, edit the PDB by hand or run a different tool.
- `add-elem-info` and `fix-altloc` are typically run together on raw
  RCSB downloads; repair element columns with `add-elem-info` first,
  then resolve altLocs with `fix-altloc`.

## See also

- `add-elem-info.md` — usually run together.
- `extract.md` — expects altloc-resolved PDB input.
- `../mlmm-structure-io/pdb.md` — col 17 (altLoc) reference.