# `add-elem-info` (repair PDB element columns)

`add-elem-info` **fills in or corrects the element symbols** (columns 77–78) of the ATOM and HETATM records in a PDB file. It infers each element from the atom name and the residue name, keeps element fields that already hold a valid symbol, and changes nothing else in the file.

---

## What it is for

* **PDB files without element columns**: structures from modeling tools or molecular dynamics (MD) that leave columns 77–78 blank.
* **Wrong element symbols**: re-infer every column with `--overwrite-elem`.
* **Preparing input for other commands**: `all` repairs blank element fields by itself with the same rules. Run `add-elem-info` before a standalone command such as `extract`, or when a symbol is wrong rather than blank.

---

## Examples

### 1. Write `<input>_add_elem.pdb`

Fill the element columns and write the result next to the input.

```bash
mlmm add-elem-info -i 1abc.pdb
```

The console prints `[add-elem-info] Wrote: 1abc_add_elem.pdb` and the counts `total atoms`, `assigned/updated`, and `kept existing`. With no `WARNING` line, every atom has an element.

### 2. Choose the output file

```bash
mlmm add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

### 3. Overwrite the input

Replace the input file itself.

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite
```

---

## How it works

1. **Reading the records**:
`add-elem-info` reads every line, follows MODEL blocks, and looks only at ATOM and HETATM records.
2. **Keeping valid symbols**:
An element field that already holds a valid symbol (or `EP` for a water virtual site) is kept. Blank or unrecognized fields are repaired; `--overwrite-elem` re-infers every field.
3. **Inferring the element**:
The element comes from the four-character atom name (columns 13–16) and the residue name. Ion residues give the element of the ion; amino acids, nucleic acids, and water follow the usual naming (H/D, C, N, O, P, S, Se); other ligands use the column alignment of the atom name (` NA ` is N, `NA  ` is Na), including LEaP halogens such as ` CL1` and hydrogen names such as `HG11`.
4. **Writing and the summary**:
Only columns 77–78 of the repaired records change. The console reports the totals, the counts per element, and the atoms that could not be assigned.

---

## Output files

* **The repaired PDB**: `<input>_add_elem.pdb` by default, the path given with `-o`, or the input file itself with `--overwrite` and no `-o`.
* **Console summary**: `total atoms`, `assigned/updated` (atoms whose element column changed), `kept existing`, and `assignment breakdown` (counts per element). Atoms that could not be assigned are left unchanged and listed after `[add-elem-info] WARNING: Could not confidently assign N atoms; left unchanged.`, up to 50 of them. Type the element symbol of each listed atom into columns 77–78 by hand, right-aligned; `--overwrite-elem` uses the same rules and will not assign them.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Input PDB file |
| `-o, --output` | path | `None` | Output PDB file; without it, `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | flag | `False` | Overwrite the input file when `-o` is not given |
| `--overwrite-elem/--no-overwrite-elem` | flag | `False` | Re-infer element columns that already hold a valid symbol |

See the [generated CLI reference](reference/commands/add_elem_info.md) for every option.

---

## Notes

* **What changes**: Every input line is preserved byte-for-byte except columns 77–78 of the ATOM and HETATM records being repaired. HEADER, REMARK, CONECT, ANISOU, and the charge columns (79–80) are written unchanged.
* **Models and chains**: ATOM and HETATM records in every model, chain, and residue are processed.
* **Special names**: deuterium labels become H, selenium (`SE*`) is recognized as Se, and halogens are recognized automatically.
* **Two different flags**: `--overwrite-elem` decides which element columns are re-inferred; `--overwrite` decides only where the file is written.
* **Writing to the input path**: when `-o` names the input file, including through a symbolic link, `--overwrite` is required; without it the run stops with an error.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [extract](extract.md) — extract the active-site model from the repaired PDB
* [all](all.md) — the full workflow, which repairs blank element fields on its own
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
