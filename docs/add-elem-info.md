# `add-elem-info`

`mlmm add-elem-info` adds or repairs PDB element symbols (columns 77-78). It infers elements from fixed-column atom names and residue context, and replaces only the element field on ATOM/HETATM records. Use it before downstream tools when element columns are missing or unreliable; valid existing fields are kept unless `--overwrite-elem` is requested.

## Examples

Command form:

```bash
mlmm add-elem-info -i INPUT [-o OUTPUT] [--overwrite] [--overwrite-elem]
```

Add or repair element columns into the non-destructive default output:

```bash
mlmm add-elem-info -i 1abc.pdb
```

This writes `1abc_add_elem.pdb`. To replace the input explicitly:

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite
```

Write the result to a separate output file:

```bash
mlmm add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

Re-infer and overwrite existing element fields:

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite-elem
```

## Workflow

1. Read raw PDB records and classify atoms with the residue definitions used
   in `extract.py` (`AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
    and whether the record is HETATM:
 - **Ion residues:** Prefers residue-derived elements; polyatomic ions
  (e.g., NH4, H3O+) are assigned per atom (H/N/O).
 - **Proteins, nucleic acids, water:** Maps H/D to H; water atoms to O/H and virtual sites to EP;
  first-letter mapping for P/N/O/S; recognizes Se; carbon labels
  (CA/CB/CG/...) to C.
 - **Ligands/cofactors:** Follows fixed-column atom-name alignment (` NA ` → N,
  `NA  ` → Na), LEaP halogens (` CL1` / ` BR1`), and hydrogen
  labels such as `HG11`.
3. Replace only columns 77–78 on ATOM/HETATM records and preserve all other
   columns and records:
 - No `-o/--out` given: writes `<input>_add_elem.pdb`.
 - `--overwrite` without `-o/--out`: replaces the input file.
 - `-o/--out` given: writes to the specified path; targeting the input requires `--overwrite`.
4. Print a summary reporting total atoms, newly assigned, kept existing,
    overwritten (when `--overwrite-elem`), per-element counts, and up to 50
    unresolved atoms (model/chain/residue/atom/serial).

## Outputs

- PDB file with element columns (77-78) populated or corrected
- Console report with totals for processed/assigned atoms, per-element counts, and up to 50 unresolved atoms

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output PDB path; a separate path takes precedence, while targeting the input requires `--overwrite`. | _None_ → `<input>_add_elem.pdb` |
| `--overwrite/--no-overwrite` | Replace the input file when `-o/--out` is omitted. | `False` |
| `--overwrite-elem/--no-overwrite-elem` | Re-infer valid existing element fields; otherwise only blank or invalid fields are repaired. | `False` |

Every input line is preserved byte-for-byte except columns 77–78 of
ATOM/HETATM records selected for repair. HEADER, REMARK, CONECT, ANISOU, and
legacy charge columns are retained.

The full flag list is in the generated [command reference](reference/commands/index.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [mm-parm](mm-parm.md) — Build AMBER topology (requires correct element columns)
- [extract](extract.md) — Extract active-site pocket from protein-ligand complex
