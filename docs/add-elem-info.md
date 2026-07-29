# `add-elem-info`

`mlmm add-elem-info` adds or repairs PDB element symbols (columns 77-78). It infers elements from fixed-column atom names and residue context, and replaces only the element field on ATOM/HETATM records. Use it before downstream tools when element columns are missing or unreliable; `--overwrite` also replaces existing element fields.

## Examples

Command form:

```bash
mlmm add-elem-info -i INPUT [-o OUTPUT] [--inplace] [--overwrite]
```

Add or repair element columns into the non-destructive default output:

```bash
mlmm add-elem-info -i 1abc.pdb
```

This writes `1abc_add_elem.pdb`. To replace the input explicitly:

```bash
mlmm add-elem-info -i 1abc.pdb --inplace
```

Write the result to a separate output file:

```bash
mlmm add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

Re-infer and overwrite existing element fields:

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite
```

## Workflow

1. Read raw PDB records and classify atoms with the residue definitions used
   in `extract.py` (`AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
    and whether the record is HETATM:
 - **Ion residues:** Prefers residue-derived elements; polyatomic ions
  (e.g., NH4, H3O+) are assigned per atom (H/N/O).
 - **Proteins, nucleic acids, water:** Maps H/D to H; water atoms to O/H;
  first-letter mapping for P/N/O/S; recognizes Se; carbon labels
  (CA/CB/CG/...) to C.
 - **Ligands/cofactors:** Uses atom-name prefixes (C*/P*, excluding CL) and
  two-letter/one-letter normalization; recognizes halogens (Cl/Br/I/F).
3. Replace only columns 77–78 on ATOM/HETATM records and preserve all other
   columns and records:
 - No `-o/--out` given: writes `<input>_add_elem.pdb`.
 - `--inplace` without `-o/--out`: replaces the input file.
 - `-o/--out` given: writes to the specified path.
4. Print a summary reporting total atoms, newly assigned, kept existing,
    overwritten (when `--overwrite`), per-element counts, and up to 50
    unresolved atoms (model/chain/residue/atom/serial).

## Outputs

- PDB file with element columns (77-78) populated or corrected
- Console report with totals for processed/assigned atoms, per-element counts, and up to 50 unresolved atoms

## CLI options

| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output PDB path; takes precedence over `--inplace`. | _None_ → `<input>_add_elem.pdb` |
| `--inplace/--no-inplace` | Replace the input file when `-o/--out` is omitted. | `False` |
| `--overwrite/--no-overwrite` | Re-infer and overwrite element fields even if already present (by default, existing values are preserved). | `False` |

Every input line is preserved byte-for-byte except columns 77–78 of
ATOM/HETATM records selected for repair. HEADER, REMARK, CONECT, ANISOU, and
legacy charge columns are retained.

The full flag list is in the generated [command reference](reference/commands/index.md).

## See Also

- [Common Error Recipes](recipes-common-errors.md) — Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) — Detailed troubleshooting guide

- [mm-parm](mm-parm.md) — Build AMBER topology (requires correct element columns)
- [extract](extract.md) — Extract active-site pocket from protein-ligand complex
