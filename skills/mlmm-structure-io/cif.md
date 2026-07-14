# mmCIF and oversized PDB inputs

Use `.cif` / `.mmcif` when a structure has multi-character chain IDs,
residue numbers of 10,000 or greater, or identifiers that do not fit PDB
fixed columns. `mlmm-toolkit` converts it to a safely reindexed temporary PDB
for computation and writes `.cif` companions that restore the original IDs.

## Selection forms

Use exact author identifiers when names or numbers repeat:

| Selection | Syntax | Example |
|---|---|---|
| center residue | `CHAIN:RESSEQ[ICODE]` | `enzyme_A:10001B` |
| center name + ID | `CHAIN:RESNAME:RESSEQ[ICODE]` | `enzyme_A:SAM:10001B` |
| scan atom | `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | `enzyme_A:SAM:10001B:CS1` |

`CHAIN:RESNAME` selects all matching residues in that chain. Unqualified
names/numbers can match several residues; use exact forms for production.

## Bridge contract

- Preserve atom order, original chain/residue/insertion identifiers, element,
  atom/residue name, occupancy, B-factor, and formal charge.
- Use the first coordinate model and select one coherent altloc per residue.
- Treat the internal PDB identifiers as temporary; consume the emitted CIF
  when reporting or selecting by original identity.
- Keep the same atoms/order across R/IM/P and in the full-system `parm7`.
- The bridge supports up to 619,938 residues and raises instead of truncating
  beyond that limit.

Most geometry commands accept CIF directly, including `all`, `extract`,
`define-layer`, `sp`, `opt`, `tsopt`, `freq`, `irc`, `dft`, scans, and path
commands. `--ref-pdb` also accepts CIF. Standalone `mm-parm` remains PDB-facing;
`all` bridges before automatic parameterization.

For the full user-facing contract and diagnostics, read `docs/cif.md`.
