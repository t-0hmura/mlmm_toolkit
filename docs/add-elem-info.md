# `add-elem-info`

## Overview

> **Summary:** Add or repair PDB element symbols (columns 77-78) using Biopython. Infers elements from atom names plus residue context (proteins, nucleic acids, water, ions, ligands).

`mlmm add-elem-info` parses the input PDB with Biopython (`PDBParser`), assigns `atom.element` using residue context and atom-name heuristics, and writes via `PDBIO` to populate columns 77-78. It supports ATOM and HETATM records across all models/chains/residues without altering coordinates.

## Minimal example

```bash
mlmm add-elem-info -i 1abc.pdb
```

## Output checklist

- PDB file with element columns (77-78) populated or corrected
- Console report with totals for processed/assigned atoms, per-element counts, and up to 50 unresolved atoms

## Common examples

1. Write to a specific output file.

```bash
mlmm add-elem-info -i 1abc.pdb -o 1abc_fixed.pdb
```

2. Re-infer and overwrite existing element fields.

```bash
mlmm add-elem-info -i 1abc.pdb --overwrite
```

## Workflow
1. Parse the input PDB with `Bio.PDB.PDBParser`, mirroring the residue
 definitions used in `extract.py` (`AMINO_ACIDS`, `WATER_RES`, `ION`).
2. For each atom, guess the element by combining the atom name, residue name,
 and whether the record is HETATM:
 - **Ion residues:** Prefers residue-derived elements; polyatomic ions
 (e.g., NH4, H3O+) are assigned per atom (H/N/O).
 - **Proteins, nucleic acids, water:** Maps H/D to H; water atoms to O/H;
 first-letter mapping for P/N/O/S; recognizes Se; carbon labels
 (CA/CB/CG/...) to C.
 - **Ligands/cofactors:** Uses atom-name prefixes (C*/P*, excluding CL) and
 two-letter/one-letter normalization; recognizes halogens (Cl/Br/I/F).
3. Write the structure through `PDBIO`:
 - No `-o/--out` given: overwrites the input file.
 - `-o/--out` given: writes to the specified path.
4. Print a summary reporting total atoms, newly assigned, kept existing,
 overwritten (when `--overwrite`), per-element counts, and up to 50
 unresolved atoms (model/chain/residue/atom/serial).

## CLI options
| Option | Description | Default |
| --- | --- | --- |
| `-i, --input PATH` | Input PDB file. | Required |
| `-o, --out PATH` | Output PDB path. When omitted, the input file is overwritten. | _None_ (overwrites input) |
| `--overwrite/--no-overwrite` | Re-infer and overwrite element fields even if already present (by default, existing values are preserved). | `False` |

## See Also

- [Common Error Recipes](recipes-common-errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [mm_parm](mm-parm.md) -- Build AMBER topology (requires correct element columns)
- [extract](extract.md) -- Extract active-site pocket from protein-ligand complex
