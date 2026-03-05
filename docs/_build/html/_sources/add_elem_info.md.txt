# `add-elem-info`

## Overview

> **Summary:** Add or repair PDB element symbols (columns 77-78) using Biopython. Infers elements from atom names plus residue context (proteins, nucleic acids, water, ions, ligands).

### At a glance
- **Use when:** Your PDB has empty element columns that downstream tools (e.g., `mm-parm`, ONIOM export) require.
- **Method:** Biopython `PDBParser` + residue-context heuristics for element assignment.
- **Outputs:** A PDB file with element columns (77-78) populated or corrected.
- **Defaults:** Preserves existing element fields unless `--overwrite` is given. Without `-o`, overwrites the input file.
- **Next step:** [mm-parm](mm_parm.md) or [extract](extract.md).

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
| `--overwrite` | Re-infer and overwrite element fields even if already present (by default, existing values are preserved). | `False` |

## Notes
- For symptom-first diagnosis, start with [Common Error Recipes](recipes_common_errors.md), then use [Troubleshooting](troubleshooting.md) for detailed fixes.

- Only columns 77-78 are modified; coordinates, occupancies, B-factors, charges, altlocs,
 insertion codes, and record ordering stay untouched.
- Existing element fields are detected by scanning the original file's ATOM/HETATM lines
 (serials 7-11, elements 77-78) to reflect the true presence and avoid parser side effects.
- Recognizes standard water/nucleic/protein residue names; treats deuterium "D" as hydrogen "H".
- Depends on Biopython (`Bio.PDB`) and Click.

---

## See Also

- [Common Error Recipes](recipes_common_errors.md) -- Symptom-first failure routing
- [Troubleshooting](troubleshooting.md) -- Detailed troubleshooting guide

- [mm_parm](mm_parm.md) -- Build AMBER topology (requires correct element columns)
- [extract](extract.md) -- Extract active-site pocket from protein-ligand complex
