# mmCIF and Large Structures

`mlmm-toolkit` accepts `.cif` / `.mmcif` structures and PDB files that exceed
fixed-column atom or residue limits. The numerical workflows still operate on
PDB internally: the input bridge writes a safely reindexed temporary PDB, keeps
the original atom-site metadata, and restores the original chain and residue
identifiers in `.cif` output companions.

## What the bridge preserves

- atom order, element, atom/residue name, coordinates, occupancy, B-factor, and
  formal charge;
- author chain ID, residue sequence ID, and insertion code, including
  multi-character chain IDs and residue numbers of 10,000 or greater;
- the first coordinate model. If an input contains several models, only the
  first is used;
- one coherent alternate-location conformer per residue, selected by occupancy.

The internal PDB identifiers are implementation details. For a bridged input,
use the emitted CIF companion when chain or residue identity matters.

## Supported workflow inputs

The structure bridge is used by `all`, `extract`, `define-layer`, `sp`, `opt`,
`tsopt`, `freq`, `irc`, `dft`, `scan`, `scan2d`, `scan3d`, `path-opt`, and
`path-search`. The same formats are accepted by `--ref-pdb` where that option
is available. Standalone `mm-parm` remains PDB-facing; `all` performs the CIF
bridge before its automatic parameterization stage.

```bash
# End-to-end run from mmCIF
mlmm all -i reactant.cif product.cif \
    -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result

# Preserve topology while using high-precision XYZ coordinates
mlmm tsopt -i hei.xyz --ref-pdb full_system.mmcif \
    --parm full_system.parm7 -q -2 -o result_tsopt
```

For mmCIF or oversized-PDB topology, coordinate conversions write an internal
PDB plus a `.cif` companion with the original identifiers. Multi-frame
trajectories are written as multi-model CIF when conversion is enabled.

## Exact selectors

Use chain-qualified selectors whenever residue names or numbers repeat:

| Context | Exact form | Example |
|---|---|---|
| `extract` / `all -c` by ID | `CHAIN:RESSEQ[ICODE]` | `enzyme_A:10001B` |
| `extract` / `all -c` by name and ID | `CHAIN:RESNAME:RESSEQ[ICODE]` | `enzyme_A:SAM:10001B` |
| scan atom | `CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` | `enzyme_A:SAM:10001B:CS1` |

`CHAIN:RESNAME` intentionally selects every matching residue in that chain.
An unqualified residue number or name may also match more than one residue, so
qualify it for production runs and inspect the selection printed by the CLI.

## Amber topology contract

An Amber `parm7` is positional. It must contain the same full-system atoms in
the same order as every coordinate input; a model PDB is only an unchanged
subset used to select the ML region. `mlmm-toolkit` validates the atom count and
the order of known elements before assigning coordinates and stops on a proven
mismatch. This guard cannot distinguish atoms with the same element, so retain
names, residue IDs, chains, and atom order when preparing R/IM/P structures.

Generate or reuse one topology for the shared atom ordering:

```bash
mlmm mm-parm -i reactant_internal.pdb -l 'SAM:1,GPP:-3' \
    --out-prefix full_system
mlmm all -i reactant.cif product.cif --parm full_system.parm7 \
    -c 'enzyme_A:SAM:10001,enzyme_A:GPP:10002' \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -o result
```

## Limits and diagnostics

- The internal bridge supports up to 619,938 residues (62 internal chains ×
  9,999 residues). It raises an error instead of truncating identifiers beyond
  that limit.
- Coordinates must fit the PDB fixed-column numeric range after conversion.
  Translate an unusually distant structure closer to the origin if requested.
- PDB files with decimal overflow or hybrid-36 serial/residue fields are
  normalized automatically.
- A missing mmCIF `_atom_site.type_symbol`, non-finite coordinates, or a
  different atom count across reaction states is a hard error.

See the “How to construct a reliable `model.pdb`” section in
[Concepts](concepts.md) for ML-region boundary construction and
[CLI Conventions](cli-conventions.md) for charge and scan-selector syntax.
