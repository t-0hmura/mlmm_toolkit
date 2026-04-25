---
name: mlmm-structure-io
description: How to read, edit, and write the structure formats mlmm-toolkit handles — PDB (with B-factor layer encoding), XYZ, GJF, and Amber parm7/rst7 — plus the charge / multiplicity decision workflow.
---

# mlmm-toolkit Structure I/O

## Overview

`mlmm-toolkit` reads four formats; each carries different information
and is preferred for different stages:

| Format | Carries | Preferred for |
|---|---|---|
| **PDB** | atom name, residue, chain, occupancy, **B-factor (layer label)**, element | Initial input; B-factor encodes ML / movable-MM / frozen layer |
| **XYZ** | element + Cartesian coordinates only | Trajectories, post-IRC outputs, single-stage exchange between subcommands |
| **GJF** | element + coords + charge / spin / route line | Round-tripping with Gaussian; `mlmm oniom-{export,import}` |
| **parm7 / rst7** | Amber topology + coordinate pair | MM region force-field parameters; output of `mlmm mm-parm` |

PDB / XYZ / GJF use Å for coordinates and conventional element symbols.
`parm7` is the Amber topology format (text but byte-aligned).

Per-format details:

| File | Topic |
|---|---|
| `pdb.md` | PDB column-by-column layout, residue selectors, **B-factor layer encoding (0.0=ML / 10.0=movable-MM / 20.0=frozen)**, link-H placement |
| `xyz.md` | XYZ format, ASE extension comment line |
| `gjf.md` | Gaussian gjf header (`%link0 → route → charge spin → coords`) |
| `parm7.md` | Amber `parm7` topology + `rst7` coordinates (mlmm-specific) |
| `charge-multiplicity.md` | Deciding `-q` and `-m` for an unfamiliar substrate (literature lookup workflow) |

## Decision tree: which format to feed `mlmm-toolkit`

```
Is the input the full enzyme + parm7 you'll run ML/MM on?
  └── PDB (full enzyme, with B-factor layer assignment) + parm7
      → all of opt / tsopt / scan / path-search / freq / irc / dft

Is the input a single TS candidate to validate?
  └── XYZ + --ref-pdb (full enzyme PDB) + --parm
      → tsopt / freq / irc / all (TS-only mode)

Is the input a Gaussian g16 ONIOM input you want to import?
  └── GJF → mlmm oniom-import → reconstructs PDB + extracts layer info

Do you need a parm7 / rst7 from a raw enzyme PDB?
  └── mlmm mm-parm → PDB + AmberTools tleap → parm7 + rst7
```

## ML/MM-aware CLI conventions

Most subcommands take `--parm FILE` (the parm7) plus one of:

- `--detect-layer` (default on) — read the layer assignment from the
  input PDB's B-factor field
- `--model-pdb FILE` — explicit ML-region PDB
- `--model-indices '1-50,75,100-110'` — explicit atom-index list

When `-i` is XYZ, also pass `--ref-pdb` so atom ordering and residue
context are recoverable.

## Editing approach (agent-side)

When an agent must edit a structure file:

1. **Read the file first** to understand current layout (residues,
   atom counts, B-factor layer assignment, charge/multiplicity if
   present).
2. **Identify the change** and confirm it does not violate format
   conventions (PDB column widths, XYZ first-line atom count,
   parm7 byte alignment).
3. For unknown charge / multiplicity values, **confirm with the user
   or do a literature lookup** before guessing — see
   `charge-multiplicity.md` for the workflow.
4. For layer-assignment changes (B-factor edits), use
   `mlmm define-layer` rather than hand-editing if possible.

> **Note**: this skill is a base layer. Subsequent rounds will add
> literature-database integration (PubChem / ChEBI / UniProt) and
> ML-based charge inference.

## Subcommand × format compatibility

| Subcommand | PDB | XYZ | GJF | parm7 |
|---|---|---|---|---|
| `extract` | ✓ (in/out) | — | — | — |
| `mm-parm` | ✓ (in) | — | — | ✓ (out) |
| `define-layer` | ✓ (in/out) | — | — | — |
| `path-search` / `path-opt` | ✓ | ✓ | ✓ | ✓ |
| `opt` / `tsopt` / `freq` / `irc` | ✓ | ✓ | ✓ | ✓ |
| `dft` | ✓ | ✓ | ✓ | ✓ |
| `scan` / `scan2d` / `scan3d` | ✓ | — | — | ✓ |
| `oniom-export` | ✓ (in) | ✓ (in) | ✓ (out) | ✓ (in) |
| `oniom-import` | ✓ (out) | — | ✓ (in) | — |
| `bond-summary` | ✓ | ✓ | — | — |

## Quick reference

```
PDB ATOM/HETATM record (cols 1-based, inclusive)
     name(13-16) altloc(17) resName(18-20) chainID(22)
     resSeq(23-26)  X(31-38)  Y(39-46)  Z(47-54)
     occupancy(55-60)  bfactor(61-66, used as layer: 0.0/10.0/20.0)
     element(77-78)

XYZ  line 1: <natoms>
     line 2: <comment, optional ASE Properties=…>
     line 3+: <element>  <x>  <y>  <z>

GJF  %nproc=...  %mem=...
     # <route line:  functional/basis  options>

     <title>

     <charge> <spin>
     <element>  <x>  <y>  <z>
     ...

parm7  Amber topology — generate with `mlmm mm-parm`; do not hand-edit.
       Pair with rst7 (coordinate snapshot).
```

Full byte-by-byte / per-keyword detail in the per-format mds.

## Charge / multiplicity defaults

- `-m 1` (singlet, closed shell) is the default for almost every
  organic / biological / metal-coordination cluster.
- Use `-m 2` for radicals, `-m 3+` for unusual high-spin metals.
- `-q` is the **ML region** charge. `-l 'RES:Q'` derives `-q` from
  per-residue charges + `mlmm`'s internal amino-acid table.
- For XYZ inputs (no header), `-q` and `-m` must be on the CLI.

If unsure about charge or spin, do **not** guess silently — follow
`charge-multiplicity.md`.

## See also

- `mlmm-cli/extract.md`, `mm-parm.md`, `define-layer.md` — pre-pipeline.
- `mlmm-cli/SKILL.md` — common flag conventions across subcommands.
- `mlmm-workflows-output/SKILL.md` — what comes out of the pipeline
  (XYZ / PDB).