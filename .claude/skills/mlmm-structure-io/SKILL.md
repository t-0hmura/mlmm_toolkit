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
| **PDB** | atom name, residue, chain, occupancy, **B-factor (= layer label)**, element | Initial input; the B-factor field encodes ML / movable-MM / frozen layer |
| **XYZ** | element + Cartesian coordinates only | Trajectories, post-IRC outputs, single-stage exchange between subcommands |
| **GJF** | element + coords + charge / spin / route line | Round-tripping with Gaussian, including `mlmm oniom-export` / `oniom-import` output |
| **parm7 / rst7** | Amber topology + coordinate pair | MM region force-field parameters; output of `mlmm mm-parm` |

All three "single-file" formats use Å for coordinates and the
conventional periodic-table element symbols. `parm7` is the binary-ish
Amber topology format (text but byte-aligned).

Per-format details are in:

| File | Topic |
|---|---|
| `pdb.md` | PDB column-by-column layout, residue selectors, **B-factor layer encoding (0.0=ML / 10.0=movable-MM / 20.0=frozen)**, link-H placement |
| `xyz.md` | XYZ format, ASE extension comment line |
| `gjf.md` | Gaussian gjf header (`%link0 → route → charge spin → coords`) |
| `parm7.md` | Amber `parm7` topology + `rst7` coordinates (mlmm-specific) |
| `charge-multiplicity.md` | Deciding `-q` and `-m` for an unfamiliar substrate (literature lookup workflow) |

## Decision tree: which format to feed `mlmm-toolkit`

```
Is the input a fresh extraction from the PDB Bank or a model from PyMOL/Maestro?
  └── PDB. mlmm reads residue names directly; use -l 'RES:Q,...'
          to assign per-residue ligand charges.

Is the input a single-segment optimized structure (TS candidate, IRC endpoint)?
  └── XYZ is fine. Pass -q TOTAL_CHARGE and -m MULT explicitly.
      Or use --ref-pdb pointing back to the original PDB so -l still works.

Is the input a Gaussian gjf (with route line, charge, spin in header)?
  └── GJF. mlmm parses the header automatically; -q and -m
      are inferred unless you override.
```

## Editing approach (agent-side)

When an agent must edit a structure file, the basic posture is:

1. **Read the file first** to understand current layout (residues, atom
   counts, charge / multiplicity if present).
2. **Identify the change** and confirm it does not violate format
   conventions (e.g. PDB column widths, XYZ first-line atom count).
3. For unknown charge / multiplicity values, **confirm with the user
   or do a literature lookup** before guessing — see
   `charge-multiplicity.md` for the workflow.
4. Make the smallest possible edit (single residue rename, single
   charge change). Avoid wholesale rewrites.

> Future expansion: this skill is intentionally a base layer.
> Subsequent rounds will add literature-database integration
> (PubChem / ChEBI / UniProt) and ML-based charge inference.

## Subcommand × format compatibility

| Subcommand | PDB | XYZ | GJF |
|---|---|---|---|
| `extract` | ✓ (input + output) | — | — |
| `path-search` | ✓ | ✓ | ✓ |
| `path-opt` | ✓ | ✓ | ✓ |
| `opt` | ✓ | ✓ | ✓ |
| `tsopt` | ✓ | ✓ | ✓ |
| `freq` | ✓ | ✓ | ✓ |
| `irc` | ✓ | ✓ | ✓ |
| `dft` | ✓ | ✓ | ✓ |
| `scan`, `scan2d`, `scan3d` | ✓ | — | — |
| `all` | ✓ | ✓ (single segment) | ✓ |
| `bond-summary` | ✓ | ✓ | — |

If you pass an XYZ to a subcommand that needs residue context (e.g.
`-l 'GLU:-1'`), supply `--ref-pdb <path>` so the residue mapping can be
recovered.

## Quick reference: which fields where

```
PDB  ATOM/HETATM record
     ┌───────────────────────────────────────────────────────┐
     │ name(13-16) altloc(17) resName(18-20) chainID(22)     │
     │ resSeq(23-26)  X(31-38)  Y(39-46)  Z(47-54)           │
     │ occupancy(55-60) bfactor(61-66) element(77-78)        │
     └───────────────────────────────────────────────────────┘

XYZ  line 1: <natoms>
     line 2: <comment, optional ASE Properties=…>
     line 3+: <element>  <x>  <y>  <z>

GJF  %nproc=...  %mem=...
     # <route line:  functional/basis  options>

     <title>

     <charge> <spin>
     <element>  <x>  <y>  <z>
     ...

     [optional connectivity, ECP blocks, ...]
```

Full byte-by-byte / per-keyword detail: see `pdb.md`, `xyz.md`, `gjf.md`.

## Charge / multiplicity defaults

- `-m 1` (singlet, closed shell) is the default for almost every
  organic / biological cluster.
- Use `-m 2` for radicals, `-m 3+` for unusual high-spin metal centers.
- `-q` (total charge) must be explicitly given for XYZ inputs (XYZ has
  no header). `-l 'RES:Q'` derives `-q` for PDB / GJF inputs from
  per-residue charges plus `mlmm-toolkit`'s amino-acid table.

If you're not sure about charge or spin, do **not** guess silently —
follow `charge-multiplicity.md`.

## See also

- `mlmm-cli/extract.md` — residue selectors and link-H caps.
- `mlmm-cli/SKILL.md` — common flag conventions across
  subcommands.
- `mlmm-workflows-output/SKILL.md` — what comes out of the
  pipeline (also XYZ / PDB).