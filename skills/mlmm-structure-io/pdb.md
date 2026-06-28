# PDB format (pdb.md)

The Protein Data Bank format is the primary input to `mlmm-toolkit`. It
is **column-based** — each field has a fixed character range — so a one-
character drift breaks everything downstream. Treat editing as a
column-aware operation, not as plain-text find-and-replace.

## Record types `mlmm-toolkit` cares about

| Record | Used for |
|---|---|
| `ATOM` | Standard amino-acid / nucleic-acid atoms (residue ≤ 3 letters, in `mlmm-toolkit`'s AMINO_ACIDS table) |
| `HETATM` | Ligand, metal, water, cofactor, link-H atoms |
| `TER` | Chain terminator — used by `extract` to identify chain breaks |
| `END`, `ENDMDL` | File terminator — informational only |
| `CRYST1` | Unit cell — read but not written by `mlmm-toolkit` (cluster model only) |

`mlmm-toolkit` ignores `MODEL`, `ANISOU`, `LINK`, `SSBOND`, etc. Strip
them with `mlmm add-elem-info` or `fix-altloc` if a downstream
step complains.

## Column layout (`ATOM` / `HETATM`)

Columns are 1-based, inclusive on both ends.

| Cols | Field | Width | Type | Example |
|---|---|---|---|---|
| 1–6 | Record name | 6 | left-just `ATOM  ` / `HETATM` | `ATOM  ` |
| 7–11 | Atom serial | 5 | right-just int | `   42` |
| 13–16 | Atom name | 4 | left-just (with 1-char element prefix) | ` CB ` |
| 17 | Alt-loc indicator | 1 | char | ` ` or `A`/`B` |
| 18–20 | Residue name | 3 | upper-case 3-letter code | `SAM` |
| 22 | Chain ID | 1 | char | `A` |
| 23–26 | Residue sequence | 4 | right-just int | `  44` |
| 27 | Insertion code | 1 | char | ` ` |
| 31–38 | X | 8 | right-just float, 3 decimals | `   4.050` |
| 39–46 | Y | 8 | right-just float, 3 decimals | `  -8.106` |
| 47–54 | Z | 8 | right-just float, 3 decimals | `   6.935` |
| 55–60 | Occupancy | 6 | right-just float, 2 decimals | `  1.00` |
| 61–66 | Temperature factor | 6 | right-just float, 2 decimals | `  0.00` |
| 77–78 | Element symbol | 2 | right-just upper-case | ` C` |
| 79–80 | Formal charge | 2 | right-just (e.g. `2+`, `1-`) | `  ` |

`mlmm.add-elem-info` repairs columns 77-78 when they are blank,
which they often are after PyMOL/Maestro export. Always run it before
`extract` if the elements are missing.

## Residue selectors (the `-c / --center` flag)

`mlmm extract` (and any subcommand that also extracts internally)
uses three forms:

```bash
# Form 1 — comma-separated residue names. Matches every residue with
# resName == SAM, GPP, or MG.
mlmm extract -i complex.pdb -c 'SAM,GPP,MG' -o cluster.pdb

# Form 2 — a separate PDB containing only the substrate residues.
mlmm extract -i complex.pdb -c substrate.pdb -o cluster.pdb

# Form 3 — chain-aware: `chainID:resSeq` (numeric resSeq is honored)
mlmm extract -i complex.pdb -c 'A:44' -o cluster.pdb
```

> **Caveat**: only `chainID:resSeq` (numeric resSeq) is parsed as
> chain-aware in `extract.py`. A token like `'A:SAM'` (chain:resName) is
> **not** parsed as chain+resName — `_parse_res_tokens` requires a numeric
> resSeq, so `'A:SAM'` is taken as a literal residue name, matches nothing,
> and raises a "Residue name not found" error. To restrict by chain you
> must supply a numeric resSeq. To select all SAM in chain A, run `extract`
> with `-c 'SAM'` first then trim chains by hand, or use the Form-2
> substrate-PDB workflow.

The pocket radius around the centers is set by `-r <Å>` (default 2.6 Å).
All residues with at least one heavy atom inside the radius are kept.

## Per-residue charge mapping (`-l / --ligand-charge`)

```bash
mlmm extract -i complex.pdb \
    -c 'SAM,GPP,MG' \
    -l 'SAM:1,GPP:-3,MG:2' \
    -o cluster.pdb
```

Standard amino acids are looked up from `mlmm-toolkit`'s internal
`AMINO_ACIDS` table — you only need to provide ligand / metal /
non-standard residues in `-l`. The total cluster charge is the sum of
all residues kept (post extraction).

If you don't know a ligand's formal charge, see
`charge-multiplicity.md` for the lookup workflow.

## Link-hydrogen capping

When `extract` cuts a covalent bond between an in-cluster atom (`A`)
and an out-of-cluster atom (`B`), it places a hydrogen (atom name `HL`)
along the `A→B` direction at 1.09 Å (standard C-H length). The link
hydrogen is written as a `HETATM` named `HL` in residue `LKH` (chain `L`);
see `_format_linkH_block` in `extract.py` for the convention.

Link hydrogens carry **no formal charge**; they do not enter the
charge sum.

Freezing the cap parents is **not** done by `extract` — `extract` only
cuts bonds and adds caps. The ML / movable-MM / frozen layer partition
(which atoms stay fixed during optimization) is assigned later by
`mlmm define-layer`, which encodes the categories in the PDB's B-factor
field for downstream consumption.

## Common edits

### Rename residue 44 in chain A from CYS to CSS

PDB columns: resName 18-20, chainID 22, resSeq 23-26. Match all three
explicitly via `awk substr` (the natural column-aware tool — sed regex
on PDB lines is error-prone because of the dot-counts):

```bash
awk 'BEGIN{OFS=""} \
  ($1=="ATOM" || $1=="HETATM") && substr($0,18,3)=="CYS" \
    && substr($0,22,1)=="A" && substr($0,23,4)+0==44 \
    { $0 = substr($0,1,17) "CSS" substr($0,21) } \
  { print }' my.pdb > my_renamed.pdb
```

For non-trivial edits use Biopython (handles altloc, ANISOU, …).

### Add element column with `add-elem-info`

```bash
mlmm add-elem-info -i my.pdb -o my_with_elem.pdb
```

### Resolve alt-loc

```bash
mlmm fix-altloc -i my.pdb -o my_clean.pdb
```

### Programmatic edits via Biopython

```python
from Bio.PDB import PDBParser, PDBIO
p = PDBParser(QUIET=True).get_structure("x", "my.pdb")
for atom in p.get_atoms():
    if atom.get_name() == "OD1" and atom.get_parent().get_resname() == "ASP":
        atom.set_bfactor(20.0)        # mark frozen, for example
io = PDBIO()
io.set_structure(p)
io.save("my_edited.pdb")
```

## Validation hooks

```bash
# atom count + residue count
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,18,3)}' my.pdb | sort -u

# any missing element columns?
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e=="  ") print NR, $0}' my.pdb

# any duplicate atom names within one residue (often breaks AMBER / parm7)?
awk '/^ATOM|^HETATM/{key=substr($0,22,5)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## See also

- `xyz.md`, `gjf.md` — alternative formats.
- `charge-multiplicity.md` — figuring out per-ligand charges.
- `mlmm-cli/extract.md` — full `extract` flag set and examples.
- `mlmm-cli/{add-elem-info,fix-altloc}.md` — utility subcommands.