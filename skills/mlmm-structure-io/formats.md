# Structure formats

Per-format layouts, edits, and checks for mlmm-toolkit inputs. Choosing a
format, selectors, layers, and the charge rules are in [SKILL.md](SKILL.md).

- [Fields at a glance](#fields-at-a-glance)
- [PDB](#pdb)
- [mmCIF and very large structures](#mmcif-and-very-large-structures)
- [XYZ](#xyz)
- [GJF](#gjf)
- [Amber parm7 and rst7](#amber-parm7-and-rst7)
- [Ligand, ion, and metal charges](#ligand-ion-and-metal-charges)

## Fields at a glance

```text
PDB ATOM/HETATM record (cols 1-based, inclusive)
     name(13-16) altloc(17) resName(18-20) chainID(22)
     resSeq(23-26)  X(31-38)  Y(39-46)  Z(47-54)
     occupancy(55-60)  bfactor(61-66, used as layer: 0.0/10.0/20.0)
     element(77-78)

XYZ  line 1: <natoms>
     line 2: <comment; mlmm writes the energy in hartree>
     line 3+: <element>  <x>  <y>  <z>

GJF  %nproc=...  %mem=...
     # <route line:  functional/basis  options>

     <title>

     <charge> <spin>
     <element>  <x>  <y>  <z>
     ...

parm7  Amber topology; generate with `mlmm mm-parm`, do not hand-edit.
       Pair with rst7 (coordinate snapshot).
```

## PDB

PDB is the main input. It is column-based: each field has a fixed character
range, so a one-character shift corrupts every later field. Edit it with
column-aware tools, not plain find-and-replace.

| Record | Use in mlmm |
|---|---|
| `ATOM` | Standard amino-acid atoms (nucleic-acid residues are treated as unknown ligands) |
| `HETATM` | Ligands, metals, water, cofactors, link H |
| `TER` | Chain terminator; `extract` finds chain breaks from the peptide C–N distance (≤ 1.9 Å), not from `TER` |
| `END`, `ENDMDL` | File terminators |

`CRYST1`, `LINK`, and `SSBOND` are not used. `fix-altloc` resolves altLocs and
drops unselected `ANISOU` records while keeping `MODEL` blocks; it does not
strip `MODEL`, `LINK`, or `SSBOND`.

Columns of `ATOM` / `HETATM`, 1-based and inclusive:

| Cols | Field | Width | Format | Example |
|---|---|---|---|---|
| 1–6 | Record name | 6 | left-justified | `ATOM  ` |
| 7–11 | Atom serial | 5 | right-justified int | `   42` |
| 13–16 | Atom name | 4 | left-justified, 1-char element prefix | ` CB ` |
| 17 | Alt-loc | 1 | char | ` ` or `A`/`B` |
| 18–20 | Residue name | 3 | upper case | `SAM` |
| 22 | Chain ID | 1 | char | `A` |
| 23–26 | Residue number | 4 | right-justified int | `  44` |
| 27 | Insertion code | 1 | char | ` ` |
| 31–38 | X | 8 | float, 3 decimals | `   4.050` |
| 39–46 | Y | 8 | float, 3 decimals | `  -8.106` |
| 47–54 | Z | 8 | float, 3 decimals | `   6.935` |
| 55–60 | Occupancy | 6 | float, 2 decimals | `  1.00` |
| 61–66 | B-factor (layer) | 6 | float, 2 decimals | `  0.00` |
| 77–78 | Element | 2 | right-justified upper case | ` C` |
| 79–80 | Formal charge | 2 | e.g. `2+`, `1-` | `  ` |

`mlmm add-elem-info` fills columns 77–78 when they are blank, as they often are
after a PyMOL or Maestro export; run it before `extract` if elements are
missing.

### Per-residue charge (-l)

```bash
mlmm extract -i complex.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o cluster.pdb
```

Standard amino acids come from the internal `AMINO_ACIDS` table and recognized
monatomic ions from `ION`; list only unknown or non-standard residues in `-l`.
`MG` is already +2: `MG:2` is accepted as a restatement, while `MG:3` does not
override it and logs an unmatched-entry warning. The model charge is the sum
over all retained residues. For an unknown ligand charge, see
[Unknown substrate charge](SKILL.md#unknown-substrate-charge).

### Link hydrogens

With `--add-linkh` (off by default), when `extract` cuts a covalent bond
between a kept atom A and a removed atom B at a carbon boundary, it places a
hydrogen along A→B at 1.09 Å. It is written as `HETATM` atom `HL` in residue
`LKH`, chain `L`, and carries no charge. This is for standalone pocket models;
a `--model-pdb` for mlmm does not need it, because the ML/MM calculator adds
link atoms from the `--parm7` topology.

`extract` only cuts bonds and adds link H; it does not freeze anything. The
ML, Movable-MM, and Frozen-MM layers are assigned later by `mlmm define-layer`
in the B-factor column.

### Common edits

Rename residue 44 of chain A from CYS to CSS by matching resName (18–20),
chainID (22), and resSeq (23–26) with `awk substr`; a sed regex must count
exact character positions and is error-prone:

```bash
awk 'BEGIN{OFS=""} \
  ($1=="ATOM" || $1=="HETATM") && substr($0,18,3)=="CYS" \
    && substr($0,22,1)=="A" && substr($0,23,4)+0==44 \
    { $0 = substr($0,1,17) "CSS" substr($0,21) } \
  { print }' my.pdb > my_renamed.pdb
```

Add element columns or resolve altLocs
([utilities](../mlmm-cli/utilities.md)):

```bash
mlmm add-elem-info -i my.pdb -o my_with_elem.pdb
mlmm fix-altloc -i my.pdb -o my_clean.pdb
```

For non-trivial edits use Biopython, which handles altLoc and `ANISOU`:

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

### Validation checks

```bash
# atom count + residue names
grep -c '^ATOM\|^HETATM' my.pdb
awk '/^ATOM|^HETATM/{print substr($0,18,3)}' my.pdb | sort -u

# any missing element columns?
awk '/^ATOM|^HETATM/{e=substr($0,77,2); if(e=="  ") print NR, $0}' my.pdb

# duplicate atom names within one residue (often breaks Amber / parm7)?
awk '/^ATOM|^HETATM/{key=substr($0,22,5)"-"substr($0,13,4); print key}' my.pdb \
    | sort | uniq -c | awk '$1>1'
```

## mmCIF and very large structures

Use `.cif` or `.mmcif` for chain IDs longer than one character, residue
numbers of 10,000 or more, or identifiers that do not fit the PDB columns.
`all`, `extract`, `define-layer`, `sp`, `opt`, `tsopt`, `freq`, `irc`, `dft`,
`scan`, `scan2d`, `scan3d`, `path-opt`, and `path-search` read them directly,
and so does `--ref-pdb`.

- mlmm reads the first coordinate model and keeps one altLoc per residue, the
  one with the highest mean occupancy.
- During the calculation the atoms carry temporary chain IDs and residue
  numbers; the `.cif` files written next to the outputs restore the original
  chain IDs, residue numbers, and insertion codes. Report and select by the
  original IDs from those files.
- Keep the same atoms and order across R/IM/P and in the full-system parm7.
- Up to 619,938 residues are handled; beyond that the run stops with an error
  instead of truncating.

`mm-parm`, `fix-altloc`, and `add-elem-info` read PDB only; `all` converts an
mmCIF input before it builds the parameters. Details and diagnostics:
[mmCIF and large structures](../../docs/cli-conventions.md#mmcif-and-large-structures).

## XYZ

XYZ holds elements and Cartesian coordinates only, with no residue, charge, or
spin. mlmm writes XYZ for trajectories, optimized stationary points, and IRC
paths.

```text
<n_atoms>
<comment line>
<element>  <x>  <y>  <z>
...
```

- Line 1: integer atom count.
- Line 2: free text. In trajectories written by mlmm, each frame's comment line
  holds its energy in hartree. In extended XYZ from other tools (a comment
  with `Properties=` or `Lattice=`), `energy=` is in eV.
- Following lines: one atom each, element symbol and coordinates in Å.
- A trajectory concatenates frames; each frame starts with its atom-count line.

Read and write with ASE, or plain Python:

```python
from ase.io import read, write
atoms = read("ts.xyz")          # single frame
trj   = read("mep.xyz", ":")    # all frames as a list
write("out.xyz", trj)           # round-trip

def read_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        comment = f.readline().rstrip()
        coords = [f.readline().split() for _ in range(n)]
    return n, comment, coords
```

An XYZ input needs `--parm7`, `--ref-pdb` for the atom order and residue
context, and the charge and multiplicity on the command line:

```bash
mlmm tsopt -i ts.xyz --parm7 real.parm7 --ref-pdb enzyme.pdb -q 0 -m 1 -b uma -o result_tsopt
mlmm dft -i ts.xyz --parm7 real.parm7 --ref-pdb enzyme.pdb -q -1 -m 1 --func-basis 'wb97m-v/def2-svp'
mlmm tsopt -i ts.xyz --parm7 real.parm7 --ref-pdb enzyme.pdb -l 'SAM:1,GPP:-3' -m 1
```

The `--ref-pdb` template also lets `-l` resolve, as in the last line. The run
stops with an error when the XYZ and `--ref-pdb` differ in atom count or in the
element order.

Common edits:

```python
from ase.io import read, write
trj = read("mep.xyz", ":")
write("ts.xyz", trj[5])         # take frame index 5

atoms = read("ts.xyz")
write("ts.pdb", atoms)          # ASE names every residue 'MOL'; overwrite if needed
```

```bash
cat reactant.xyz ts.xyz product.xyz > rts.xyz   # one trajectory of stationary points
```

Concatenation works because each frame starts with its own atom-count line;
`mlmm trj2fig` plots the result ([utilities](../mlmm-cli/utilities.md)).

```bash
# atom count consistent with line 1?
awk 'NR==1{n=$1; expected=n+2} END{if(NR!=expected) print "BAD: line count " NR " expected " expected}' ts.xyz

# any non-element symbols?
awk 'NR>2 && !/^[A-Z][a-z]?[ ]/{print "weird element: " $0}' ts.xyz

# frame count = lines / (natoms+2)
awk 'NR==1{n=$1; per=n+2} END{print "frames:", NR/per}' mep.xyz
```

## GJF

In mlmm, `.gjf` / `.com` is the Gaussian ONIOM exchange format, not a geometry
input. `mlmm oniom-export` writes it from a parm7 and a layered PDB, and
`mlmm oniom-import` reads it (or an ORCA `.inp`) back into an XYZ and a
layer-encoded PDB. The geometry commands do not read gjf; convert it with
`oniom-import` first. Both commands are in [oniom](../mlmm-cli/oniom.md).

```text
%nproc=8
%mem=8GB
%chk=run.chk

# wB97X-D/def2-svp opt freq

  Title (one line, free text)

0 1
  C   0.000   0.000   0.000
  H   0.000   1.090   0.000
  ...

[blank line]
[optional: connectivity table, ECP, basis set, ...]
```

| Block | Lines | Content |
|---|---|---|
| `%` section | 0+ | Link0 commands: `%nproc`, `%mem`, `%chk` |
| Route line | 1 | Starts with `#`; method, basis, job type |
| (blank) | 1 | Required separator |
| Title | 1 | Free text |
| (blank) | 1 | Required separator |
| Charge / spin | 1 | `<charge> <multiplicity>` |
| Coordinates | n | `<element>  <x>  <y>  <z>` (Å), or `<element> -1 <x> <y> <z>` with a frozen flag |
| (blank) | 1 | Terminator |
| (optional) | 0+ | Connectivity, ECP, custom basis |

Charge, spin, and frozen atoms in a calculation come from the command line, not
from a gjf: `-q` / `-m`, and `--freeze-atoms` (1-based indices) or
`geom.freeze_atoms` in a `--config` YAML. A `-1` frozen flag in a Gaussian file
is read only by `oniom-import`.

`--convert-files` (on by default) writes PDB copies of XYZ and trajectory
outputs; it takes no format value, and there is no `--convert-files gjf`. Only
`oniom-export` writes gjf.

## Amber parm7 and rst7

mlmm-toolkit uses the Amber parm7 (topology) and rst7 (coordinates) pair for
the MM part. `mlmm mm-parm` writes one of each from a PDB through AmberTools `tleap`.
Regenerate with `mm-parm` instead of editing a parm7 by hand.

| File | Holds |
|---|---|
| `<name>.parm7` | Atom types, bonds, angles, dihedrals, charges, masses, residue table |
| `<name>.rst7` | Coordinates (optionally velocities) at one geometry |

Commands take the structure with `-i` (PDB/mmCIF, or XYZ with `--ref-pdb`) and
the topology with `--parm7`:

```bash
mlmm opt -i complex.pdb --parm7 complex.parm7 -q 0 -m 1 -o result_opt
```

**Atom order.** The parm7 is positional: its atoms must match the full-system
input one for one. A `model.pdb` is only an unchanged subset that selects the
ML atoms, not a replacement topology. Keep atom names, residues, chains,
insertion codes, and order the same across R/IM/P. Before assigning
coordinates, mlmm compares the atom count and, atom by atom, the element, atom
name, residue name, and residue position with the parm7, and stops with
`Atom-order mismatch between input structure and parm7` on a difference.
Build or reuse the topology from the same ordering; the check does not repair
a reordered structure.

Inspect a parm7, from quick to detailed:

```bash
parmed -p complex.parm7 -i <(echo "summary"; echo "go")            # atom and residue counts
parmed -p complex.parm7 -i <(echo "printDetails @CA"; echo "go")   # atom list
python -c "
import parmed
p = parmed.load_file('complex.parm7', xyz='complex.rst7')
print(f'atoms     = {len(p.atoms)}')
print(f'residues  = {len(p.residues)}')
print(f'bonds     = {len(p.bonds)}')
print(f'box       = {p.box}')
for r in p.residues[:5]:
    print(f'  {r.idx:4d} {r.name:5s}  {len(r.atoms):3d} atoms')
"
```

**Force field.** `mm-parm --ff-set ff19SB` (default) uses OPC3 water;
`--ff-set ff14SB` uses TIP3P, to match an earlier ff14SB/TIP3P setup. The water
model follows the set, and other Amber protein force fields are not offered.
For non-standard residues, `mm-parm` runs `antechamber` to derive GAFF2
parameters; `--keep-temp` keeps the build directory so you can inspect the
generated `.frcmod`.

**Change one residue's layer** without rebuilding the parm7, which does not
change; only the B-factors in the PDB move. This moves residue 44 to Frozen-MM:

```bash
awk 'BEGIN{OFS=""} /^ATOM|^HETATM/{if(substr($0,23,4)+0==44){$0=substr($0,1,60)" 20.00"substr($0,67)}} {print}' \
    complex.pdb > complex_edited.pdb
mlmm opt -i complex_edited.pdb --parm7 complex.parm7 -q 0 -m 1 -o result_opt
```

For larger changes, use `mlmm define-layer`.

Pitfalls:

- `Could not find unit "GPP"`: the ligand is not in the standard Amber
  libraries. `mm-parm` runs `antechamber` on it when the ligand is in the PDB.
- `mismatching atom counts` between parm7 and rst7: the rst7 belongs to a
  different system; regenerate both with `mm-parm`.
- No MM layer is read: the B-factors hold ML atoms only, or are all zero, so
  they are not a layer assignment. A Frozen-MM-only environment is valid;
  re-run `define-layer` only when every MM class is missing.
- The parm7 charge sum differs from the `-l` total: `tleap` rounded charges;
  sum the `charge` column of `printDetails *` in `parmed`.
- `parmed` not on `PATH`: install AmberTools
  ([ambertools.md](../mlmm-install-backends/ambertools.md)).

## Ligand, ion, and metal charges

Common values; always confirm against the mechanism.

| Ligand | Resname | Charge at pH 7 |
|---|---|---|
| S-Adenosylmethionine | `SAM` | +1 |
| S-Adenosylhomocysteine | `SAH` | 0 |
| Geranyl pyrophosphate | `GPP` | −3 |
| ATP | `ATP` | −4 |
| ADP | `ADP` | −3 |
| GTP | `GTP` | −4 |
| NADH | `NAI` / `NDH` | −2 |
| NAD⁺ | `NAD` | −1 |
| FAD | `FAD` | −2 |
| Pyridoxal phosphate | `PLP` | −2 |
| Heme | `HEM` | From the oxidation state, axial ligands, and propionate protonation |
| Phosphate ion | `PO4` | −2 to −3 |

Monatomic ions are summed from the internal `ION` table. `-l` applies only to
residues outside `AMINO_ACIDS`, `ION`, and the water set, so `-l 'MG:3'` or
`-l 'FE:2'` matches no unknown residue, logs a warning, and is ignored.

| Ion | Resname | `ION` value |
|---|---|---|
| Mg²⁺ | `MG` | +2 |
| Zn²⁺ | `ZN` | +2 |
| Mn²⁺ | `MN` | +2 |
| Fe³⁺ | `FE` | +3 |
| Fe²⁺ | `FE2` | +2 |
| Cu²⁺ / Cu⁺ | `CU` / `CU1` | +2 / +1 |
| Na⁺ / K⁺ | `NA` / `K` | +1 |
| Cl⁻ | `CL` | −1 |

When a deposited residue name does not match the intended oxidation state,
rename the residue in the model or give the verified ML-region total with `-q`.

| Metal | Common spin S | Multiplicity (2S+1) |
|---|---|---|
| Mn²⁺ (d⁵) | 5/2 | 6 |
| Fe²⁺ (d⁶), high-spin tetrahedral / weak field | 2 | 5 |
| Fe²⁺ (d⁶), low-spin octahedral / strong field | 0 | 1 |
| Fe³⁺ (d⁵), high-spin | 5/2 | 6 |
| Co²⁺ (d⁷), high-spin | 3/2 | 4 |
| Cu²⁺ (d⁹) | 1/2 | 2 |
| Zn²⁺ (d¹⁰) | 0 | 1 |

## See also

- [SKILL.md](SKILL.md): format choice, selectors, layers, charge and multiplicity.
- [extract](../mlmm-cli/extract.md), [mm-parm](../mlmm-cli/mm-parm.md),
  [define-layer](../mlmm-cli/define-layer.md): the commands that write these files.
- [oniom](../mlmm-cli/oniom.md): gjf export and import.
- [utilities](../mlmm-cli/utilities.md): `add-elem-info`, `fix-altloc`, `trj2fig`.
- [Outputs](../mlmm-overview/outputs.md): where a run writes its XYZ, PDB, and CIF files.
- [ambertools.md](../mlmm-install-backends/ambertools.md): AmberTools for `mm-parm` and `parmed`.
