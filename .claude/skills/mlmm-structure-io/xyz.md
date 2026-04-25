# XYZ format (xyz.md)

XYZ is the lingua franca of computational chemistry: a simple
plain-text dump of element + Cartesian coordinate, with no residue or
charge information embedded. `mlmm-toolkit` writes XYZ for trajectories,
optimized stationary points, and IRC paths; it accepts XYZ as input
when `-q` / `-m` (or `--ref-pdb`) is supplied.

## Layout

```
<n_atoms>
<comment line — free text, often used by ASE for energies and properties>
<element>  <x>  <y>  <z>
<element>  <x>  <y>  <z>
...
```

- Line 1: integer atom count (no decimals).
- Line 2: arbitrary string. ASE encodes per-frame energies / properties
  here in a `key=value` syntax (see below).
- Lines 3 ... 2 + n_atoms: one atom each. Element is the periodic-table
  symbol (`H`, `C`, `N`, `Mg`, `Zn`, …); coordinates are Å.
- Multiple frames (a trajectory) are concatenated: each frame restarts
  with the atom-count line.

## ASE Properties extension (line 2)

ASE writes XYZ trajectories with a structured comment line:

```
Lattice="..." Properties=species:S:1:pos:R:3 energy=-12345.67 pbc="F F F"
```

Key fields you'll see in `mlmm-toolkit` output:

| Key | Meaning |
|---|---|
| `Properties=species:S:1:pos:R:3` | Standard "element + 3 coords per atom" |
| `energy=...` | Total energy (Hartree by default in mlmm) |
| `forces:R:3` | Per-atom forces appended after coords |
| `Lattice="..."` | Cell vectors, if periodic (cluster models are not periodic) |
| `pbc="F F F"` | Periodic-boundary flags (always F for cluster models) |

`mlmm-toolkit` writes `*_trj.xyz` files this way; the comment line of
each frame carries the energy at that point along the IRC, MEP, or scan.

## Reading and writing in Python

```python
from ase.io import read, write
atoms = read("ts.xyz")          # single frame
trj   = read("mep.xyz", ":")    # all frames as a list
write("out.xyz", trj)           # round-trip
```

Pure-text reading without ASE:

```python
def read_xyz(path):
    with open(path) as f:
        n = int(f.readline())
        comment = f.readline().rstrip()
        coords = [f.readline().split() for _ in range(n)]
    return n, comment, coords
```

## Charge / multiplicity for XYZ inputs

XYZ has **no header for charge or spin**. When passing XYZ to
`mlmm-toolkit`, supply both:

```bash
mlmm tsopt -i ts.xyz -q 0 -m 1 -b uma -o result_tsopt
mlmm dft   -i ts.xyz -q -1 -m 1 --func-basis 'wb97m-v/def2-svp'
```

Or, if a corresponding PDB exists, point at it with `--ref-pdb` so the
residue-based charge mapping (`-l`) still resolves:

```bash
mlmm tsopt -i ts.xyz --ref-pdb cluster.pdb -l 'SAM:1,GPP:-3' -m 1
```

If unsure about `-q` / `-m`, see `charge-multiplicity.md` before
running.

## Common edits

### Strip a frame from a trajectory

```python
from ase.io import read, write
trj = read("mep.xyz", ":")
write("ts.xyz", trj[5])         # frame index 5
```

### Convert XYZ → PDB so that residue-aware tools can consume it

```python
from ase.io import read, write
atoms = read("ts.xyz")
# ASE's PDB writer assigns residue 'MOL' to all atoms — overwrite if needed.
write("ts.pdb", atoms)
```

### Combine multiple stationary points into one trajectory

```bash
cat reactant.xyz ts.xyz product.xyz > rts.xyz
```

This works because each frame already starts with its own atom-count
line; `mlmm trj2fig` consumes the result for energy plots.

## Validation hooks

```bash
# atom count consistent with line 1?
awk 'NR==1{n=$1; expected=n+2} END{if(NR!=expected) print "BAD: line count " NR " expected " expected}' ts.xyz

# any non-element symbols?
awk 'NR>2 && !/^[A-Z][a-z]?[ ]/{print "weird element: " $0}' ts.xyz

# multiple frames? frame count = lines / (natoms+2)
awk 'NR==1{n=$1; per=n+2} END{print "frames:", NR/per}' mep.xyz
```

## See also

- `pdb.md` — the format you'll typically convert XYZ back to when
  you need residue / chain context.
- `charge-multiplicity.md` — how to fill `-q` and `-m` when the XYZ
  doesn't say.
- `mlmm-cli/{trj2fig,energy-diagram}.md` — plotting an XYZ
  trajectory directly.
- `mlmm-workflows-output/SKILL.md` — where the toolkit's
  output XYZ files live.