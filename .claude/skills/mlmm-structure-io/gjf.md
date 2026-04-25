# Gaussian gjf format (gjf.md)

`gjf` (or `.com`) is the input format for Gaussian. `mlmm-toolkit`
reads it as an alternative to PDB / XYZ when you have a Gaussian-style
input that already encodes charge, spin, and a route specification.

## Structure

```
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

Layout:

| Block | Lines | Content |
|---|---|---|
| `%`-section | 0+ | Link0 commands: `%nproc`, `%mem`, `%chk`, … |
| Route line | 1 | Starts with `#`, declares method, basis, jobtype |
| (blank) | 1 | Required separator |
| Title | 1 | Free text |
| (blank) | 1 | Required separator |
| Charge / spin | 1 | Two integers separated by space: `<charge> <multiplicity>` |
| Coordinates | n | `<element>  <x>  <y>  <z>` (Å), or with frozen flag `<element> -1 <x> <y> <z>` |
| (blank) | 1 | Terminator |
| (optional) | 0+ | Connectivity, ECP block, custom basis, etc. |

`mlmm-toolkit` reads the geometry block plus the charge/spin line. It
**ignores** the route line (the MLIP backend / DFT engine selected on
the CLI overrides whatever's in the gjf).

## What `mlmm-toolkit` parses

| gjf field | Where it goes |
|---|---|
| Element + x/y/z | Geometry, fed to the calculator |
| Charge | `-q` (only if you don't override on CLI) |
| Multiplicity | `-m` (only if you don't override on CLI) |
| Frozen flag (`-1` in column 2) | `freeze_atoms` |

## CLI usage

```bash
mlmm tsopt -i ts.gjf -b uma -o result_tsopt
```

No need to pass `-q`/`-m` when they're in the header. To override:

```bash
mlmm tsopt -i ts.gjf -q -1 -m 2 -b uma -o result_tsopt
```

`mlmm dft` reads gjf the same way, but its `--func-basis` flag
takes precedence over the route-line method:

```bash
mlmm dft -i ts.gjf --func-basis 'wb97m-v/def2-tzvpd' --engine gpu
```

## Frozen atoms

A `-1` in the second column of a coordinate line marks the atom as
frozen during optimization. Gaussian convention:

```
  C  -1   1.234   5.678   9.012
```

`mlmm-toolkit` reads this as: this atom is not allowed to move during
geometry optimization. The `freeze_atoms` set is otherwise populated
by `extract` (link-atom parents) or by the `--freeze-links` CLI flag.

## Generating gjf from `mlmm-toolkit` output

```bash
# Convert a PDB or XYZ stationary point to gjf
mlmm extract -i input.pdb -o cluster.pdb \
    --convert-files gjf      # writes cluster.gjf alongside cluster.pdb
```

Or via ASE:

```python
from ase.io import read, write
atoms = read("ts.xyz")
atoms.info["charge"] = 0
atoms.info["multiplicity"] = 1
write("ts.gjf", atoms, format="gaussian-in",
      properties=["forces"], extra="# wB97X-D/def2-svp opt freq")
```

## Validation

```bash
# Are charge and spin present?
awk '/^[ ]*-?[0-9]+ +[1-9][0-9]*$/{print "charge,spin:", $0; exit}' ts.gjf

# Atom count
awk 'BEGIN{at=0} /^[ ]*[A-Z][a-z]?[ ]+(-?[0-9]+|[0-9.]+)/{at++} END{print "atoms:", at}' ts.gjf

# Route line
grep -m1 '^#' ts.gjf
```

## Common gotchas

| Symptom | Cause | Fix |
|---|---|---|
| `mlmm-toolkit` says "unknown element symbol Mn"  | Two-letter elements written without correct casing (`MN`) | Use `Mn`, not `MN`, in the element column |
| Charge is read wrong | A blank line missing between title and charge/spin block | Insert the missing blank line |
| Frozen flag ignored | `-1` is in the wrong column or there's a tab instead of space | Use space-separated columns; `-1` must be the 2nd whitespace-separated token |
| Connectivity block confuses the parser | Trailing blocks (after coords) are not parsed by `mlmm-toolkit` and are discarded | Strip them or just ignore the warning |

## See also

- `pdb.md`, `xyz.md` — the other two formats.
- `charge-multiplicity.md` — when the gjf header is missing or wrong.
- `mlmm-cli/dft.md` — gjf is the most natural input for
  `mlmm dft` because it carries charge / spin already.