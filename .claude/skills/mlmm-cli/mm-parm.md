# `mlmm mm-parm`

## Purpose

Generate Amber `parm7` (topology) + `rst7` (coordinates) + a LEaP-
exported PDB from an input PDB by driving AmberTools (`tleap`,
`pdb4amber`, `antechamber`/`parmchk2` for non-standard ligands). The
output is consumed by every downstream subcommand that needs MM
gradients (`opt`, `tsopt`, `freq`, …).

## Synopsis

```bash
mlmm mm-parm -i complex.pdb \
    [-l 'GPP=-3,SAM:1'] \
    [--ligand-mult 'HEM=1'] \
    [--ff-set ff19SB|ff14SB] \
    [--add-h --ph 7.0] \
    [--add-ter / --no-add-ter] \
    [--keep-temp] \
    [--out-prefix complex]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB |
| `-l, --ligand-charge` | str | none | Per-residue charges, e.g. `'GPP=-3,SAM:1'`. Both `=` and `:` separators accepted. |
| `--ligand-mult` | str | none | Per-residue multiplicities, e.g. `'HEM=1,NO:2'` |
| `--ff-set` | choice | `ff19SB` | Force-field set: `ff19SB` (with OPC water) or `ff14SB` (with TIP3P) |
| `--add-h / --no-add-h` | flag | `--no-add-h` | Add hydrogens via PDBFixer at `--ph` |
| `--ph` | float | `7.0` | pH for PDBFixer hydrogen placement |
| `--add-ter / --no-add-ter` | flag | `--add-ter` | Insert TER records before/after target residues |
| `--keep-temp / --no-keep-temp` | flag | `--no-keep-temp` | Keep the LEaP temp directory for debugging |
| `--out-prefix` | str | input PDB stem | Output prefix (`<prefix>.parm7`, `<prefix>.rst7`, optional `<prefix>_parm.pdb` if `--add-h`) |
| `--help-advanced` | flag | — | Reveal advanced flags |

There is **no `--force-field`, `--water`, `--frcmod`, `--charge-method`
flag** (those names belong to a different toolchain). Force field is a
single `--ff-set` choice; antechamber/GAFF2 invocation for non-standard
residues is automatic when needed.

## Examples

### Standard residues only

```bash
mlmm mm-parm -i complex.pdb --out-prefix complex
# → complex.parm7, complex.rst7  in CWD
```

### Non-standard ligand requiring antechamber/GAFF2

```bash
mlmm mm-parm -i complex.pdb \
    --ligand-charge 'GPP=-3,SAM=1' \
    --ff-set ff14SB \
    --add-h --ph 7.0 \
    --out-prefix complex
```

`mm-parm` invokes `antechamber -at gaff2 -c bcc` and `parmchk2`
automatically for residues in `--ligand-charge` whose names are not in
the Amber library, producing GAFF2-typed parameters in the temp dir.
With `--keep-temp` you can inspect the generated `.frcmod` and
`tleap.in` afterwards.

## Output

Files are written **directly to the current working directory** (no
`mm_parm/` subdirectory):

```
complex.parm7              # Amber topology
complex.rst7               # Coordinates (and box if water)
complex_parm.pdb           # LEaP-exported PDB (only when --add-h is True)
```

There is **no `result.json` output**. To verify the parm:

```bash
parmed -p complex.parm7 -i <(echo "summary"; echo "go")
```

## Caveats

- AmberTools must be on PATH (see `../mlmm-install-backends/ambertools.md`).
- The `--ligand-charge` syntax accepts both `=` and `:` (`'GPP=-3,SAM:1'`
  is fine). For consistency with the rest of the toolkit's `-l` flags,
  use `:`.
- Use `--add-h` only if your input PDB lacks hydrogens at the desired
  protonation. Otherwise the PDB is fed to LEaP as-is.
- The B-factor layer encoding (ML / movable-MM / frozen) is **not**
  carried into `parm7`; downstream subcommands pair `parm7` + `rst7`
  with the layer-encoded PDB via `--ref-pdb` or `--detect-layer`.

## See also

- `../mlmm-structure-io/parm7.md` — what `parm7` / `rst7` actually contain.
- `../mlmm-install-backends/ambertools.md` — install AmberTools.
- `define-layer.md` — assign ML / movable-MM / frozen labels on the
  PDB after `mm-parm` (the parm7 itself is layer-agnostic).
- `extract.md` — extracts a binding pocket; orthogonal to parm
  generation.