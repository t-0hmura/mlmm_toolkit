# `mlmm mm-parm`

## Purpose

Generate Amber `parm7` (topology) + `rst7` (coordinates) + a LEaP-
exported PDB from an input PDB by driving AmberTools (`tleap`,
`antechamber`/`parmchk2` for non-standard ligands). The
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
| `--ff-set` | choice | `ff19SB` | Force-field set: `ff19SB` (with OPC3 water) or `ff14SB` (with TIP3P) |
| `--add-h / --no-add-h` | flag | `--no-add-h` | Add hydrogens via PDBFixer at `--ph` |
| `--ph` | float | `7.0` | pH for PDBFixer hydrogen placement |
| `--add-ter / --no-add-ter` | flag | `--add-ter` | Insert TER records before/after target residues |
| `--keep-temp / --no-keep-temp` | flag | `--no-keep-temp` | Keep the LEaP temp directory for debugging |
| `--out-prefix` | str | input PDB stem | Output prefix (`<prefix>.parm7`, `<prefix>.rst7`, `<prefix>.pdb` when a prefix is set; `<stem>_parm.pdb` only when `--out-prefix` is omitted with `--add-h`) |
| `--help-advanced` | flag | — | Reveal advanced flags |

Force field selection is a single `--ff-set` choice; antechamber/GAFF2
invocation for non-standard residues is automatic when needed.

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
complex.pdb                # LEaP-exported PDB
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

## Antechamber pitfalls (sulfonium / phosphate / odd-electron residues)

`antechamber -c bcc` calls `sqm` for AM1-BCC charge fitting. `sqm` requires
**closed-shell** electron count for `--ligand-mult 1`. If your `-l` value
disagrees with the actual protonation H-count, `sqm` aborts with:

```
Info: Total number of electrons: <N>; net charge: <q>
Info: The number of electrons is odd (<N>).
... Fatal Error! Cannot properly run "sqm".
```

Pre-flight (≤ 30 s, login node):

```bash
awk '$4=="<RES>" && substr($0,77,2)~/H/' input.pdb | wc -l   # H count
# Σ(Z_atoms) - q must be even for closed-shell sqm
```

Common collisions:

| Residue | H count → likely net charge |
|---|---|
| SAM (sulfonium S+) | 22 H → 0 (NH2/COO⁻/S⁺ zwitterion); 23 H → +1 (NH3⁺/COO⁻/S⁺ canonical biological) |
| DMAPP (diphosphate −3) | 9 H → −3 closed-shell; otherwise odd |
| ATP / ADP | 12-13 H → −4 / −3 canonical; verify per dataset |

Do not auto-protonate to "fix" parity — respect the dataset's protonation choice and adjust `-l`.

Copy-paste recovery (match `-l` to the PDB's actual SAM H count):

```bash
# 22 H (neutral zwitterion; e.g. Kulik 2016 COMT model):
mlmm mm-parm -i in.pdb -l 'SAM:0,...' --out-prefix system
# 23 H (canonical biological cation; e.g. bezA / Tsutsumi 2022):
mlmm mm-parm -i in.pdb -l 'SAM:1,...' --out-prefix system
```

(`mlmm mm-parm` / `mlmm all` now run an electron-count preflight that names
the residue and the SAM 22 H = 0 / 23 H = +1 rule before antechamber.)

## See also

- `../mlmm-structure-io/parm7.md` — what `parm7` / `rst7` actually contain.
- `../mlmm-install-backends/ambertools.md` — install AmberTools.
- `define-layer.md` — assign ML / movable-MM / frozen labels on the
  PDB after `mm-parm` (the parm7 itself is layer-agnostic).
- `extract.md` — extracts a binding pocket; orthogonal to parm
  generation.