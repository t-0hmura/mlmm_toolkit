# Amber `parm7` / `rst7` (parm7.md)

`mlmm-toolkit` uses Amber's `parm7` (topology) + `rst7` (coordinates)
pair as the canonical MM-region representation. The pair is produced
by `mlmm mm-parm` (which wraps AmberTools `tleap`) and consumed by
every subcommand that needs MM gradients (`opt`, `tsopt`, `freq`,
`irc`, `path-search`, `path-opt`, `all`, …).

You don't usually edit `parm7` by hand — regenerate via `mm-parm`
when residue / parameter changes are needed. This file is the
reference for **reading** and **diagnosing** them.

## File pair

| File | Purpose |
|---|---|
| `<name>.parm7` | Topology: atom types, bonds, angles, dihedrals, charges, masses, residue table |
| `<name>.rst7` | Coordinates (and optionally velocities) for the system at one geometry |

The names match: a single `mm-parm` invocation produces one `.parm7`
and one `.rst7` file. Downstream subcommands take both:

```bash
mlmm opt -i complex.parm7 complex.rst7 ...
```

## When you need to inspect a `parm7`

Three options, in increasing depth:

```bash
# 1. Quick atom/residue count
parmed -p complex.parm7 -i <(echo "summary"; echo "go")

# 2. Detailed atom list
parmed -p complex.parm7 -i <(echo "printDetails @CA"; echo "go")

# 3. Programmatic via parmed in Python
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

## Layer assignment is **NOT** in `parm7`

`parm7` carries MM parameters. The ML / movable-MM / frozen partition
is encoded in the **PDB's B-factor field** (see `pdb.md`). When a
subcommand takes a `parm7` + `rst7` pair, it pairs them with the
**original PDB** via `--ref-pdb`:

```bash
mlmm opt -i complex.parm7 complex.rst7 --ref-pdb complex.pdb -b uma -o result_opt
```

Without `--ref-pdb`, the toolkit defaults to "all-MM" (no ML region).

If you need to update layer labels, edit the PDB's B-factor column
(see `pdb.md` § "B-factor layer encoding") or use
`mlmm define-layer`.

## Force-field choices

`mlmm mm-parm` accepts:

| `--force-field` | Use case |
|---|---|
| `amber14sb` | Default for protein backbone + side chains |
| `amber99sb-ildn` | Older but widely cited |
| `amber19sb` | Newer; check parameter availability |

For non-standard residues / ligands, `mm-parm` invokes `antechamber`
to derive GAFF2 parameters automatically. Inspect the produced
`.frcmod` file alongside `.parm7` to see what was added.

## Common gotchas

| Symptom | Cause / fix |
|---|---|
| `Could not find unit "GPP"` | Ligand not in standard Amber libraries; `antechamber` must run on it. `mm-parm` does this automatically if the ligand is in the PDB. |
| `mismatching atom counts` between `parm7` and `rst7` | `rst7` was written for a different molecule; regenerate with `mm-parm`. |
| MM-evaluating subcommand reports "no MM region found" | The PDB passed via `--detect-layer` (or `--ref-pdb`) does not have B-factor 10.0 atoms. Re-run `define-layer` to assign movable-MM atoms before passing to `opt`/`tsopt`/etc. (`mlmm extract` itself does **not** consume `parm7` — only PDB. Run `extract` → `mm-parm` → `define-layer` separately.) |
| Charge sum in `parm7` ≠ `-l` total | `tleap` rounded charges; check `parmed` output: `printDetails *` and sum the `charge` column. |
| `parmed` not on PATH | Install via AmberTools (see `../mlmm-install-backends/ambertools.md`). |

## Editing on the fly

If you need to flip a single residue's layer **without rebuilding** the
parm7 (the parm7 is unchanged — only the layer assignment in the PDB
moves):

```bash
# Move residue 44 from movable-MM (10.0) to frozen (20.0):
awk 'BEGIN{OFS=""} /^ATOM|^HETATM/{if(substr($0,23,4)+0==44){$0=substr($0,1,60)" 20.00"substr($0,67)}} {print}' \
    complex.pdb > complex_edited.pdb
mlmm opt -i complex.parm7 complex.rst7 --ref-pdb complex_edited.pdb ...
```

For larger reassignments, use `mlmm define-layer` instead of by-hand
edits.

## See also

- `pdb.md` — B-factor layer encoding (the *partition* that `parm7`
  pairs with).
- `mlmm-cli/mm-parm.md` — full subcommand flag reference.
- `mlmm-cli/define-layer.md` — programmatic layer reassignment.
- `mlmm-install-backends/ambertools.md` — install AmberTools so
  `mm-parm` and `parmed` are available.