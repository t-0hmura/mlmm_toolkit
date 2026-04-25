# `mlmm extract`

## Purpose

Cuts an active-site cluster from a PDB around a substrate selection.
Severed covalent bonds are capped with hydrogens (link-H), residue
charges are summed (`-l 'RES:Q'` for non-standard residues), and the
extracted cluster is written as a PDB ready for any other subcommand.

## Synopsis

```bash
mlmm extract -i complex.pdb -c <substrate-spec> [-l 'RES:Q,...'] \
    [-r <radius_A>] [-o cluster.pdb] [--multi-model]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required | Protein–substrate complex PDB(s); multi-input requires identical atom counts |
| `-c, --center` | str | required | Substrate selector: PDB path, `'RES1,RES2'`, `'A:44,B:SAM'`, or residue-name list |
| `-r, --radius` | float | 2.6 | Pocket radius (Å) around `-c` atoms |
| `-l, --ligand-charge` | str | none | Per-residue charges (amino acids derived from internal table) |
| `-q, --charge` | int | derived | Override total cluster charge |
| `-m, --multiplicity` | int | 1 | Spin multiplicity |
| `-o, --output` | path | `model.pdb` | Output PDB path; multiple inputs → per-file or multi-MODEL |
| `--multi-model` | flag | off | Write all inputs into one multi-MODEL PDB |
| `--freeze-links / --no-freeze-links` | flag | `--freeze-links` | Mark link-H parents as frozen (B-factor) |
| `--convert-files` | str | none | Also write `.xyz` / `.gjf` alongside `.pdb` |
| `--show-config` / `--dry-run` | flag | off | Print resolved config, no extraction |
| `--help-advanced` | flag | — | Hidden flags (residue rename, custom AA table, …) |

## Examples

### Minimal — extract around two residues

```bash
mlmm extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 -o cluster.pdb
```

### Multiple structures, identical ordering

```bash
mlmm extract -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o "1.R_clu.pdb" "3.P_clu.pdb"
```

### Substrate-only PDB driving extraction

```bash
mlmm extract -i complex.pdb -c substrate.pdb -r 3.5 -o cluster.pdb
```

## Output

```
cluster.pdb                # extracted cluster, link-H caps applied
result.json (if --out-json) # extraction stats: total charge, n_atoms, n_residues, freeze list
```

The PDB B-factor field marks **frozen atoms** (link-H parents) — read
this back with any subsequent subcommand via `freeze_atoms`.

```python
import json
d = json.load(open("result.json"))
print(d["charge"], d["n_atoms"], d["residues"])
print(d["freeze_atoms"])         # indices kept fixed in optimization
```

## Caveats

- Atom names must match exactly (case-sensitive) when using `'A:44:CA'`-style
  selectors. `add-elem-info` and `fix-altloc` should run before extract
  if the PDB came out of PyMOL / Maestro.
- `-r` < 2 Å usually leaves the cluster missing essential coordinating
  atoms; 3.0–4.5 Å is typical.
- Ligand charges come **only** from `-l`; the internal table covers
  standard amino acids and a small list of common cofactors. For
  uncommon ligands, see `mlmm-structure-io/charge-multiplicity.md`.

## See also

- `mlmm-structure-io/pdb.md` — PDB column layout, residue selectors.
- `add-elem-info.md`, `fix-altloc.md` — pre-clean a raw PDB.
- Defaults: `import mlmm.defaults as d; print(d.GEOM_KW_DEFAULT)`