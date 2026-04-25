# `mlmm extract`

## Purpose

Extracts a **binding pocket** from a protein-substrate complex PDB
around the specified substrate residues, with biochemically aware
truncation and optional carbon-only link-H. Default output is
`pocket.pdb`. Multi-input mode produces one PDB per file (or one
multi-MODEL PDB) for use with `path-search`-style workflows.

This is mlmm-toolkit's pocket extractor; for layer assignment use
`define-layer.md` (extract does not write B-factor labels).

## Synopsis

```bash
mlmm extract -i complex.pdb -c <substrate-spec> [-l 'RES:Q,...'] \
    [-r 2.6] [-o pocket.pdb] [--add-linkh] [--include-h2o] \
    [--exclude-backbone] [--out-json]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path(s) | required | Protein-substrate complex PDB(s); multi-input requires identical atom count + ordering |
| `-c, --center` | str | required | Substrate selector: PDB path, residue-ID list (`'123,124'`, `'A:123,B:456'`), or residue-name list (`'GPP,SAM'`) |
| `-r, --radius` | float | `2.6` | Cutoff (Å) around substrate atoms for pocket inclusion |
| `--radius-het2het` | float | `0` | Cutoff (Å) for substrate-protein hetero-atom proximity (non-C/H); 0 disables |
| `-l, --ligand-charge` | str | none | Per-residue charges, e.g. `'GPP:-3,SAM:1'` |
| `-o, --output` | path(s) | `pocket.pdb` (single) / `pocket_<filename>.pdb` (multi) | Output PDB(s) |
| `--include-h2o / --no-include-h2o` | flag | `--include-h2o` | Include waters (HOH/WAT/H2O/DOD/TIP/TIP3/SOL) |
| `--exclude-backbone / --no-exclude-backbone` | flag | `--no-exclude-backbone` | Delete main-chain atoms from non-substrate amino acids |
| `--add-linkh / --no-add-linkh` | flag | `--no-add-linkh` | Add carbon-only link-H at 1.09 Å along cut-bond directions |
| `--selected-resn` | str | none | Comma/space-separated residue IDs to force-include |
| `--modified-residue` | str | none | Residue names (+ optional charge) to treat as amino acids, e.g. `'HD1:0,SEP:-2'` |
| `-v, --verbose / --no-verbose` | flag | `--verbose` | INFO-level logging |
| `--out-json / --no-out-json` | flag | `--no-out-json` | Write `result.json` next to the output PDB |
| `--help-advanced` | flag | — | Reveal advanced flags |

There is **no `-q/--charge`, `--multi-model`, `--freeze-links`,
`--convert-files`, `--show-config`, or `--dry-run` flag** in mlmm
extract (these belong to `pdb2reaction extract`). For mlmm, total
charge is derived from `-l` only.

## Examples

### Minimal — extract around two residues

```bash
mlmm extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 -o pocket.pdb
```

### With link-H and JSON metadata

```bash
mlmm extract -i complex.pdb -c 'A:123,A:124' -l 'GPP:-3' -r 3.5 \
    --add-linkh --out-json -o pocket.pdb
# → pocket.pdb + pocket.json (when --out-json)
```

### Multi-input (atom-ordering must match)

```bash
mlmm extract -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o "1.R_pocket.pdb" "3.P_pocket.pdb"
```

## Output

```
pocket.pdb                  # extracted pocket (default name)
pocket.json (if --out-json) # extraction stats: charge, n_atoms, n_residues
```

The PDB does **not** contain B-factor layer encoding by default — run
`define-layer` afterwards to assign ML / movable-MM / frozen labels.

```python
import json
d = json.load(open("pocket.json"))
print(d["charge"], d["n_atoms"], d["residues"])
```

## Caveats

- `--add-linkh` is **off by default** in mlmm (in contrast to pdb2reaction).
- `--include-h2o` is **on by default**; turn off with `--no-include-h2o`
  for dry pockets.
- `--exclude-backbone` is **off by default**; turn on for cluster-style
  truncated backbones.
- For the standard ML/MM workflow, follow `extract → mm-parm →
  define-layer → opt/tsopt/...`.
- Atom names must match exactly (case-sensitive). Run
  `mlmm add-elem-info` / `fix-altloc` first if the PDB came out of
  PyMOL or Maestro.

## See also

- `../mlmm-structure-io/pdb.md` — PDB column layout, residue selectors.
- `mm-parm.md` — parm7 generation (typically run after extract).
- `define-layer.md` — assign ML / movable-MM / frozen labels (B-factor).
- `add-elem-info.md`, `fix-altloc.md` — pre-clean a raw PDB.