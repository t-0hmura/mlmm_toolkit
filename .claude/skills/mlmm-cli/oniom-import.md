# `mlmm oniom-import`

## Purpose

Read a Gaussian g16 ONIOM input (or an ORCA ONIOM input) and
reconstruct an `mlmm-toolkit`-shaped PDB with B-factor layer labels.
The reverse of `oniom-export.md`.

Use it to:

- Validate a third-party ONIOM input by round-tripping through the
  toolkit.
- Adapt an existing Gaussian ONIOM workflow into an `mlmm` pipeline.
- Migrate from a hand-built `oniom(...)` calculation to MLIP-driven
  ONIOM.

## Synopsis

```bash
mlmm oniom-import -i oniom.gjf [--engine g16|orca] \
    [--ref-pdb reference.pdb] \
    [-o reconstructed.pdb]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Gaussian `.gjf` (or ORCA `.inp`) ONIOM input |
| `--engine` | str | auto | `g16` or `orca`; auto-detected from the file extension if omitted |
| `--ref-pdb` | path | none | Reference PDB used to recover residue / chain names (the gjf only carries elements) |
| `-o, --output` | path | `<input>_imported.pdb` | Output PDB path |
| `--also-write-xyz` | flag | off | Also write the geometry as `.xyz` |

## Examples

### Minimal (Gaussian gjf, no reference PDB)

```bash
mlmm oniom-import -i oniom.gjf -o reconstructed.pdb
```

The output PDB has each atom assigned residue `MOL` chain `A`
(generic) plus B-factor 0.0 / 10.0 / 20.0 from the gjf's H/M/L
layer suffixes.

### Recover original residue / chain identities

```bash
mlmm oniom-import -i oniom.gjf --ref-pdb original.pdb -o reconstructed.pdb
```

The reference PDB must have the **same atom ordering** as the gjf;
the importer matches by index.

### ORCA input

```bash
mlmm oniom-import -i oniom.inp --engine orca -o reconstructed.pdb
```

## Output

```
reconstructed.pdb    # PDB with ATOM / HETATM and B-factor encoding
result.json          # what was imported (n atoms, n in each layer, charge / spin)
```

`result.json` keys:

```python
import json
d = json.load(open("result.json"))
print(d["status"])
print(d["n_atoms"], d["charge"], d["spin"])
print(d["layers"])              # {"ml": 89, "movable_mm": 412, "frozen": 5104}
print(d["high_method"], d["low_method"])
```

## Common gotchas

| Symptom | Cause / fix |
|---|---|
| All atoms end up in residue `MOL` | No `--ref-pdb` was given; supply one to recover residue identities. |
| Atom count mismatch with `--ref-pdb` | The reference PDB has a different ordering; re-export from the original PDB and try again. |
| B-factors all 0.0 in the output | The gjf was a non-ONIOM Gaussian input; the importer cannot infer layers without the H/M/L tags. |
| Charge doesn't sum correctly | The gjf header's charge/spin pair is for the high layer only; total may differ. Check `result.json`'s `charge` field, which is the **total** layer-summed charge. |

## See also

- `oniom-export.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what
  the importer writes.