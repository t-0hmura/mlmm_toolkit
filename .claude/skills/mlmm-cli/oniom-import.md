# `mlmm oniom-import`

## Purpose

Read a Gaussian g16 ONIOM input (or an ORCA ONIOM input) and
reconstruct an `mlmm-toolkit`-shaped XYZ + B-factor-layered PDB. The
reverse of `oniom-export.md`.

Use it to:

- Validate a third-party ONIOM input by round-tripping through the toolkit.
- Adapt an existing Gaussian ONIOM workflow into an `mlmm` pipeline.
- Migrate from a hand-built `oniom(...)` calculation to MLIP-driven ONIOM.

## Synopsis

```bash
mlmm oniom-import -i oniom.gjf [--mode g16|orca] \
    [--ref-pdb reference.pdb] \
    [-o reconstructed]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Gaussian `.gjf`/`.com` (g16) or ORCA `.inp` ONIOM input |
| `--mode` | choice | inferred | `g16` or `orca`; falls back to input-file suffix |
| `--ref-pdb` | path | none | Reference PDB used to recover residue / chain / atom-name metadata (atom count must match) |
| `-o, --out-prefix` | path | input stem | Output prefix (writes `<prefix>.xyz` and `<prefix>.pdb`) |
| `--help-advanced` | flag | — | Reveal advanced flags |

There is **no `--engine` flag** (use `--mode`) and **no
`--also-write-xyz` flag** (XYZ is always written).

## Examples

### Gaussian gjf, no reference PDB

```bash
mlmm oniom-import -i oniom.gjf -o reconstructed
# → reconstructed.xyz, reconstructed.pdb
```

The output PDB has each atom assigned generic residue/chain
identities; B-factors carry the H/M/L → 0.0/10.0/20.0 layer mapping
from the gjf.

### Recover original residue / chain identities

```bash
mlmm oniom-import -i oniom.gjf --ref-pdb original.pdb -o reconstructed
```

The reference PDB must have the **same atom count and ordering** as
the gjf; the importer matches by index.

### ORCA input

```bash
mlmm oniom-import -i oniom.inp --mode orca -o reconstructed
```

## Output

```
<out_prefix>.xyz    # element + coords, ASE-compatible
<out_prefix>.pdb    # PDB with B-factor layer encoding (0.0/10.0/20.0)
```

## Common gotchas

| Symptom | Cause / fix |
|---|---|
| All atoms end up in residue `MOL` | No `--ref-pdb` was given; supply one to recover residue identities. |
| Atom-count mismatch with `--ref-pdb` | Reference PDB has a different ordering; re-export from the original PDB. |
| B-factors all 0.0 in output | Gjf was a non-ONIOM Gaussian input; the importer cannot infer layers without H/M/L tags. |
| Mode guess wrong | Pass `--mode g16` or `--mode orca` explicitly. |

## See also

- `oniom-export.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what the
  importer writes.