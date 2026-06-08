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
| `-o, --out-prefix` | path | input stem | Output prefix (writes `<prefix>.xyz` and `<prefix>_layered.pdb`) |
| `--help-advanced` | flag | — | Reveal advanced flags |

Use `--mode` to select the engine. XYZ is always written.

## Examples

### Gaussian gjf, no reference PDB

```bash
mlmm oniom-import -i oniom.gjf -o reconstructed
# → reconstructed.xyz, reconstructed_layered.pdb
```

The output PDB has each atom assigned generic residue/chain
identities; B-factors carry the three-layer encoding (ML / movable-MM /
frozen-MM → 0.0/10.0/20.0) derived from each atom's `H`/`L` marker plus
its per-atom freeze flag in the gjf.

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
<out_prefix>_layered.pdb    # PDB with B-factor layer encoding (0.0/10.0/20.0)
```

## Common gotchas

| Symptom | Cause / fix |
|---|---|
| All atoms end up in residue `MOL` | No `--ref-pdb` was given; supply one to recover residue identities. |
| Atom-count mismatch with `--ref-pdb` | Reference PDB has a different ordering; re-export from the original PDB. |
| `ClickException` on a coordinate row | Gjf was a non-ONIOM Gaussian input; without `H`/`L` layer markers the importer cannot parse the coordinate rows and raises an error rather than emitting layers. |
| Mode guess wrong | Pass `--mode g16` or `--mode orca` explicitly. |

## See also

- `oniom-export.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what the
  importer writes.