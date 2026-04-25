# `mlmm oniom-export`

## Purpose

Export an `mlmm-toolkit` system (PDB with B-factor layer labels +
parm7) as a Gaussian g16 ONIOM input file (`.gjf`) or an ORCA ONIOM
input. Useful for hand-comparing to a reference DFT/MM calculation,
re-running with a third-party engine, or feeding into a pipeline
that expects a Gaussian-style input.

The reverse direction is `oniom-import.md`.

## Synopsis

```bash
mlmm oniom-export -i complex.pdb \
    [-p complex.parm7] \
    --engine g16|orca \
    [--high-method 'wB97X-D/def2-svp'] \
    [--low-method 'amber=gaff2'] \
    [-o oniom.gjf]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | PDB with B-factor layer encoding |
| `-p, --parm7` | path | none | Amber parm7 (carries MM parameters; required for `--engine orca`) |
| `--engine` | str | `g16` | Output engine: `g16` (Gaussian) or `orca` |
| `--high-method` | str | `wB97X-D/def2-svp` | Method for the ML/QM region |
| `--low-method` | str | `amber=gaff2` | Method for the MM region (Gaussian convention) |
| `-q` / `-l` / `-m` | — | — | Total charge / per-residue charges / spin |
| `--solvent` | str | none | Implicit solvent (engine-specific naming) |
| `-o, --output` | path | `<input>_oniom.gjf` (or `.inp` for ORCA) | Output path |
| `--two-layer / --three-layer` | flag | `--two-layer` | Use 2-layer (ML+MM) or 3-layer (ML + movable-MM + frozen) ONIOM |

## Examples

### Two-layer ONIOM for Gaussian

```bash
mlmm oniom-export -i complex_layered.pdb -p complex.parm7 \
    --engine g16 \
    --high-method 'wB97X-D/def2-svp' \
    --low-method 'amber=gaff2' \
    -o complex_oniom.gjf
```

### Three-layer for ORCA

```bash
mlmm oniom-export -i complex_layered.pdb -p complex.parm7 \
    --engine orca \
    --high-method 'B3LYP/def2-svp' \
    --three-layer \
    -o complex_oniom.inp
```

## Output

Single text input file the chosen engine can run directly:

- **g16 `.gjf`**: standard Gaussian ONIOM input (`%nproc`, route line
  with `oniom(...)`, charges, coordinates with H/L/M layer suffix).
- **ORCA `.inp`**: ORCA-style multi-layer block (`%multilayer ...`).

The exporter also writes:

```
result.json    # what was exported, layer counts, charge / spin sums
```

## Round-trip

```bash
mlmm oniom-export -i complex.pdb -p complex.parm7 -o oniom.gjf
mlmm oniom-import -i oniom.gjf -o reconstructed.pdb
diff complex.pdb reconstructed.pdb        # should be near-empty
```

The round-trip is lossy in two ways:

1. The PDB element / chain / serial fields are reconstructed from
   the gjf coordinates — atom ordering matches but numbering may
   differ.
2. Custom Amber type assignments are lost (the engine only sees
   element + atomic charge in the gjf).

See `oniom-import.md` for the reverse path.

## Caveats

- The MM force field is **passed by name** in the route line; the
  receiving engine must recognize it (`amber=gaff2` is supported by
  Gaussian g16 ≥ 16.B.01).
- For systems with non-standard ligands, the `parm7`-derived `frcmod`
  parameters get baked into the gjf — review the produced file
  before submitting.
- Three-layer mode requires PDB B-factors to use all three values
  (0.0 / 10.0 / 20.0). If only two are present, the exporter falls
  back to two-layer regardless of `--three-layer`.

## See also

- `oniom-import.md` — reverse direction.
- `../mlmm-structure-io/gjf.md` — Gaussian gjf format reference.
- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — what
  layer label drives the export.
- `../mlmm-structure-io/parm7.md` — parm7 contents needed by exporter.