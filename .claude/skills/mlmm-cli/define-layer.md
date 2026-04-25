# `mlmm define-layer`

## Purpose

Assign 3-layer ML / movable-MM / frozen labels by writing the
appropriate B-factor values (0.0 / 10.0 / 20.0) to a PDB. The ML
region is supplied either as a separate PDB or as an atom-index list;
movable-MM is everything within `--radius-freeze` of the ML region
(non-ML atoms inside the radius), and the rest becomes frozen.

Use it before any ML/MM-evaluating subcommand if you want explicit
layer control.

## Synopsis

```bash
mlmm define-layer -i full_system.pdb \
    (--model-pdb model.pdb | --model-indices '1-50,75,100-110') \
    [--radius-freeze 8.0] \
    [-o full_system_layered.pdb]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Full-system PDB |
| `--model-pdb` | path | none | PDB defining the ML-region atoms |
| `--model-indices` | str | none | Comma/range-separated atom indices, e.g. `'1,2,3'` or `'1-10,15,20-25'`. Takes precedence over `--model-pdb`. |
| `--radius-freeze` | float | `8.0` | Distance cutoff (Å) from ML region. Atoms beyond are **frozen** (B-factor 20.0); inside but not ML are **movable-MM** (10.0). |
| `--radius-partial-hessian` | float | `0.0` | Deprecated in 3-layer mode (ignored) |
| `--one-based / --zero-based` | flag | `--one-based` | Interpret `--model-indices` (1- vs 0-based) |
| `-o, --output` | path | `<input>_layered.pdb` | Output PDB with B-factor layer encoding |

`--model-pdb` and `--model-indices` are alternatives; supply exactly
one. `--model-indices` overrides `--model-pdb` if both are given.

## Examples

### Auto-expand around an ML model PDB

```bash
mlmm define-layer -i complex.pdb --model-pdb ml_model.pdb \
    --radius-freeze 8.0 \
    -o complex_layered.pdb
```

### Atom-index list (1-based)

```bash
mlmm define-layer -i complex.pdb --model-indices '1-50,75,100-110' \
    --radius-freeze 6.0 \
    -o complex_layered.pdb
```

### Atom-index list (0-based)

```bash
mlmm define-layer -i complex.pdb --model-indices '0-49,74,99-109' \
    --zero-based -o complex_layered.pdb
```

## Output

A single PDB with the B-factor field set per layer:

```
ATOM      1  CB  TYR A  44       4.050  -8.106   6.935  1.00  0.00           C
                                                              ▲     ▲
                                                              │     └── B-factor: 0.00 (ML)
                                                              └── occupancy unchanged
```

ML atoms get `0.00`, movable-MM (within `--radius-freeze` of any ML
atom) get `10.00`, the rest get `20.00`.

## Caveats

- The radius logic is **freezing-threshold**, not an expansion radius.
  Atoms **beyond** `--radius-freeze` of any ML atom are frozen; atoms
  **inside** but not ML are movable-MM. Increase the radius to free
  more environment, decrease to lock more.
- Layer assignment lives in the **PDB B-factor**, not the parm7. The
  parm7 / rst7 pair is unaffected; downstream subcommands pair the
  layer-encoded PDB with parm7 via `--detect-layer` (default).
- `mlmm extract` does **not** auto-call `define-layer`; if you only
  ran `extract`, run `define-layer` afterwards or pass
  `--model-pdb` / `--model-indices` directly to the consuming
  subcommand.
- `--radius-partial-hessian` is a deprecated 4-layer remnant; ignore
  it in 3-layer mode.

## See also

- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — the
  semantics this command writes.
- `mm-parm.md` — typically run before `define-layer` (parm7 is
  layer-agnostic).
- `extract.md` — extraction of a binding pocket; orthogonal to
  layer assignment.
- Related YAML: `bfactor_ml`, `bfactor_movable_mm`, `bfactor_frozen`,
  `bfactor_tolerance` (live: `import mlmm.defaults as d; print(d.BFACTOR_ML, d.BFACTOR_MOVABLE_MM, d.BFACTOR_FROZEN, d.BFACTOR_TOLERANCE)`).