# `mlmm define-layer`

## Purpose

Assign or update the ML / movable-MM / frozen layer labels on a PDB
by writing the appropriate B-factor values (0.0 / 10.0 / 20.0). Use it
as a standalone "set up the partitioning" step, or after `mm-parm` to
adjust the labels without rebuilding the topology.

## Synopsis

```bash
mlmm define-layer -i complex.pdb -o complex_layered.pdb \
    [--ml '<selector>'] \
    [--movable-mm '<selector>' --mm-radius 4.0] \
    [--frozen-everything-else]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input PDB (existing B-factors are overwritten) |
| `-o, --output` | path | `<input>_layered.pdb` | Output PDB |
| `--ml` | str | required | ML region selector — e.g. `'SAM,GPP'`, `'A:44,A:SAM'`, or a separate PDB path |
| `--movable-mm` | str | `--ml` neighbors | Movable-MM selector (defaults to atoms within `--mm-radius` of `--ml`) |
| `--mm-radius` | float | 4.0 | Pocket radius for auto movable-MM expansion |
| `--frozen-everything-else / --no-frozen-everything-else` | flag | `--frozen-everything-else` | Atoms not in ML or movable-MM become frozen (B-factor 20.0) |
| `--keep-existing` | flag | off | Don't overwrite B-factors that are already 0.0 / 10.0 / 20.0 |
| `--show-summary` | flag | off | Print per-residue layer counts after assignment |

## Examples

### Minimal — auto-expand around ligand

```bash
mlmm define-layer -i complex.pdb \
    --ml 'SAM,GPP,MG' \
    --mm-radius 4.0 \
    -o complex_layered.pdb
```

After this, the PDB has B-factor 0.0 on SAM/GPP/MG atoms,
B-factor 10.0 on residues with at least one heavy atom within 4 Å of
the ML region, and B-factor 20.0 elsewhere.

### Hand-curated layer with explicit movable-MM list

```bash
mlmm define-layer -i complex.pdb \
    --ml 'SAM,GPP,MG' \
    --movable-mm 'A:42,A:43,A:44,A:96,A:97' \
    --no-frozen-everything-else \
    -o complex_layered.pdb
```

Anything not in `--ml` or `--movable-mm` keeps its existing B-factor.

### Promote existing frozen residues to movable-MM

```bash
mlmm define-layer -i complex.pdb \
    --movable-mm 'A:120,A:121' \
    --keep-existing \
    -o complex_layered.pdb
```

## Output

A single PDB with B-factor field set per layer:

```
ATOM      1  CB  TYR A  44       4.050  -8.106   6.935  1.00  0.00           C
                                                              ▲     ▲
                                                              │     └── B-factor: 0.00 (ML)
                                                              └── occupancy unchanged
```

If `--show-summary` is on, stderr lists:

```
Layer counts after define-layer:
  ML (0.0)         : 89 atoms /  4 residues
  Movable MM (10.0): 412 atoms / 28 residues
  Frozen (20.0)    : 5104 atoms / 327 residues
  Total            : 5605 atoms / 359 residues
```

## Caveats

- B-factor is overwritten for matched atoms by default. Use
  `--keep-existing` to preserve manually set values.
- Selector strings use the same grammar as `extract -c`; see
  `../mlmm-structure-io/pdb.md` § "Residue selectors".
- The toolkit's BFACTOR constants (`BFACTOR_ML=0.0`,
  `BFACTOR_MOVABLE_MM=10.0`, `BFACTOR_FROZEN=20.0`,
  `BFACTOR_TOLERANCE=1.0`) are exported:
  ```bash
  python -c "import mlmm.defaults as d; print(d.BFACTOR_ML, d.BFACTOR_MOVABLE_MM, d.BFACTOR_FROZEN, d.BFACTOR_TOLERANCE)"
  ```
- Editing B-factor with `define-layer` does **not** change the
  `parm7`. The MM topology is unaffected; only which atoms the
  toolkit treats as ML / movable-MM / frozen.

## See also

- `../mlmm-structure-io/pdb.md` § "B-factor layer encoding" — the
  semantics this command writes.
- `mm-parm.md` — usually run *before* `define-layer`, but order can
  swap if you regenerate parm7 after relabeling.
- `extract.md` — combines extraction + layering in one shot.