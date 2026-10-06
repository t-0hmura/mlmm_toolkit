# `mlmm define-layer`

Writes the three ML/MM layers into the B-factor column of a full-system PDB:
0 for the ML region, 10 for Movable-MM, 20 for Frozen-MM. Run
`mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb`.
Success is a `Layer Summary` with the atom count you expect in each layer.

## When to use

- Before any ML/MM calculation command, when you want explicit control of
  the layers.
- The ML region comes from a model PDB (`--model-pdb`) or an atom-index list
  (`--model-indices`). Non-ML atoms within `--movable-cutoff` (default
  8.0 Å) of the ML region become Movable-MM; the rest become Frozen-MM.

## Minimal run

Around an ML model PDB:

```bash
mlmm define-layer -i complex.pdb --model-pdb ml_region.pdb \
    --movable-cutoff 8.0 -o complex_layered.pdb
```

From an atom-index list, 1-based by default (`--zero-based` for 0-based):

```bash
mlmm define-layer -i complex.pdb --model-indices '1-50,75,100-110' \
    --movable-cutoff 6.0 -o complex_layered.pdb
mlmm define-layer -i complex.pdb --model-indices '0-49,74,99-109' \
    --zero-based -o complex_layered.pdb
```

Without `-o`, the output is `<input>_layered.pdb` next to the input.

## Judge success

- The console prints `Layer Summary` with `Layer 1 (ML, B=0):`,
  `Layer 2 (Movable MM, B=10):`, `Layer 3 (Frozen MM, B=20):`, and
  `Total atoms:`. Check the ML count against the region you meant.
- The output keeps every record of the input and changes only the B-factor:
  `0.00` for ML atoms, `10.00` for Movable-MM, `20.00` for the rest. The
  occupancy is unchanged:

  ```
  ATOM      1  CB  TYR A  44       4.050  -8.106   6.935  1.00  0.00           C
  ```

- Color the output by B-factor in a viewer to see the three layers.

## Pitfalls and recovery

- Give at least one of `--model-pdb` and `--model-indices`; without either,
  the command exits with code 2 and
  `ERROR: Either --model-pdb or --model-indices must be provided.` When both
  are given, `--model-pdb` is used.
- A residue without ML atoms goes into one layer as a whole, by its closest
  atom; in a residue that holds ML atoms, each non-ML atom is assigned on its
  own. Only the first MODEL of a multi-MODEL PDB is used.
- The layers live in the PDB B-factors, not in the `parm7`; the calculation
  commands read them from the PDB (`--detect-layer`, on by default). Use the
  PDB written by `mm-parm` as `-i` so that the atoms match the `parm7`.
- `mlmm extract` does not call `define-layer`. After `extract`, run
  `define-layer`, or pass `--model-pdb` or `--model-indices` directly to the
  calculation command, with `--movable-cutoff` to set the movable shell by
  distance (this replaces the B-factor layers).
- What the model PDB must contain, where to cut the ML/MM boundary, and how
  to choose the cutoff:
  [Set the layers](../mlmm-model-setup/SKILL.md#set-the-layers) and
  [Check the boundary and the charge](../mlmm-model-setup/SKILL.md#check-the-boundary-and-the-charge).

## Next step

- [ML region and layers](../mlmm-structure-io/SKILL.md#ml-region-and-layers):
  what the B-factor values mean to the other commands.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md): grow or trim the ML
  region and the movable shell.
- [mm-parm.md](mm-parm.md): usually run before `define-layer`; the `parm7`
  holds no layers.
- [extract.md](extract.md): cut a binding pocket.
