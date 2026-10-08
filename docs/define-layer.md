# `define-layer` (assign the ML and MM layers)

`define-layer` divides the full system into three layers around the ML region and writes the layer of each atom into the B-factor column of a PDB. The calculation commands read the layers back from the B-factors.

| Layer | B-factor | Atoms | In the calculation |
| --- | --- | --- | --- |
| ML | 0.0 | the ML region | MLIP energy, forces, and Hessian |
| Movable-MM | 10.0 | MM atoms within `--movable-cutoff` (default 8.0 Å) of the ML region | MM, free to move |
| Frozen-MM | 20.0 | MM atoms farther away | MM, coordinates fixed; still part of the MM energy |

## What it is for

* **Preparing the input of the calculation commands**: write the layered full-system PDB that `opt`, `tsopt`, `freq`, and the other commands take with `--parm7`.
* **Changing the movable shell**: widen or narrow the Movable-MM layer with `--movable-cutoff`.
* **Checking the layer sizes**: read the atom count of each layer before a calculation.

`all` runs `define-layer` for you when you pass `-c`.

---

## Examples

### 1. ML region from a model PDB

Give the full system and an ML-region PDB written by `extract` or `all`.

```bash
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb -o labeled.pdb
```

The console prints the atom count of each layer under `Layer Summary`.

### 2. ML region from atom indices

List the ML atoms by index.

```bash
mlmm define-layer -i system.pdb --model-indices "0,1,2,3,4" --zero-based -o labeled.pdb
```

### 3. A wider movable shell

Make every MM residue within 10.0 Å of the ML region movable.

```bash
mlmm define-layer -i system.pdb --model-pdb ml_region.pdb \
    --movable-cutoff 10.0 -o labeled.pdb
```

---

## How it works

1. **ML region**: the atoms of `--model-pdb` are matched to the input by chain, residue number, insertion code, residue name, and atom name; a model atom with a blank chain matches by the other fields and stops with an error when they fit atoms in more than one chain. Without `--model-pdb`, the atoms listed in `--model-indices` form the ML region.
2. **Distances**: for every atom outside the ML region, the shortest distance to an ML atom is computed.
3. **Assignment**: a residue without ML atoms goes into one layer as a whole, by the distance of its closest atom: Movable-MM within `--movable-cutoff`, Frozen-MM beyond. In a residue that contains ML atoms, each non-ML atom is assigned on its own by its distance.
4. **Output**: the input is written again with only the B-factor column changed to 0, 10, or 20, and a Layer Summary is printed.

---

## Output files

```text
./
├─ <input>_layered.pdb   # next to the input when -o is not given
└─ <input>_layered.cif   # mmCIF input, or PDB input too large for the PDB columns
```

The PDB keeps every record of the input and changes only the B-factors. An `-o` path that ends in `.cif` or `.mmcif` is written as a PDB with the same name and `.pdb`. Under `Layer Summary`, the console prints the atom count of each layer on the lines `Layer 1 (ML, B=0):`, `Layer 2 (Movable MM, B=10):`, and `Layer 3 (Frozen MM, B=20):`, and then `Total atoms:`. Color the output by B-factor in a viewer to see the three layers.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | Full-system PDB or mmCIF |
| `--model-pdb` | path | `None` | PDB or mmCIF of the ML-region atoms |
| `--model-indices` | text | `None` | ML atom indices, such as `'1,2,3,4'` or `'1-10,15,20-25'`; 1-based unless `--zero-based`. Used when `--model-pdb` is not given |
| `--movable-cutoff` | float | `8.0` | Distance (Å) from the ML region within which MM atoms are Movable-MM; farther atoms are Frozen-MM |
| `-o, --output` | path | `<input>_layered.pdb` | Output PDB |
| `--one-based/--zero-based` | flag | `--one-based` | How to read `--model-indices` |

See the [generated CLI reference](reference/commands/define_layer.md) for every option.

---

## Notes

* **The ML region is required**: without `--model-pdb` or `--model-indices`, the command stops with exit code 2 and `ERROR: Either --model-pdb or --model-indices must be provided.`
* **`--model-pdb` wins**: when both are given, `--model-pdb` is used, as in the calculation commands.
* **Use the topology-matched PDB**: give the PDB that `mm-parm` writes as `-i`, so that the layered PDB has the same atoms in the same order as the `parm7` ([mm-parm example 4](mm-parm.md#examples)). An atom of `--model-pdb` that is not in the input stops the command with an error.
* **Multi-MODEL input**: only the first MODEL is used, with a warning.
* **Choosing the cutoff**: a smaller `--movable-cutoff` makes the calculation cheaper, and a larger one lets more of the environment relax. A `--movable-cutoff` given to a calculation command replaces the B-factor layers. See [Building the ML region and layers](model-setup.md).

---

## See also

* [Building the ML region and layers](model-setup.md) — choose the ML region, the movable shell, and the Hessian range
* [extract](extract.md) — cut the ML region that `--model-pdb` takes
* [mm-parm](mm-parm.md) — build the topology and the matching PDB to layer
* [all](all.md) — the full workflow; runs `define-layer` with `-c`
* [opt](opt.md) — optimize the layered system
* [Troubleshooting](troubleshooting.md) — layer and atom-order errors
