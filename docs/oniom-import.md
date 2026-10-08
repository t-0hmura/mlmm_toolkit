# `oniom-import` (ONIOM input back to XYZ and a layered PDB)

`oniom-import` **rebuilds an XYZ file and a layered PDB from a Gaussian ONIOM or ORCA QM/MM input file**.

---

## What it is for

* **Inputs edited outside mlmm**: bring a Gaussian or ORCA QM/MM input that you prepared or edited back into mlmm
* **Round trip with names**: restore an input written by [`oniom-export`](oniom-export.md) with the atom and residue names of the original PDB through `--ref-pdb`

---

## Examples

### 1. ORCA input

Rebuild the structure from an ORCA QM/MM input.

```bash
mlmm oniom-import -i ts_guess.inp -o ts_guess_imported
```

The console prints the atom and layer counts and then two `[oniom-import] wrote:` lines for `ts_guess_imported.xyz` and `ts_guess_imported_layered.pdb`.

### 2. Gaussian input, mode from the suffix

The `.gjf` and `.com` suffixes select Gaussian mode, and `.inp` selects ORCA mode.

```bash
mlmm oniom-import -i model.gjf -o model_imported
```

### 3. Set the mode explicitly

`--mode` takes precedence over the suffix and is needed when the file name does not end in `.gjf`, `.com`, or `.inp`.

```bash
mlmm oniom-import -i model.inp --mode orca -o model_imported
```

### 4. Keep names from a reference PDB

Copy the atom and residue names of the PDB that went into `oniom-export` onto the imported coordinates.

```bash
mlmm oniom-import -i model.inp --ref-pdb complex_layered.pdb -o model_imported
```

The log line `[oniom-import] ref_order=identity-verified` shows that the atom order matched the marker written by `oniom-export`.

---

## How it works

1. **Mode**:
`--mode`, or else the suffix of `-i`.
2. **Coordinates and layers**:
In Gaussian mode, each coordinate row reads `<atom> <0|-1> x y z H|L`; `H` atoms are QM, `L` atoms with `0` are movable MM, and `L` atoms with `-1` are frozen MM. In ORCA mode, `QMAtoms` and `ActiveAtoms` come from the `%qmmm` block and the coordinates from the `* xyz` block. QM atoms always count as movable.
3. **QM charge and multiplicity**:
Gaussian: the third and fourth numbers of the six-integer charge and multiplicity line (the QM region at the high level). ORCA: the `* xyz <charge> <multiplicity>` line. Both go into the XYZ comment line as `q=<charge> m=<multiplicity>`.
4. **XYZ**:
`<out_prefix>.xyz` is written with every atom of the input.
5. **Layered PDB**:
`<out_prefix>_layered.pdb` is written with B-factors 0, 10, and 20. Without `--ref-pdb`, each atom is named by its element in residue `MOL` 1 of chain A. With `--ref-pdb`, the reference lines are kept and only the coordinates and B-factors change. The order is checked first: the elements must agree atom by atom; an input with the `oniom-export` marker must also match its hash (`identity-verified`); an input without the marker is accepted when no element appears twice (`element-verified`); otherwise `--allow-unverified-ref-order` is needed (`unverified-opt-in`).

---

## Output files

* **`<out_prefix>.xyz`**: the coordinates of every atom; the comment line reads `mode=<mode> atoms=N qm=… movable=… q=<charge> m=<multiplicity>`.
* **`<out_prefix>_layered.pdb`**: the same coordinates with the layers in the B-factors.
* **Log**: `[oniom-import] mode=…`, the counts `atoms=N, qm=…, movable=…, frozen=…`, and the two `wrote:` lines. With `--ref-pdb`, the log also records `ref_order=identity-verified`, `element-verified`, or `unverified-opt-in`.

---

## Main options

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path | (required) | ONIOM input: `.gjf` / `.com` (g16) or `.inp` (ORCA) |
| `--mode` | `g16` / `orca` | from the `-i` suffix | Format of the input |
| `-o, --out-prefix` | path | the input stem in the current directory | Prefix of the output files |
| `--ref-pdb` | path | `None` | PDB whose atom and residue names are copied to the output; same atoms in the same order |
| `--allow-unverified-ref-order/--no-allow-unverified-ref-order` | flag | `False` | Map `--ref-pdb` by position when neither the marker nor unique elements can confirm the order |

See the [generated CLI reference](reference/commands/oniom_import.md) for every option.

---

## Notes

* **Input files, not output logs**: `oniom-import` reads input files; it does not take optimized geometries from Gaussian or ORCA output. To continue in mlmm, write the final geometry from the external program in the topology atom order and use it with the original parm7 and ML region.
* **Accepted layouts**: Gaussian inputs need two-layer ONIOM coordinate rows with the movable-flag column, as `oniom-export` writes them. ORCA inputs need `QMAtoms {…} end` and `ActiveAtoms {…} end`, each on one line inside `%qmmm`, and a `* xyz <charge> <multiplicity>` block. Other layouts stop with an error.
* **Marker checks**: the marker written by `oniom-export` is checked strictly. A malformed or repeated marker, or a hash that does not match the reference, stops the run even with `--allow-unverified-ref-order`.
* **`--allow-unverified-ref-order`**: use it only for an input without the marker whose repeated elements keep the order from being checked, after you have checked the order yourself. It requires `--ref-pdb`.
* **Exit codes**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [oniom-export](oniom-export.md) — write Gaussian ONIOM or ORCA QM/MM input
* [define-layer](define-layer.md) — the layer B-factors
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
