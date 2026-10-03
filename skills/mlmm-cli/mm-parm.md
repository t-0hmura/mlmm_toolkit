# `mlmm mm-parm`

Builds the Amber topology of the full system with AmberTools (tleap, plus
antechamber and parmchk2 for unknown residues). Run
`mlmm mm-parm -i complex.pdb -l 'GPP:-3,SAM:1' --out-prefix system`. Success
is `[mm-parm] Wrote: …` for `system.parm7`, `system.rst7`, and `system.pdb`.

## When to use

- Every calculation command except `all` needs `--parm7`; `all` runs
  `mm-parm` itself unless you pass one.
- Build it once per system and reuse it for R, intermediates, and P, which
  have the same atoms in the same order.

## Minimal run

Standard residues only:

```bash
mlmm mm-parm -i complex.pdb --out-prefix system
```

A non-standard ligand, with hydrogens added at pH 7 and ff14SB:

```bash
mlmm mm-parm -i complex.pdb --ligand-charge 'GPP:-3,SAM:1' \
    --ff-set ff14SB --add-h --ph 7.0 --out-prefix system
```

Each residue that tleap does not know is parameterized with
`antechamber -at gaff2 -c bcc` and `parmchk2`, using its charge from `-l` (0
if not given) and its multiplicity from `--ligand-mult` (1 if not given).
`--ff-set` is `ff19SB` (default, OPC3 water) or `ff14SB` (TIP3P water).
`--keep-temp` keeps the working directory `parm7build_*`, with the `.frcmod`
and the tleap logs.

## Judge success

- The files are written to the current directory, with no subdirectory and no
  `result.json`: `<prefix>.parm7` (topology), `<prefix>.rst7` (coordinates,
  and the box if there is water), and `<prefix>.pdb`.
- `<prefix>.pdb` has the same atoms in the same order as the `parm7`, with
  element columns filled. It is written when `--out-prefix` is given, or as
  `<input>_parm.pdb` with `--add-h` and no prefix; otherwise only `parm7` and
  `rst7` are written.
- `<prefix>` defaults to the input name; choose another so that `<prefix>.pdb`
  does not replace the input.
- To check the topology:

  ```bash
  parmed -p system.parm7 -i <(echo "summary"; echo "go")
  ```

## Pitfalls and recovery

- AmberTools must be on `PATH`: [ambertools.md](../mlmm-install-backends/ambertools.md).
- `-l` accepts both `=` and `:` (`'GPP=-3,SAM:1'`). Use `:` to match the rest
  of the toolkit.
- Use `--add-h` only when the PDB lacks hydrogens at the protonation you want;
  otherwise the PDB goes to tleap as it is.
- The `parm7` holds no layers. Downstream commands read the ML, Movable-MM,
  and Frozen-MM layers from the PDB B-factors written by `define-layer`.
- Charge and hydrogens must agree. AM1-BCC runs `sqm`, which needs a
  closed-shell electron count at multiplicity 1. `mm-parm` and `all` check
  this before antechamber and stop with
  `[<RES>] electron-count check failed before antechamber`, naming the residue.
  Without that check, `sqm` aborts with `The number of electrons is odd`.
- Count the hydrogens before a run (Σ Z − q must be even for closed shell):

  ```bash
  awk '$4=="<RES>" && substr($0,77,2)~/H/' input.pdb | wc -l
  ```

- SAM has 22 H at charge 0 (NH2/COO⁻/S⁺ zwitterion) and 23 H at +1
  (NH3⁺/COO⁻/S⁺, the usual biological form). A −3 diphosphate such as DMAPP
  has 9 H. ATP / ADP: 12–13 H for the usual −4 / −3, depending on the
  dataset; count the H in your file. Do not add or remove
  protons to fix the parity; keep the dataset's protonation and change `-l`:

  ```bash
  mlmm mm-parm -i in.pdb -l 'SAM:0,...' --out-prefix system   # 22 H
  mlmm mm-parm -i in.pdb -l 'SAM:1,...' --out-prefix system   # 23 H
  ```

## Next step

- [define-layer.md](define-layer.md): assign layers on `system.pdb`.
- [extract.md](extract.md): cut the ML region from `system.pdb`.
- [Amber parm7 and rst7](../mlmm-structure-io/formats.md#amber-parm7-and-rst7):
  what the two files contain.
