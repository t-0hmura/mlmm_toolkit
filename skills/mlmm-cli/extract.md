# `mlmm extract`

Cuts a binding pocket around the `-c` residues out of a protein-substrate
complex, with residue-aware truncation and optional carbon-only link H. Run
`mlmm extract -i complex.pdb -c 'A:SAM:44' -l 'SAM:1' -o pocket.pdb`. Success
is `pocket.pdb` and a `Total active site model charge` line with the charge
you expect.

## When to use

- Build the active-site model, or check what `all -c` would put in the ML
  region, before a long run.
- Several inputs with the same atoms give one pocket per file (or one
  multi-MODEL PDB) with the same atoms, for `path-search`-style runs.
- `extract` writes no B-factor layers; assign them with
  [define-layer.md](define-layer.md).

## Minimal run

```bash
mlmm extract -i 1abc.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3' -r 4.0 -o pocket.pdb
```

`-c` takes a PDB path, residue IDs (`'A:44,B:321'`), or names (`'GPP,SAM'`);
with chain IDs, `A:SAM:44` is the safest form. `-l` takes the total charge or
a per-residue mapping. With link H and `result.json`:

```bash
mlmm extract -i complex.pdb -c 'A:123,A:124' -l 'GPP:-3' -r 3.5 \
    --add-linkh --out-json -o pocket.pdb
```

Several inputs (same atom count and order):

```bash
mlmm extract -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -o 1.R_pocket.pdb 3.P_pocket.pdb
```

## Judge success

- The console prints `[extract] Atoms after truncation: N` and
  `Total active site model charge: …`. N excludes the cap H. With several
  inputs, the line is `[extract:multi] Atoms after truncation (model k): N`
  for each model.
- `pocket.pdb` exists (default name; `pocket_<filename>.pdb` for several
  inputs without `-o`).
- With `--out-json`, `result.json` has `total_charge`, `n_atoms_extracted`
  (N, without cap H), and `n_link_hydrogens` (M). The model computes N + M
  atoms.

```python
import json
d = json.load(open("result.json"))
print(d["total_charge"], d["n_atoms_extracted"], d["n_link_hydrogens"])
```

## Pitfalls and recovery

- Defaults: `--add-linkh` off, `--include-h2o` on (`--no-include-h2o` for a
  dry pocket), `--exclude-backbone` off (on for cluster-style truncated
  backbones).
- For a reusable standalone ML/MM workflow, follow `mm-parm → extract →
  define-layer → opt/tsopt/...`. Request that PDB with a distinct prefix, for
  example `mlmm mm-parm -i input.pdb --out-prefix system`, then run the latter
  two commands on `system.pdb`. LEaP may change hydrogens; the exported PDB has
  the same atom identity/order as `system.parm7`, and `mm-parm` fills missing
  element columns.
- Atom names match exactly (case-sensitive). Run `mlmm add-elem-info` first
  on a PDB from PyMOL or Maestro. `extract` keeps one altLoc per residue
  itself; run `fix-altloc` only when you need the cleaned file.
- `--add-linkh` is for standalone pockets, not for an mlmm `--model-pdb`. The
  calculator caps the ML/MM boundary with link H from the `--parm7`
  topology; extract's link H is distance-based and can misfire on unusual
  topologies.
- Automatic extraction versus a hand-built ML region, and which charge
  options apply to each: [Build the ML region](../mlmm-model-setup/SKILL.md#build-the-ml-region).

## Next step

- [define-layer.md](define-layer.md): assign ML, Movable-MM, and Frozen-MM.
- [mm-parm.md](mm-parm.md): the topology and the matching PDB to extract from.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md): grow or trim the ML region.
- [Selecting residues and atoms](../mlmm-model-setup/SKILL.md#selecting-residues-and-atoms)
  and [PDB](../mlmm-model-setup/formats.md#pdb).
- [add-elem-info](utilities.md#add-elem-info) and
  [fix-altloc](utilities.md#fix-altloc): pre-clean a raw PDB.
