---
name: mlmm-structure-io
description: "Input structures for mlmm-toolkit: PDB, mmCIF, XYZ, Gaussian gjf, and Amber parm7/rst7, with residue selectors, ML-region and layer encoding (B-factor 0/10/20 and `model.pdb`), atom-order checks between the PDB and parm7, and the charge/multiplicity decision. `formats.md` holds per-format details. TRIGGER on inspecting or editing a structure, choosing `-q` / `-l` / `-m`, reading or checking B-factor layers, or a coordinate/topology mismatch. SKIP for choosing which atoms go into the ML region or layers (mlmm-model-setup), subcommand syntax, output parsing, install, or HPC questions."
---

# mlmm-toolkit structure input

Give `-i` the whole system as PDB (or XYZ with `--ref-pdb`), `--parm7` from `mm-parm`, the ML region by B-factor or `--model-pdb`, and the ML-region charge with `-q` or `-l`.

```bash
mlmm opt -i complex_layered.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' -m 1 -o result_opt
```

A structure is ready when a calculation starts without an atom-count or
`Atom-order mismatch` error and the ML-region charge, given with `-q` or
derived from `-l`, matches your own count. Byte-level layouts are in
[formats.md](formats.md).

## Which format

| Format | Carries | Use for |
|---|---|---|
| PDB | Atom and residue names, chain, occupancy, B-factor (layer), element | The normal input; the B-factor holds the ML, Movable-MM, and Frozen-MM layers |
| mmCIF | The same, without the one-character chain and four-digit residue limits | Long chain IDs, residue numbers of 10,000 or more, oversized structures |
| XYZ | Element and Cartesian coordinates only | Trajectories, single TS candidates, exchange between subcommands |
| GJF | Coordinates with charge, spin, and route line | Gaussian ONIOM exchange through `oniom-export` and `oniom-import` |
| parm7 / rst7 | Amber topology and coordinates | MM parameters of the full system, written by `mm-parm` |

PDB, mmCIF, XYZ, and GJF use Å and ordinary element symbols.

```text
Full enzyme you will run ML/MM on?
  └── PDB or mmCIF with B-factor layers + parm7
      → opt / tsopt / scans / path commands / freq / irc / dft / all

A single TS candidate to validate?
  └── XYZ + --ref-pdb (full-system PDB/mmCIF) + --parm7
      → tsopt / freq / irc / all (TS-only mode)

A Gaussian or ORCA ONIOM input to bring in?
  └── GJF/INP → mlmm oniom-import → XYZ + layer-encoded PDB

A raw enzyme PDB that needs a parm7?
  └── mlmm mm-parm → parm7 + rst7 (AmberTools tleap)
```

## Which subcommand reads which format

| Subcommand | PDB/mmCIF | XYZ | GJF | parm7 |
|---|---|---|---|---|
| `extract` | ✓ in/out | — | — | — |
| `mm-parm` | PDB in | — | — | ✓ out |
| `define-layer` | ✓ in/out | — | — | — |
| `path-search` / `path-opt` | ✓ | ✓ with `--ref-pdb` | — | required |
| `sp` / `opt` / `tsopt` / `freq` / `irc` / `dft` | ✓ | ✓ with `--ref-pdb` | — | required |
| `scan` / `scan2d` / `scan3d` | ✓ | ✓ with `--ref-pdb` | — | required |
| `oniom-export` | PDB in | — | ✓ out | required |
| `oniom-import` | PDB out | XYZ out | ✓ in | — |

For an mmCIF input, `extract` and `define-layer` also write `.cif` files that
restore the original chain IDs and residue numbers.

## Selecting residues and atoms

`-c/--center` on `extract` and `all` names the residues at the center of the
model. Write the chain first:

```bash
mlmm extract -i complex.pdb -c 'A:SAM:44' -o cluster.pdb     # chain + name + number: one residue
mlmm extract -i complex.pdb -c 'A:44' -o cluster.pdb         # chain + number (a trailing letter is the insertion code)
mlmm extract -i complex.pdb -c 'SAM,GPP,MG' -o cluster.pdb   # names: every residue with that name
mlmm extract -i complex.pdb -c substrate.pdb -o cluster.pdb  # residues matching a separate PDB
```

`A:SAM` selects every SAM in chain A and logs a warning when more than one
matches; add the number (`A:SAM:44`) when the match must be unique. Names or
numbers without a chain can match several residues, so use the exact forms
for production runs. A PDB with an empty chain column, such as the bundled
examples, takes only the name or number forms (`-c 'SAM,GPP,MG'`).

Long chain IDs and residue numbers of 10,000 or more use the same forms through
mmCIF, for example `enzyme_A:SAM:10001B`
([mmCIF](formats.md#mmcif-and-very-large-structures)).

Single atoms, as in scan lists, take four fields
`CHAIN:RESNAME:RESSEQ[ICODE]:ATOM` (`A:SAM:320:CS1`). On a PDB with an empty
chain column, use three fields without the chain (`SAM,320,CS1`).

## ML region and layers

The B-factor column carries the layers:

| B-factor | Layer |
|---|---|
| 0 | ML |
| 10 | Movable-MM |
| 20 | Frozen-MM |

Values within ±1.0 count. Each command takes the ML region from the first of:

1. `--model-pdb FILE`
2. `--model-indices '1-50,75,100-110'`
3. the B-factor-0 atoms of the input PDB under `--detect-layer` (on by default)

With an explicit ML region, `--detect-layer` still reads the Movable-MM and
Frozen-MM layers from the B-factors. The B-factors count as layers only when at
least one atom is ML, at least one is MM, and at least 80% of the atoms carry
0, 10, or 20; an all-zero PDB is not a layer assignment. `--movable-cutoff`
turns off `--detect-layer` and sets the MM layers by distance instead.

The parm7 carries only MM parameters, never the layers. A command reads the
layered structure from `-i` and the topology from `--parm7`:

```bash
mlmm opt -i complex.pdb --parm7 complex.parm7 -q 0 -m 1 -b uma -o result_opt
```

When `-i` is XYZ, also pass the full-system PDB or mmCIF to `--ref-pdb`; it
supplies the atom order and residue context and nothing else. To change
layers, run `mlmm define-layer` ([define-layer](../mlmm-cli/define-layer.md));
a one-residue edit is in [formats.md](formats.md#amber-parm7-and-rst7). Which
atoms belong in the ML region and in each layer is in
[mlmm-model-setup](../mlmm-model-setup/SKILL.md).

## Charge and multiplicity

`-q` is the charge of the ML region, not of the whole system. A wrong charge or
multiplicity can silently give a chemically wrong trajectory.

For PDB/mmCIF input, give `-l 'RES:Q'` for unknown or non-standard residues
only and let mlmm sum the ML-region total. Standard amino acids and recognized
ions come from internal tables; waters and link atoms are neutral. Recheck the
reported breakdown whenever the ML region, residue naming, protonation state,
or oxidation state changes.

The charge is taken from the first of:

1. `-q/--charge`
2. with `-l/--ligand-charge`, the sum of standard residues, ions, and your ligand charges in the ML region
3. `calc.model_charge` from `--config`
4. otherwise, the run stops with an error

Use `-q` when the input has no residue metadata or to override the derived
value; in `mlmm all`, `-q` wins and the workflow reports the value it would
have derived. With `--model-indices`, `-l` cannot derive the charge: give `-q`,
or define the ML region with `--model-pdb` or B-factor layers. XYZ input gets
its residue context from `--ref-pdb` and follows the same rules.

Recognized monatomic ions keep their table value. Listing one in `-l` with the
same value is accepted (`MG:2`); a different value (`MG:3`) is ignored with a
warning. For another oxidation state, use the matching residue name in the
model or give the verified total with `-q`. To see the tables:

```bash
python -c "from mlmm.core.residue_data import AMINO_ACIDS, ION; print(dict(AMINO_ACIDS)); print(dict(ION))"
```

Without `-m`, mlmm uses `calc.model_mult` and then 1. The default 1 is an input
default, not a scientific assignment: use it only when the modeled electron
count and state are known to be closed-shell, not merely because the system is
biological or metal-bound. Determine the multiplicity from composition,
oxidation states, experiment, or explicit state comparisons; metals and
radicals need particular care.

| `-m` | Examples |
|---|---|
| 1 | Closed shell |
| 2 | Radicals, unpaired-electron TSs (radical SAM enzymes, low-spin Fe(III)) |
| 3 | O₂, some carbenes, high-spin Ni(II) (d⁸) in tetrahedral or weak-field octahedral sites |
| 4 | High-spin Co²⁺ (d⁷), Cr³⁺ / V²⁺ (d³) |
| 5 | Mn(III), high-spin Fe(II) |
| 6 | High-spin Mn(II), S=5/2 ferric |

These are examples, not a spin-state calculator. For metals, radicals,
antiferromagnetically coupled centers, or uncertain protonation or oxidation
states, derive charge and multiplicity from the modeled mechanism and primary
literature. Common ligand, ion, and metal values are in
[formats.md](formats.md#ligand-ion-and-metal-charges).

## Unknown substrate charge

When a ligand's formal charge is unknown:

1. Check the primary paper. Most mechanism papers state the substrate charge
   state in the Methods; the PDB entry page links to the reference.
2. Look it up. [PubChem](https://pubchem.ncbi.nlm.nih.gov) lists `Formal Charge`
   under Computed Properties (search by three-letter code or name);
   [ChEBI](https://www.ebi.ac.uk/chebi) often has the state used in published
   mechanisms; the RCSB ligand page (`https://www.rcsb.org/ligand/SAM`) shows
   the SMILES and charge of the deposited model.
3. Derive it from the SMILES:

   ```python
   from rdkit import Chem
   mol = Chem.MolFromSmiles("CC(=O)[O-]")     # acetate
   print(sum(a.GetFormalCharge() for a in mol.GetAtoms()))    # → -1
   ```

4. Check the protonation state at pH 7. Typical contributions:

   | Group | At pH 7 | Charge |
   |---|---|---|
   | Carboxylate | deprotonated | −1 each |
   | Phosphate monoester | mostly `-OPO₃²⁻` | −2 |
   | Phosphate diester | mostly `-OPO₂⁻` | −1 |
   | Triphosphate (ATP) | fully deprotonated | −4 |
   | Sulfonium (SAM) | quaternary | +1 |
   | Lys / Arg side chain | protonated | +1 |
   | Asp / Glu side chain | deprotonated | −1 |
   | His | mostly neutral | 0 or +1 |

   Mechanisms sometimes invoke an unusual protonation state; check the
   literature for the model you build.

5. Check the total that mlmm reads. `extract --out-json` writes `result.json`
   next to the output PDB with `total_charge` and its breakdown
   (`protein_charge`, `ligand_total_charge`, `ion_total_charge`); the terminal
   prints `Total active site model charge`.

   ```bash
   mlmm extract -i complex.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' -o cluster.pdb --out-json
   python -c "import json; print(json.load(open('result.json'))['total_charge'])"
   ```

   For a layered full system, a one-cycle optimization records the charge it
   used:

   ```bash
   mlmm opt -i complex_layered.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' -m 1 --max-cycles 1 -o check_opt --out-json
   python -c "import json; print(json.load(open('check_opt/result.json'))['charge'])"
   ```

Use sources in this order: the mechanism's primary paper or deposited structure
documentation, then PubChem, ChEBI, or the RCSB CCD. Cite the source and state
the modeled protonation and oxidation state. If the sources do not settle one
state, ask rather than defaulting a metal or radical model to `-q 0 -m 1`.

## Editing approach

When you edit a structure file:

1. Read it first: residues, atom counts, B-factor layers, and any charge or
   multiplicity.
2. Confirm the change keeps the format: PDB column widths, the XYZ atom-count
   line, the parm7 layout.
3. For an unknown charge or multiplicity, confirm with the user or follow
   [Unknown substrate charge](#unknown-substrate-charge) before guessing.
4. For layer changes, prefer `mlmm define-layer` over editing B-factors by hand.

## Next step

- [extract](../mlmm-cli/extract.md), [mm-parm](../mlmm-cli/mm-parm.md), and
  [define-layer](../mlmm-cli/define-layer.md): the preparation commands.
- [mlmm-model-setup](../mlmm-model-setup/SKILL.md): choosing the ML region and layers.
- [Outputs](../mlmm-overview/outputs.md): the XYZ, PDB, and CIF files a run writes.
- [formats.md](formats.md): per-format layouts, edits, and checks.
