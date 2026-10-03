# Building the ML region and layers

mlmm-toolkit computes the whole enzyme: the ML region around the substrate with an MLIP, and the rest of the protein with the Amber force field (MM). The cost of a run is set by three things: the number of atoms in the ML region, the number of MM atoms that move (Movable-MM), and the range of atoms in the Hessian. This page shows how to build the ML region and the layers, how to make the model smaller so the calculation is lighter, how to make it larger when residues are missing, and how to freeze atoms and restrain distances.

## Quick guide

| Goal | What to do | Section |
| --- | --- | --- |
| Build the default model | Give the substrate, cofactors, and metals to `-c` of `all` | [Build the default model](#build-the-default-model) |
| See the atom counts and the charge before a long run | Read the Layer Summary of `define-layer` or the log of `all` | [Check the model](#check-the-model) |
| Make the calculation lighter | Shrink the ML region, thin the movable MM shell, or use an ML-only Hessian | [Make the model smaller](#make-the-model-smaller) |
| A residue of the reaction is missing | Raise `-r`, or add the residue to `-c` or `--selected-resn` | [Make the model larger](#make-the-model-larger) |
| Use a model you built yourself | Pass `--parm7` and `--model-pdb` | [Use a model you built yourself](#use-a-model-you-built-yourself) |
| Write `model.pdb` by hand | Follow the checklist | {ref}`How to construct a reliable model.pdb <model-pdb-selection>` |
| Freeze atoms or restrain a distance | The Frozen-MM layer, `--freeze-atoms`, `--distance-restraint` | {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>` |

## Build the default model

### The ML region

Give the substrate to `-c` of `all`, and put the cofactors and metals in the same list.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'
```

This is the bundled example in [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza). `all` cuts the ML region out with `extract`: a residue joins when any of its atoms lies within `-r` (default 2.6 Å) of an atom in `-c`, and waters are included. The full rules are in [extract](extract.md#how-it-works). `all` writes the ML region of the first input to `ml_region.pdb` in the output directory; pass it with `--model-pdb` to use the same ML region in a later run.

(mm-layers)=
### The MM layers

The atoms outside the ML region form two MM layers. The layer of each atom is stored in the B-factor column of the PDB.

| Layer | B-factor | What it does |
| --- | --- | --- |
| **ML** | 0 | The reactive region; energy, forces, and Hessian from the MLIP |
| **Movable-MM** | 10 | MM atoms that move during optimizations |
| **Frozen-MM** | 20 | MM atoms whose coordinates stay fixed; they still take part in the MM energy |

[`define-layer`](define-layer.md) puts the MM atoms within `--movable-cutoff` (default 8.0 Å) of the ML atoms in Movable-MM and the rest in Frozen-MM. With `-c`, `all` runs `define-layer` at 8.0 Å on each input and writes the layered PDBs to `layered/` in the output directory.

## Check the model

Check the atom counts and the charge before a long run.

- **Layers**: `define-layer` prints a Layer Summary with the lines `Layer 1 (ML, B=0):`, `Layer 2 (Movable MM, B=10):`, `Layer 3 (Frozen MM, B=20):`, and `Total atoms:`. `all` prints `[all] define-layer [i]: … (ML=…, MovableMM=…, FrozenMM=…)` for each input.
- **Link hydrogens**: `all` prints `[all] ML structure with link H (N + M; …)`, where N is the number of ML atoms and M the number of link hydrogens.
- **Charge**: `extract`, and `all` with `-c`, print the charge of the ML region on the line `Total active site model charge`. Check the multiplicity too.
- **By eye**: open the layered PDB in a viewer, color it by B-factor, and check that the residues of the reaction are in the ML region.

## Make the model smaller

Each of the three costs has its own setting.

| What it sets | Why it costs | In `all` | In the individual commands | Default |
| --- | --- | --- | --- | --- |
| ML region | Computed with the MLIP (with DFT under `-b dft`) at every step | `-r`, `--exclude-backbone`, `--no-include-h2o`, `--selected-resn` | `--model-pdb` | `-r 2.6` Å |
| Movable-MM | The atoms that an optimization moves | No option (8.0 Å with `-c`) | `define-layer --movable-cutoff`, or `--movable-cutoff` of `opt`, `tsopt`, `freq`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, `sp` | 8.0 Å |
| Hessian range | A dense matrix over the ML region and Movable-MM | No option | `--hessian-cutoff` of `opt`, `tsopt`, `freq`, `sp` | ML and all of Movable-MM |

### Shrink the ML region

For the bundled example (`examples/beza/1.R.pdb` with `-c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3'`), these options give ML regions of the following sizes. The counts were taken with `extract --add-linkh` and include the cap hydrogens.

| How to build the ML region | Options | Atoms |
| --- | --- | --- |
| Default | — | 632 |
| Move main-chain atoms to MM | `--exclude-backbone` | 397 |
| …and waters | `--exclude-backbone --no-include-h2o` | 367 |
| Pick residues yourself | `-r 0 --selected-resn '44,63,186'` | 129 |

- **`--exclude-backbone`** and **`--no-include-h2o`** move main-chain atoms and waters out of the ML region; they stay in the system as MM.
- **`-r 0 --selected-resn`** adds no residues by distance: the ML region holds the `-c` residues and the residues you list. In the bundled example, residues 44, 63, and 186 are the three closest to the methyl carbon of SAM (CS1); for your own system, pick the residues that take part in the reaction.

`all` takes the same options; see example 1 of [Refine an MLIP TS with DFT](dft-backend.md#examples), and [Keep the ML region under about 300 atoms](dft-backend.md#keep-the-ml-region-under-about-300-atoms) for a DFT/MM run. A smaller ML region often has a different charge, so check it each time.

### Thin the movable MM shell

A smaller `--movable-cutoff` puts fewer atoms in Movable-MM and more in Frozen-MM, so optimizations move fewer atoms. Assign the layers again with `define-layer --movable-cutoff` and pass the new PDB to the individual commands, or give `--movable-cutoff` to one of the commands in the table above.

### Narrow the Hessian

`--hessian-cutoff` keeps in the Hessian only the Movable-MM atoms within that distance of the ML region. The vibrational analysis of `freq`, and the one at the end of `tsopt`, needs a Hessian that covers its atoms (`--active-dof-mode`; the default `partial` is the ML region and all of Movable-MM). When the Hessian is narrower, the run stops with an error before the calculation, so narrow the analysis together with it: for an ML-only Hessian, pass `--hessian-cutoff 0.0 --active-dof-mode ml-only`. `all` has no `--hessian-cutoff`. How the MM part of the Hessian is computed is described in [MM Hessian](mlmm-calc.md#mm-hessian).

With many movable MM atoms, {ref}`microiteration <microiteration>` (on by default in `tsopt` and in `opt --opt-mode hess`) relaxes the MM atoms with the force field alone between the ML steps, so the MLIP is called fewer times.

## Make the model larger

Make the ML region larger when a residue, water, or cofactor of the reaction lies outside it. The rest of the protein is already in MM, so add to the ML region the atoms whose bonds, protonation, or charge change, together with their covalent partners.

- **Raise `-r`** (default 2.6 Å).
- **`--selected-resn`** adds residues without starting a distance search from them.
- **`--radius-het2het`** (default 0, off) adds a second cutoff measured only between atoms other than C and H, on both the center side and the neighbor side. It picks up close N and O partners without enlarging the whole radius.
- **Add the residue to `-c`** to keep it whole: amino acids in `-c` start their own distance search and, without `--exclude-backbone`, keep all their atoms.

To let more of the environment relax, widen Movable-MM instead of the ML region, for example with `define-layer --movable-cutoff 10.0`. If you narrowed the Hessian to the ML region, drop `--hessian-cutoff` to go back to the default; see [Tips for studying reaction mechanisms](mechanism-tips.md).

The radius is a parameter for checking that the result has converged with the size of the ML region. A larger ML region costs more and does not always improve accuracy, so compare energies, forces, and barriers over a few chemically sensible ML regions. Use the same ML region for the reactant, intermediates, and product: pass `ml_region.pdb` with `--model-pdb`.

## Use a model you built yourself

The individual commands need `--parm7` and an ML region, and `all` builds both. To build the model by hand, run three commands: `mm-parm` builds the topology and a PDB with the same atoms in the same order, `extract` cuts the ML region out of that PDB, and `define-layer` writes the layers into it.

```bash
mlmm mm-parm -i input.pdb -l 'LIG:0' --out-prefix system
mlmm extract -i system.pdb -c LIG -l 'LIG:0' -o model.pdb
mlmm define-layer -i system.pdb --model-pdb model.pdb -o system_layered.pdb
```

Pass `system_layered.pdb` with `--parm7 system.parm7` and `--model-pdb model.pdb` to the individual commands, or to `all` without `-c`. Each command takes the ML region from `--model-pdb`, then `--model-indices`, then the B-factor layers ({ref}`ML/MM options <mlmm-options>`). Give the charge of the ML region with `-q`, or the charge of each residue name with `-l` (PDB/mmCIF input). When you cut the region or changed protonation by hand, the charge cannot be derived from the residue names; give it with `-q`.

(model-pdb-selection)=
### How to construct a reliable `model.pdb`

`model.pdb` is an **atom-selection file**, not an independently rebuilt cluster.
Every atom must be an unchanged subset of the full PDB/`parm7` topology: preserve
atom names, residue names/numbers, chain IDs, and full-system atom order. Do not
renumber, reorder, add link hydrogens, or export a separately hydrogenated model.

- Include the complete reactive center, covalent cofactors/partners, and any
  atoms whose protonation or bonding changes along the path.
- For retained protein-backbone fragments, choose the span so both main-chain
  ends terminate consistently at alpha carbons (`CA`), then let the ML/MM link
  treatment satisfy boundary valences.
- At side-chain/ligand/cofactor boundaries, place the ML/MM cut on an aliphatic
  **C–C single bond** whenever possible (`CA–CB` or farther from the reactive
  center). Avoid peptide C–N, polar C–N/C–O, aromatic/conjugated, disulfide,
  and metal-coordination cuts; include the bonded partner or move the boundary.
- Use the identical full-system atom set/order and the identical `model.pdb`
  selection for R/IM/P. A model built independently for each state invalidates
  atom mapping and controlled barrier comparisons.
- Visually inspect every boundary and verify the ML-region charge/multiplicity
  before production. `define-layer` assigns layers; it does not repair a
  chemically poor boundary.

When you save `model.pdb` from PyMOL, tick **Original atom order** in the export dialog.

(freeze-atoms-and-restraints)=
## Freeze atoms and restrain distances

Freezing keeps an atom in place by setting its force to zero. A restraint pulls the distance between two atoms toward a target with a harmonic potential. In ML/MM the Frozen-MM layer already holds the outer protein in place, so give `--freeze-atoms` only the extra atoms you want to hold.

### Three ways to freeze atoms

- **The Frozen-MM layer** (B-factor 20) is frozen automatically; see {ref}`The MM layers <mm-layers>`.
- **`--freeze-atoms 'i,j,k'`** takes 1-based atom numbers of the full system. It works in `all`, `opt`, `tsopt`, `irc`, `freq`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `sp`.
- **YAML `geom.freeze_atoms`** (passed with `--config`) suits a long list, or a list you keep with the rest of the settings.

```yaml
geom:
  freeze_atoms: [12, 15, 28, 29, 42]   # 1-based
```

A run freezes the union of the three; none of them replaces another.

### What freezing does

- Frozen atoms get zero force, so they do not move.
- Frozen atoms are left out of the Hessian, so `freq` runs a partial Hessian vibrational analysis (PHVA) on the other atoms. For the rigid motions that are removed, see [Rigid modes with frozen boundaries](freq.md#rigid-modes-with-frozen-boundaries).
- `--mep-mode dmf` (Direct Max Flux) in `path-opt` and `path-search` holds frozen atoms with a harmonic restraint instead (k = 300 eV/Å², YAML `dmf.k_fix`), so they can move slightly; see the [`path-opt` notes](path-opt.md#notes).

### Restrain a distance

`--distance-restraint` of `opt` adds a harmonic restraint between two atoms. Give `(i, j, target)` with the target in Å, or `(i, j)` to keep the starting distance.

```bash
mlmm opt -i system_layered.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
    -q 0 -m 1 --distance-restraint '[(12,45,2.20)]' --restraint-k 20.0 \
    --out-dir ./result_opt_rest
```

- Atom numbers are 1-based; `--zero-based` counts from 0. `--restraint-k` sets the force constant (default 300 eV/Å²).
- In `scan`, `scan2d`, and `scan3d`, `--restraint-k` is the force constant of each step (eV/Å² for distances, eV/rad² for angles) and takes priority over YAML `bias.k`. Each stage restrains only its own coordinates. In `all`, pass the value with `--scan-restraint-k`. For the format of the scan lists, see {ref}`Scan-list spec <scan-list-spec>`.
- Reported energies leave out the energy of the restraint.

## Notes

- **`--movable-cutoff` replaces the B-factor MM layers** in every command that takes it. In `opt`, `tsopt`, `freq`, `path-opt`, and `path-search` it also turns off `--detect-layer`, so give the ML region with `--model-pdb` or `--model-indices`.
- **`--model-pdb` sets only the ML region**: with `--detect-layer`, Movable-MM and Frozen-MM still come from the B-factors of the input.
- **Change the layers with `define-layer`**: run it again instead of editing B-factors by hand.
- **The bundled PDB has an empty chain column**; in such a PDB, give residues by name or by number. In a PDB with chains, write each residue as chain:name:number, such as `-c 'A:TYR:44'`; see [Residue selectors](cli-conventions.md#residue-selectors).
- **Link hydrogens are added for you**: the calculator puts one on each `parm7` bond that has exactly one end in the ML region, so `model.pdb` does not contain them. See [Link-atom redistribution](mlmm-calc.md#link-atom-redistribution).
- **Boundary bonds**: link hydrogens go only on C–C, C–N, and N–C bonds. Any other `parm7` bond across the boundary stops the run with `Unsupported ML/MM boundary bond in parm7`; move the boundary to a supported bond.
- **Metals, glycans, and MD snapshots**: build the `parm7` yourself and pass it with `--parm7`; see the [`mm-parm` notes](mm-parm.md#notes).

## See also

- [`extract`](extract.md) — extraction options, the residue selectors of `-c`, and non-standard residue names
- [`define-layer`](define-layer.md) — assign the ML, Movable-MM, and Frozen-MM layers
- [`mm-parm`](mm-parm.md) — build the Amber topology and the matching PDB
- [`all`](all.md) — the full workflow; builds the ML region and the layers with `-c`
- [`opt`](opt.md) — optimization with distance restraints
- [`scan`](scan.md) — staged scans with restraints
- [`freq`](freq.md) — PHVA and rigid modes with frozen atoms
- [Tips for studying reaction mechanisms](mechanism-tips.md) — when to enlarge the model
- [Refine an MLIP TS with DFT](dft-backend.md) — keeping the ML region small enough for DFT
- [ML/MM Calculator](mlmm-calc.md) — link atoms, microiteration, and the MM Hessian
- {ref}`ML/MM options <mlmm-options>` — `--parm7`, `--model-pdb`, `--detect-layer`, and `--movable-cutoff`
- [Device Configuration & HPC Setup](device-hpc.md) — GPU memory and the size of the model
- [Troubleshooting](troubleshooting.md) — extraction and layer errors
