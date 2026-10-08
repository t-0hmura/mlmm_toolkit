# `path-search` (recursive MEP through two or more structures)

`path-search` builds one continuous minimum-energy path (MEP) through **two or more** layered enzyme structures given in reaction order (R → … → P), using the ML/MM calculator on the whole system. It refines the path recursively, only in the regions where covalent bonds change, and builds each piece with GSM (growing string method, the default) or DMF (direct max flux).

## What it is for

* **Splitting R → P into reactive segments**: when you do not know whether the reaction has one step or several, find the regions where bonds change.
* **A multistep path through intermediates**: give known intermediates between R and P and get one stitched path.
* **TS candidates per segment**: each reactive segment gets its own HEI (highest-energy image), `hei_seg_NN.xyz`, to optimize with [`tsopt`](tsopt.md).

The ML region is computed with **UMA** (Meta) by default, and the rest of the system is computed with the Amber force field of `--parm7`. For exactly two endpoints without recursive refinement, [`path-opt`](path-opt.md) is simpler.

---

## Examples

### 1. Two endpoints

Give the reactant and the product after one `-i`, with the charge of the ML region and the spin multiplicity. `reactant.pdb` and `product.pdb` hold the whole system that matches `real.parm7` (the Amber topology), and `ml_region.pdb` selects its ML region.

```bash
mlmm path-search -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-dir ./result_path_search
```

When the run finishes, open `summary.log` (section `[2] Segment-level MEP summary`), or read `summary.json`. `scientific_status` is `success` when the pre-optimizations and every path run converged, otherwise `partial` or `failed`. `segments` lists each segment with its `index`, `tag`, `kind`, `bond_changes`, `converged`, and `barrier_kcal`. Each reactive segment also has its TS candidate, `hei_seg_NN.xyz`.

### 2. Add intermediates for a multistep path

List the structures in reaction order after one `-i`; each adjacent pair is searched and the pieces are stitched into one path.

```bash
mlmm path-search -i R.pdb IM1.pdb IM2.pdb P.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q -1 -m 1 --out-dir ./result_path_search_multi
```

### 3. Lighter pass without pre-optimization or alignment

Skip the pre-optimization and the alignment of the inputs and use fewer movable images, for inputs that are already optimized and superimposed.

```bash
mlmm path-search -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --no-preopt --no-align --max-nodes 8 --out-dir ./result_path_search_fast
```

---

## How it works

Before the search, each input is pre-optimized (`--preopt`) and aligned to the one before it (`--align`), both by default, with frozen atoms matched step by step while the other atoms relax.

1. **A coarse MEP for each pair**:
Between each pair of adjacent inputs (A → B), GSM or DMF builds a coarse MEP and finds its HEI.
2. **Relaxing around the HEI**:
`--refine-mode peak` optimizes the images on either side of the HEI (HEI ± 1); `minima` searches outward from the HEI for the nearest local minimum on each side. The result is two nearby minima, End1 and End2. When `--refine-mode` is omitted, GSM uses `peak` and DMF uses `minima`.
3. **Kink or reactive segment**:
If no covalent bond changes between End1 and End2, the region is a *kink*: `path-search` inserts a few linear nodes and optimizes each one. Otherwise the region is a *reactive segment*, and a new GSM or DMF path between End1 and End2 sharpens its barrier.
4. **Recursing where bonds still change**:
The parts A → End1 and End2 → B are checked for bond changes, and only parts that still have them are searched again, down to `--max-depth` levels.
5. **Stitching**:
The pieces are joined into one path. Duplicate endpoints are dropped; where the ends of two neighboring pieces still differ in bonds, that gap is searched as a new segment, and any other gap is filled with a short connecting path.

Bond changes are judged with the thresholds in the YAML `bond` section, by the same rules as in {ref}`scan <section-bond>`.

---

## Reading the segments

| What you see | Meaning | Next step |
| --- | --- | --- |
| A segment with bond changes, with its `hei_seg_NN.xyz` | A TS candidate for that step | Optimize it with [`tsopt`](tsopt.md), check for one imaginary mode, then run [`irc`](irc.md) |
| A segment whose `tag` is `seg_NNN_maxdepth` or `seg_NNN_kinklimit` | Splitting stopped there, at the depth limit (`_maxdepth`) or after consecutive kinks (`_kinklimit`) | It may hold more than one step; check it as above, raise `--max-depth`, or give intermediates |
| Only segments whose `tag` ends in `_kink`, or the warning `HEI is at an endpoint` | No bond change was found, or the path has no peak between its ends | Check the inputs, or give intermediates (example 2) |

The segmentation is a guide based on bond-distance criteria. One segment is not guaranteed to be one elementary step or to contain exactly one TS. A successful TS optimization gives one imaginary mode along the reaction coordinate. Confirm every HEI with `tsopt` (n_imag = 1) and IRC before you read it as a step of the mechanism.

---

## Output files

`path-search` writes these files to `--out-dir` (default `./result_path_search/`):

```text
result_path_search/
├─ mep_trj.xyz               # The whole stitched MEP, energies on the comment lines
├─ mep_trj.pdb               # Same path as PDB
├─ mep_plot.png              # ΔE profile along the path (kcal/mol, relative to the reactant)
├─ energy_diagram_MEP.png    # State-energy diagram of the MEP (relative to the reactant)
├─ summary.json              # Barrier and classification summary for every segment
├─ summary.log               # The same summary as text
├─ mep_seg_NN_trj.xyz        # Path of reactive segment NN (PDB: mep_seg_NN.pdb)
├─ hei_seg_NN.xyz            # HEI of reactive segment NN, the TS candidate (PDB: hei_seg_NN.pdb)
├─ hei_mode_seg_NN.*         # Reaction-direction guess at that HEI; all passes it to tsopt
├─ align_refine/             # Alignment and relaxation files of the inputs (--align)
├─ initNN_*_opt/             # Pre-optimization of each input (--preopt)
└─ seg_NNN_*/                # Working files of each GSM/DMF run and HEI-side optimization
```

`summary.json` is always written and has its own structure, unlike the `result.json` of the other commands; see the section `summary.json (path-search / all)` of the [JSON Output Reference](json-output.md). Only segments with bond changes get `mep_seg_NN_*` and `hei_seg_NN.*` files. NN is the segment's `index` in `summary.json` (counted from 01 along the final path), while NNN in a `seg_NNN` tag or directory counts the GSM/DMF runs from 000, so the two numbers differ. mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers (see {ref}`mmCIF input <mmcif-input>`); `--no-convert-files` writes only the `.xyz` files.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `path-search`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | paths | (required) | Two or more structures in reaction order (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`), after one `-i` (`-i` may also be repeated for each file) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-b, --backend` | text | `uma` | Backend of the ML region (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `-o, --out-dir` | path | `./result_path_search/` | Output directory |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | Path method: growing string method / direct max flux |
| `--dmf-backend` | `gpu` / `cpu` | `gpu` | DMF compute backend (`--mep-mode dmf` only): PyTorch on CUDA / NumPy |
| `--refine-mode` | `peak` / `minima` | `peak` for GSM, `minima` for DMF | How the region around each HEI is relaxed: HEI ± 1 / nearest local minima |
| `--max-depth` | integer | `10` | Maximum levels of recursive subdivision; `0` turns subdivision off |
| `--max-nodes` | integer | `20` | Movable images per segment; a segment has `max_nodes + 2` images |
| `--preopt/--no-preopt` | flag | `True` | Pre-optimize each input before the search |
| `--align/--no-align` | flag | `True` | Align each input to the one before it before the search |
| `--freeze-atoms` | text | `None` | Comma-separated 1-based atom indices to freeze, added to YAML `geom.freeze_atoms` and the Frozen-MM layer (see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`) |
| `--climb/--no-climb` | flag | `True` | Run the GSM climbing-image search on the reactive segments; connecting paths never climb |

See the [generated CLI reference](reference/commands/path_search.md) for every option.

> **Note:** In YAML (`--config`), `search.max_depth` sets the depth limit when `--max-depth` is not given, `search.kink_max_nodes` (default `3`) sets the number of nodes inserted in a kink, and `bond.bond_factor` (default `1.20`) scales the covalent radii used to decide whether a bond has changed. Every key is listed under [`search`](yaml-reference.md#search) and [`bond`](yaml-reference.md#bond) in the YAML Reference; [`stopt`](yaml-reference.md#stopt) also lists `stopt.lbfgs` and `stopt.rfo`, which set the single-structure optimizers as `opt.lbfgs` and `opt.rfo` do.

---

## Notes

* **Inputs**: give at least two structures, all with the same atoms in the same order as `--parm7`; fewer than two stops with an error.
* **Templates for XYZ inputs**: `--ref-pdb` takes one full-system PDB per input, in the same order as `-i`; an `.xyz` input without its template stops with an error.
* **No climbing between segments**: `--climb` applies to the reactive segments; the short paths that connect neighboring pieces always run without climbing.
* **Inputs are protected**: if a fixed output name (`mep_trj.*`, `mep_plot.png`, `energy_diagram_MEP.png`, `summary.json`, `summary.log`) would replace an input file, `path-search` stops before writing anything.
* **Conflicting optimizer settings in YAML**: setting the same key to different values in `opt:` and in the section of the optimizer that runs (`lbfgs:`, `opt.lbfgs:`, `stopt.lbfgs:`, or the `rfo` equivalents) stops the run with an error.
* **Frozen atoms move slightly with DMF**: DMF holds frozen atoms with a harmonic restraint (k = 300 eV/Å², YAML `dmf.k_fix`), so they can drift a little; see [path-opt](path-opt.md#notes) and {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`.
* **DMF needs `cyipopt` and `pydmf`**: neither is installed with `mlmm-toolkit`; install them before you run `--mep-mode dmf` (see [path-opt](path-opt.md#notes)).
* **Complex mechanisms** may need adjusted intermediates, scan settings, or convergence thresholds.
* **Option priority**: default < YAML < command line (see {ref}`Configuration precedence <configuration-precedence>`).

---

## See also

* [path-opt](path-opt.md) — single-pass MEP between two structures
* [scan](scan.md) — drive a bond step by step to make a path or a TS candidate
* [tsopt](tsopt.md) — optimize each segment HEI into a TS
* [Building the ML region and layers](model-setup.md) — make the full-system PDB, `real.parm7`, and the ML region used as inputs
* [all](all.md) — the full workflow; `all --refine-path` runs `path-search` for its MEP step
* [YAML Reference](yaml-reference.md) — every `search`, `bond`, `gs`, and `dmf` setting
* [Glossary](glossary.md) — MEP, GSM, DMF, HEI, kink, and other terms
* [Troubleshooting](troubleshooting.md) — when a run fails
* {ref}`Exit codes <exit-codes>` — what each exit status means
