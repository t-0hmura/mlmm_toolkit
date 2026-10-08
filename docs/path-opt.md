# `path-opt` (MEP between two structures)

`path-opt` finds a minimum-energy path (MEP) between two layered enzyme structures, a reactant and a product, in one pass with GSM (growing string method, the default) or DMF (direct max flux), using the ML/MM calculator on the whole system. It writes the highest-energy image (HEI) as a TS candidate.

---

## What it is for

* **A first MEP from R and P**: get a path and its energy profile from two endpoint structures, without recursive refinement.
* **A TS candidate for `tsopt`**: `hei.pdb` (or `hei.xyz`) is the starting structure for [`tsopt`](tsopt.md).
* **Comparing GSM and DMF**: run the same pair with `--mep-mode gsm` and `--mep-mode dmf` and compare the paths.

The ML region is computed with **UMA** (Meta) by default, and the rest of the system is computed with the Amber force field of `--parm7`. For two or more structures with automatic refinement of the reactive region, use [`path-search`](path-search.md).

---

## Examples

### 1. Two endpoints

Give the reactant and the product after one `-i`, with the charge of the ML region and the spin multiplicity. `reactant.pdb` and `product.pdb` hold the whole system that matches `real.parm7` (the Amber topology), and `ml_region.pdb` selects its ML region.

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --out-json --out-dir ./result_path_opt
```

A console line that starts with `[write] Wrote '…/hei.xyz'` shows that the TS candidate was written. In `result.json`, `scientific_status` is `success` when every requested stage (endpoint pre-optimization and the MEP) converged, otherwise `partial` or `failed`. `barrier_kcal` is the HEI energy relative to the first image, and `hei_index` tells you where the HEI sits on the path (see [Reading the HEI](#reading-the-hei)).

### 2. Set the endpoint pre-optimization limit

Both endpoints are pre-optimized by default; `--preopt-max-cycles` caps each pass. Use `--no-preopt` when the endpoints are already optimized.

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --preopt-max-cycles 20000 --out-dir ./result_path_opt_preopt
```

### 3. DMF instead of GSM

DMF needs `cyipopt` and `pydmf` (see [Notes](#notes)); here we also use fewer movable images.

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --mep-mode dmf --max-nodes 12 --out-dir ./result_path_opt_dmf
```

### 4. Quick pass: no climbing, fewer nodes

Skip the climbing-image search and use fewer movable images for a fast first look at the path.

```bash
mlmm path-opt -i reactant.pdb product.pdb --parm7 real.parm7 --model-pdb ml_region.pdb \
  -q 0 -m 1 --no-climb --max-nodes 8 --out-dir ./result_path_opt_quick
```

---

## How it works

1. **Preparing the endpoints**:
Each endpoint is pre-optimized (L-BFGS by default). The product is then rigidly aligned to the reactant by a Kabsch fit on the frozen atoms (on all atoms if there are none), and the frozen atoms are moved step by step onto their reactant positions while the other atoms relax. The frozen atoms are the Frozen-MM layer together with any `--freeze-atoms`.
2. **Growing and refining the path**:
GSM grows a string of `--max-nodes` movable images between the two endpoints and optimizes it to `--thresh-gsm`. With `--climb` (on by default), a climbing-image search then pushes the highest image toward the saddle. DMF instead builds an interpolated path and optimizes it with IPOPT (an interior-point optimizer) to `--dmf-tol`.
3. **Writing the HEI**:
The image with the highest energy on the final path becomes the HEI. `hei.xyz` holds it with its energy on the comment line.

---

## Reading the HEI

| Where the HEI is | Meaning | Next step |
| --- | --- | --- |
| Inside the path (`hei_index` between `1` and `n_images − 2`) | A TS candidate | Optimize it with [`tsopt`](tsopt.md), then run [`irc`](irc.md) |
| At an endpoint (`hei_index` is `0` or `n_images − 1`) | Not a TS candidate: no image between the endpoints lies above the higher endpoint | Check the endpoints, or get a candidate another way (see {ref}`When a TS search fails <ts-search-fails>`) |

The HEI is the top of an approximate path, not a TS. A successful TS optimization gives one imaginary mode along the reaction coordinate. The HEI becomes a TS once `tsopt` gives n_imag = 1 and IRC from that TS reaches the intended R and P.

---

## Output files

`path-opt` writes these files to `--out-dir` (default `./result_path_opt/`):

```text
result_path_opt/
├─ final_geometries_trj.xyz   # Final path, every image, energies on the comment lines
├─ final_geometries.pdb       # Same path as PDB
├─ hei.xyz                    # HEI, the TS candidate, with its energy on the comment line
├─ hei.pdb                    # Same HEI as PDB
├─ preopt/endNN/              # Endpoint pre-optimization files (--preopt)
├─ align_refine/              # Endpoint alignment and relaxation files
├─ dmf_initial_trj.xyz        # Interpolated starting path (DMF only)
├─ result.json                # Summary (--out-json)
└─ summary.json               # Same content as result.json (--out-json)
```

Open `final_geometries_trj.xyz` to watch the path. The PDB files mark the layers in the B-factor column (ML 0, Movable-MM 10, Frozen-MM 20), so pass `hei.pdb` to `tsopt` with the same `--parm7` and `--model-pdb`; with `hei.xyz`, add `--ref-pdb`. mmCIF input, and PDB input too large for the PDB columns, also get `.cif` files that keep the original identifiers (see {ref}`mmCIF input <mmcif-input>`); `--no-convert-files` writes only the `.xyz` files. With DMF, the IPOPT logs `dmf_fbenm_ipopt.out` and `dmf_ipopt.out` are also written. `--dump` also keeps the optimizer trajectories.

The console prints the MEP progress cycle by cycle, with timings.

---

## Main options

The options shared by every ML/MM calculation command are explained once in {ref}`ML/MM options <mlmm-options>`; the table below lists only the options specific to `path-opt`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | 2 paths | (required) | Reactant and product, in that order, after one `-i` (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`) |
| `-q, --charge` | integer | `None` | Charge of the ML region. Required unless `-l` is given |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) |
| `-l, --ligand-charge` | text | `None` | Total ligand charge (for example `-1`) or a charge per residue name (for example `'GPP:-3,SAM:1'`), used to derive the ML-region charge when `-q` is omitted (PDB/mmCIF input or `--ref-pdb`) |
| `-b, --backend` | text | `uma` | Backend of the ML region (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `-o, --out-dir` | path | `./result_path_opt/` | Output directory |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | Path method: growing string method / direct max flux |
| `--dmf-backend` | `gpu` / `cpu` | `gpu` | DMF compute backend (`--mep-mode dmf` only): PyTorch on CUDA / NumPy |
| `--max-nodes` | integer | `20` | Movable images between the endpoints; the path has `max_nodes + 2` images |
| `--preopt/--no-preopt` | flag | `True` | Pre-optimize each endpoint before alignment |
| `--preopt-max-cycles` | integer | `100000` | Maximum cycles of each endpoint pre-optimization |
| `--opt-mode` | `grad` / `hess` | `grad` | Optimizer of the endpoint pre-optimization: L-BFGS / RFO |
| `--thresh-gsm` | preset | `gau_loose` | Convergence criteria of the GSM string (`gau_loose`, `gau`, `gau_tight`, `gau_vtight`, `baker`, `never`) |
| `--dmf-tol` | text | `tight` | IPOPT tolerance of the DMF path: `tight` (0.04), `middle` (0.10), `loose` (0.20), or a positive number; alias `--thresh-dmf` |
| `--fix-ends/--no-fix-ends` | flag | `True` | Keep the endpoints fixed while the GSM string is optimized (not used by DMF) |
| `--climb/--no-climb` | flag | `True` | Run the GSM climbing-image search after the path is grown (not used by DMF) |
| `--freeze-atoms` | text | `None` | Comma-separated 1-based atom indices to freeze in every image, added to YAML `geom.freeze_atoms` and the Frozen-MM layer (see {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`) |
| `--out-json/--no-out-json` | flag | `False` | Write a summary to `result.json` ([JSON Output Reference](json-output.md)) |

See the [generated CLI reference](reference/commands/path_opt.md) for every option.

> **Note:** In YAML (`--config`), the [`gs`](yaml-reference.md#gs) section sets the GSM string, [`dmf`](yaml-reference.md#dmf) the DMF path, and [`stopt`](yaml-reference.md#stopt) the string optimizer. `stopt.lbfgs` and `stopt.rfo` also set the endpoint optimizers, as `opt.lbfgs` and `opt.rfo` do.

---

## Notes

* **Frozen atoms move slightly with DMF**: DMF holds frozen atoms with a harmonic restraint (k = 300 eV/Å², YAML `dmf.k_fix`) instead of fixing them, so they can drift a little from their reference positions; GSM keeps them fixed. See {ref}`Freeze atoms and restrain distances <freeze-atoms-and-restraints>`.
* **DMF needs `cyipopt` and `pydmf`**: neither is installed with `mlmm-toolkit`. Install them before you run `--mep-mode dmf`: `conda install -c conda-forge cyipopt -y`, then `pip install 'pydmf[torch]>=1.2'` for the GPU backend or `pip install 'pydmf>=1.2'` for the CPU backend. The default `--dmf-backend gpu` stops with an error when CUDA is unavailable; use `--dmf-backend cpu` then, or after a GPU out-of-memory error.
* **Options DMF ignores**: `--climb` and `--fix-ends` are accepted but not used by DMF.
* **One template for XYZ endpoints**: `--ref-pdb` takes a single full-system PDB, applied to both `.xyz` endpoints.
* **Conflicting optimizer settings in YAML**: setting the same key to different values in `opt:` and in the section of the optimizer that runs (`lbfgs:`, `opt.lbfgs:`, `stopt.lbfgs:`, or the `rfo` equivalents) stops the run with an error.
* **Option priority**: default < YAML < command line (see {ref}`Configuration precedence <configuration-precedence>`).
* **Exit status**: see {ref}`Exit codes <exit-codes>`.

---

## See also

* [path-search](path-search.md) — MEP through two or more structures, refined where bonds change
* [tsopt](tsopt.md) — optimize the HEI into a TS
* [irc](irc.md) — check that the TS connects the intended R and P
* [all](all.md) — the full workflow; its MEP step uses `path-opt`, and `--refine-path` switches it to `path-search` (see [Main options](all.md#main-options))
* [YAML Reference](yaml-reference.md) — every `gs`, `dmf`, and `stopt` setting
* [Glossary](glossary.md) — MEP, GSM, DMF, HEI, and other terms
* [Troubleshooting](troubleshooting.md) — when a run fails
* {ref}`Exit codes <exit-codes>` — what each exit status means
