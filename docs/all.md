# `all` (end-to-end workflow)

`all` runs the whole ML/MM workflow in one command: it selects the ML region around the active site, builds the Amber topology and the three layers of the full system, and finds the minimum energy path (MEP). When asked, it optimizes the transition state (TS) of each reaction step and runs the intrinsic reaction coordinate (IRC), frequency, and DFT calculations on it.

Without `--tsopt`, the run ends with TS candidates: the highest-energy image (HEI) of each MEP segment. The ML region is computed by **UMA**, Meta's pretrained [machine-learning interatomic potential (MLIP)](backends.md), by default; `-b/--backend` also selects **ORB**, **MACE**, **AIMNet2**, or [DFT](dft-backend.md) (`dft`). The rest of the enzyme is computed with the Amber force field, and ONIOM combines the two.

---

## What it is for

What you pass selects the mode:

* **Path and energy diagram from R and P (Endpoint mode)**: give two or more full structures in reaction order (reactant, intermediates, product); `all` finds the MEP between each neighbouring pair and draws the energy diagram.
* **Path from a reactant alone (Scan-list mode)**: give one structure and the bonds to form or break with `-s`; a staged scan makes the intermediates, and the MEP search runs through them.
* **Check one TS candidate (TS-only mode)**: give one structure with `--tsopt` and no `-s`; `all` optimizes the TS and runs IRC from it. The TS is confirmed when n_imag = 1 and the IRC ends at the intended R and P.

---

## Examples

The examples use the GPP C6-methyltransferase BezA ([Tsutsumi et al., *Angew. Chem. Int. Ed.* 2022, 61, e202111217](https://doi.org/10.1002/anie.202111217)); the full scripts are in [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza). `1.R.pdb` (reactant), `2.IM.pdb` (intermediate), and `3.P.pdb` (product) are full enzyme structures with every hydrogen; your own structures need [hydrogens](getting-started.md) too. Examples 1–3 are walked through, with how to check the results, in [Quickstart: `all`](quickstart-all.md), [Quickstart: `--scan-lists`](quickstart-scan.md), and [Quickstart: TS-only mode](quickstart-tsopt.md).

### 1. MEP with TS optimization, thermochemistry, and DFT

`-c` names the residues the ML region is built around, and `-l` gives the charges of the non-standard residues. `--refine-path` splits the path where bonds change, so each chemical step becomes its own segment.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --refine-path --tsopt --thermo --dft --out-dir ./result_mep
```

Every requested stage finished when the console prints `[Imaginary modes] n=1 (...)` for each TS and `Scientific status: success` under the last `====== Pipeline summary ======`. Then check the endpoints as in [Reading the run status](#reading-the-run-status). The optimized structures are in `result_mep/segments/seg_NN/`.

### 2. Path from the reactant by a staged scan

Stage 1 brings the methyl carbon of SAM (CS1) to C7 of GPP (1.50 Å) and away from SD of SAM (3.30 Å), and stage 2 moves H11 of GPP away from C7 (2.90 Å) onto OE2 of Glu186 (1.00 Å).

```bash
mlmm all -i 1.R.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -s '[("CS1 SAM 320","C7 GPP 321",1.50),("CS1 SAM 320","SD SAM 320",3.30)]' \
       '[("C7 GPP 321","H11 GPP 321",2.90),("OE2 GLU 186","H11 GPP 321",1.00)]' \
    --tsopt --thermo --out-dir ./result_scan
```

The targets inside one literal move together in one stage. Literals given in a row run as successive stages, each starting from the end of the one before, and the stage ends become the inputs of the MEP search. Give `-s` once and list every literal after it. To decide how to split a reaction, see [Tips for studying reaction mechanisms](mechanism-tips.md). In a PDB with an empty chain field, write an atom as its residue name, residue number, and atom name in any order (`"CS1 SAM 320"`); with chains, write `A:SAM:320:CS1`. All accepted forms are in [Common options and selectors](cli-conventions.md).

### 3. Check a TS candidate (TS-only mode)

One input with `--tsopt` and no `-s` skips the MEP search.

```bash
mlmm all -i TS_candidate.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --tsopt --thermo --dft
```

The optimized R, TS, and P are written to `result_all/segments/seg_01/`.

### 4. Resume post-processing from a segment

To redo the post-processing from segment N, repeat the original command with the same inputs, extraction, layer, path, and calculator options and the same `--out-dir`, and add `--resume-segment N`. Post-processing options such as `--tsopt-max-cycles` may change.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    --refine-path --tsopt --thermo --dft \
    --resume-segment 1 --out-dir ./result_mep
```

The segments before N are kept; the post-processing from segment N onward, the summary, and the diagrams are written again.

---

## How it works

```text
Full structure(s) (PDB / mmCIF, or XYZ with --ref-pdb)
  ├─ (with -c) ML-region selection: extract
  │   └─ ml_region.pdb
  ├─ Amber topology of the full system: mm-parm (skipped with --parm7)
  │   └─ mm_parm/<input>.parm7
  ├─ three layers in the B-factors: define-layer
  │   └─ layered/<input>_layered.pdb
  ├─ (one structure with -s) staged scan: scan
  │   └─ stage ends as intermediates
  ├─ MEP search: path-opt (default) or path-search (--refine-path)
  │   └─ mep_trj.xyz and energy_diagram_MEP.png
  └─ (with --tsopt) TS optimization and IRC: tsopt → irc
      ├─ (with --thermo) frequencies and thermochemistry: freq
      └─ (with --dft) DFT single points: dft
```

1. **Preparing the input and the ML region**: for a PDB with blank element columns, `all` fills them in. With `-c`, it cuts out the active-site model around the given residues; the model of the first input becomes the ML region and is written to `ml_region.pdb`.
2. **Building the topology and the layers**: `mm-parm` builds the Amber topology of the full system from the first input with AmberTools (ff19SB, GAFF2 for non-standard residues); `--parm7` skips this step. `define-layer` then writes the three layers into the B-factors: ML (0), Movable-MM residues within 8 Å of the ML region (10), and Frozen-MM (20). Every later calculation runs on the full system with ML/MM, and the bonds cut at the ML/MM boundary are capped with link hydrogens.
3. **Building the path**: the input structures are optimized first (`--preopt`). With `-s`, the staged scan makes the intermediates. `path-opt` then finds the MEP between each neighbouring pair by GSM (growing string method, the default) or DMF (direct max flux); with `--refine-path`, the recursive `path-search` refines the path and splits it into steps where bonds change. The HEI of each step is its TS candidate.
4. **Optimizing the TS and following the IRC** (`--tsopt`): each HEI is optimized by RS-P-RFO (restricted-step partitioned rational function optimization) by default, and the final Hessian of the ML and movable MM atoms (PHVA, partial Hessian vibrational analysis) gives n_imag. From the TS, the IRC is traced in both directions with EulerPC (an Euler predictor–corrector integrator), and both ends are optimized to minima. These become the R and P of the segment.
5. **Thermochemistry and DFT**: `--thermo` runs `freq` on R, TS, and P for the ML/MM Gibbs energy, and `--dft` computes the ML region of the same structures with DFT and combines it with the MM energy. Each adds its own energy diagram.

---

## Reading the run status

A successful TS optimization gives one imaginary mode along the reaction coordinate (n_imag = 1). `all` continues from the TS to IRC only when the TS optimization converged, its final Hessian was computed, and n_imag ≥ 1:

| TS result | What `all` does |
| --- | --- |
| Converged, n_imag = 1 | Runs IRC and optimizes both IRC ends. |
| Converged, n_imag ≥ 2 | Warns, then runs IRC along the imaginary mode that best matches the MEP direction (the lowest one when none matches). The result is `partial`. |
| Converged, n_imag = 0 | Stops before IRC. |
| Stopped at the cycle limit | Computes no final Hessian and stops before IRC. |
| Stopped by `--stop-plateau` (stalled) | Computes the final Hessian, reports n_imag, and stops before IRC. |
| `--skip-final-freq`, or a failed Hessian | Stops before IRC. |

The full table of how a TS optimization can end is in [`tsopt` → Reading the TS result](tsopt.md#reading-the-ts-result).

Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

Read the outcome in three places:

* **Console**: each TS with one imaginary mode prints `[Imaginary modes] n=1 (...)` with its frequency. The `====== Pipeline summary ======` block prints `Execution status:` and `Scientific status:`. When the result is not `success`, `RESULT WARNING:` lines give the reasons.
* **`summary.log`**: the header shows `Pipeline mode` (`MEP`, `Scan`, or `TS-only`) and both statuses. Section [1] is the MEP overview; [2] lists the barrier ΔE‡, the reaction energy ΔE, and the bond changes of each segment on the MEP; [3] gives the post-processing of each segment, with n_imag under `TS imaginary freq:`; [4] tabulates the energy diagrams; [5] shows the output tree.
* **`summary.json`**: `scientific_status` is `success` when every requested stage converged, otherwise `partial` or `failed` with the reasons in `scientific_status_reasons`. n_imag of each TS is `post_segments[].tsopt.n_imaginary_modes`. The status fields are described in [Execution and requested-stage completion](json-output.md#execution-and-requested-stage-completion).
  * **Barriers**: with `--tsopt`, the barrier of each segment is `post_segments[].mlip.barrier_kcal`, the ML/MM energy of the optimized TS minus R. With `--thermo` and `--dft`, the barriers of the other methods are in the same form under `gibbs_mlip`, `dft`, and `gibbs_dft_mlip`. `segments[].barrier_kcal` is the barrier on the MEP before TS optimization; in TS-only mode it is TS − R.

Whether the endpoints are the intended R and P is for you to check: compare the bond changes in section [2] of `summary.log` and the structures `segments/seg_NN/reactant.*` and `product.*` with the R and P you intended. If n_imag is not 1 or the IRC ends are not the intended ones, see {ref}`When the TS search fails <ts-search-fails>`.

---

## Output files

`all` writes these files to `--out-dir`:

```text
result_all/
├─ summary.log                  # Text summary
├─ summary.json                 # Machine-readable results (always written; all has no --out-json)
├─ mep_trj.xyz                  # MEP trajectory over all segments
├─ mep_trj.pdb                  # Same trajectory as PDB
├─ mep_plot.png                 # ML/MM energy along the MEP trajectory
├─ energy_diagram_MEP.png       # MEP energy profile over all segments
├─ energy_diagram_*_all.png     # R → TS → P diagrams over all segments (--tsopt, --thermo, --dft)
├─ irc_plot_all.png             # IRC profiles over all segments (--tsopt)
├─ ml_region.pdb                # ML region (reusable with --model-pdb)
├─ ml_region_without_linkH.xyz  # ML region without and with the link hydrogens
├─ ml_region_with_linkH.xyz     #   (.pdb versions too for PDB input)
├─ mm_parm/                     # Amber topology <input>.parm7 and .rst7 (reusable with --parm7)
├─ layered/                     # Full structures with the three layers in the B-factors
├─ segments/
│  └─ seg_NN/                   # One reaction step: seg_01, seg_02, ...
│     ├─ reactant.*             # Optimized R, TS, and P in the input format (--tsopt)
│     ├─ ts.*
│     ├─ product.*
│     ├─ energy_diagram_*.png   # R → TS → P diagrams of this step
│     ├─ ts/                    # TS optimization; vib/imag_*_trj.xyz animates the imaginary modes
│     ├─ irc/                   # IRC trajectories and irc_plot.png
│     ├─ endpoint_opt/          # Endpoint optimizations (kept with --dump or when an endpoint did not converge)
│     ├─ freq/{R,TS,P}/         # Frequencies and thermochemistry (--thermo)
│     └─ dft/{R,TS,P}/          # DFT single points (--dft)
└─ _work/                       # Intermediate files, including the TS candidates (HEI)
   ├─ pockets/                  # Extracted models, pocket_<input>.pdb (with -c)
   ├─ scan/                     # Staged scan (with -s)
   └─ path_opt/                 # MEP search and hei_seg_NN.* (path_search/ with --refine-path)
```

* **Structures to report**: cite `segments/seg_NN/reactant.*`, `ts.*`, and `product.*`. The subdirectories of `seg_NN/` hold the files of each stage.
* **`.cif` files**: for mmCIF input and for PDB input too large for the PDB columns, `all` also writes `.cif` files that keep the original identifiers; see {ref}`mmCIF and large structures <mmcif-input>`.
* **TS-only mode**: there is no MEP search, so the MEP files and `_work/path_opt/` are absent; R, TS, and P go to `segments/seg_01/`.

The energy diagrams are named by method:

| File | Written when | Content |
| --- | --- | --- |
| `energy_diagram_MEP.png` | The MEP search finishes | MEP energy profile over all segments |
| `energy_diagram_MLIP.png` | `--tsopt` | R → TS → P, ML/MM energy |
| `energy_diagram_G_MLIP.png` | `--thermo` | R → TS → P, ML/MM Gibbs energy |
| `energy_diagram_DFT.png` | `--dft` | R → TS → P, DFT energy of the ML region on the ML/MM geometries |
| `energy_diagram_G_DFT_plus_MLIP.png` | `--dft` and `--thermo` | R → TS → P, ML(DFT)/MM energy plus the ML/MM thermal correction |
| `energy_diagram_*_all.png` | Same as the diagram without `_all` | The same diagram over all segments, at the top of the output directory |
| `irc_plot.png` (in `seg_NN/irc/`), `irc_plot_all.png` | `--tsopt` | IRC energy profile of one segment, and of all segments |

Energies in the diagrams are in kcal/mol relative to the first state (the reactant).

---

## Main options

`all` builds the Amber topology and the layers itself, so the shared {ref}`ML/MM options <mlmm-options>` are needed only to reuse the topology and `ml_region.pdb` of a previous run; the table below lists only the options specific to `all`.

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `-i, --input` | path(s) | (required) | Two or more full structures in reaction order, or one structure with `-s` or `--tsopt` (`.pdb`, `.cif`, `.mmcif`, or `.xyz` with `--ref-pdb`). Give several files after one `-i`, or repeat `-i` |
| `-c, --center` | text | `None` | Residues the ML region is built around, normally the substrate and catalytic residues: residue names (`'SAM,GPP'`), residue IDs (`'123,124'`, `'A:123,B:456'`), or a PDB file. Omit to take the ML region from the input B-factors or `--model-pdb` |
| `-l, --ligand-charge` | text | `None` | Charges of non-standard residues (e.g. `'SAM:1,GPP:-3'`), or their total charge as one number. Used for both the ML-region charge and the topology |
| `-q, --charge` | integer | `None` | Charge of the ML region. It is derived from the ML region; an explicit value overrides it with a warning |
| `-m, --multiplicity` | integer | `1` | Spin multiplicity (2S+1) of the ML region |
| `-b, --backend` | text | `uma` | Backend of the ML region (`uma`, `orb`, `mace`, `aimnet2`, `dft`) |
| `-r, --radius` | float | `2.6` | Extraction cutoff (Å) around the center atoms. `0` keeps only the `-c` and `--selected-resn` residues |
| `--selected-resn` | text | `""` | Residues to include without radius expansion: IDs (`'123'`, `'A:123A'`), names (`'SAM'`), or chain-qualified names (`'A:SAM'`, `'A:SAM:123'`) |
| `--auto-mm-ff-set` | `ff19sb` / `ff14sb` | `ff19sb` | Force field of the topology: ff19SB with OPC3 water, or ff14SB with TIP3P |
| `--auto-mm-add-ter/--no-auto-mm-add-ter` | flag | `True` | Insert TER records around ligand, water, and ion blocks and between disconnected peptides before building the topology |
| `--auto-mm-disulfide/--no-auto-mm-disulfide` | flag | `True` | Bond cysteines found by their SG–SG distance and rename them CYX; when off, only residues already named CYX are bonded |
| `--auto-mm-ligand-mult` | text | `None` (1 for every ligand) | Spin multiplicities of the ligands for the topology (e.g. `'GPP:2,SAM:1'`) |
| `--auto-mm-keep-temp` | flag | `False` | Keep the temporary directory of the topology build |
| `-s, --scan-lists` | text | `None` | Staged scan targets for one input, one literal per stage (e.g. `'[("A:SAM:320:CS1","A:GPP:321:C7",1.50)]'`; format in {ref}`Scan-list spec <scan-list-spec>`) |
| `--tsopt/--no-tsopt` | flag | `False` | Optimize the TS of each segment and run IRC |
| `--thermo/--no-thermo` | flag | `False` | Frequencies and ML/MM thermochemistry on R, TS, and P (needs `--tsopt`) |
| `--dft/--no-dft` | flag | `False` | DFT single points of the ML region on R, TS, and P (needs `--tsopt`) |
| `--refine-path/--no-refine-path` | flag | `False` | Run the recursive `path-search` instead of one `path-opt` per pair |
| `--mep-mode` | `gsm` / `dmf` | `gsm` | MEP method: GSM or DMF |
| `--opt-mode` | `grad` / `hess` | `grad` | Optimizer for the single-structure optimizations and the scan: `grad` = L-BFGS, `hess` = RFO. When given and `--opt-mode-post` is not, it is also used for the TS and the endpoints |
| `--opt-mode-post` | `grad` / `hess` | `hess` (or `--opt-mode` when that is given) | Optimizer for the TS and the endpoints after IRC: `grad` = Dimer for the TS and L-BFGS for the endpoints, `hess` = RS-P-RFO for the TS and RFO for the endpoints |
| `--preopt/--no-preopt` | flag | `True` | Optimize the input structures before the scan and the MEP search |
| `--flatten/--no-flatten` | flag | `False` | Remove extra imaginary modes left after the TS optimization |
| `--stop-plateau/--no-stop-plateau` | flag | `False` | Stop an optimization when the energy stops changing before convergence; the run is reported as stalled, not converged. The MM micro-iterations are not stopped this way |
| `--tsopt-max-cycles` | integer | `100000` | Cycle limit of the TS optimization |
| `--resume-segment` | integer | `None` | Redo the post-processing from segment N, reusing the MEP in `--out-dir` (example 4) |
| `--dry-run/--no-dry-run` | flag | `False` | Prepare the input and run the preflight checks in a temporary directory, print the plan, and skip the calculations |
| `-o, --out-dir` | path | `./result_all/` | Output directory |

For every option, run `mlmm all --help-advanced` or see the [generated CLI reference](reference/commands/all.md).

> **Note:** In YAML (`--config`), you can set what the options above do not cover; options given on the command line override the file. See [YAML Reference](yaml-reference.md) for the sections and keys.

---

## Notes

* **`--dft` and `-b dft`**: they cannot be used together, and the run stops with an error at startup. To add DFT single points after a `-b dft` run, run `mlmm sp -b dft` as a separate job.
* **Cost of `--dft`**: memory use depends on the size of the ML region, basis, functional, precision, and software stack. Try a representative structure on the target node and watch the peak memory. For a large ML region, finish the MLIP run first and run the DFT single points as a separate job.
* **R and P in TS-only mode**: the higher-energy IRC end is named the reactant. The names, the file names, the barrier, and the reaction energy follow this energy order, not a known chemical direction; the barrier from P is `barrier_kcal − delta_kcal`. `summary.json` records the rule under `endpoint_assignment`, with `chemical_direction_known: false`.
* **`summary.log` in TS-only mode**: [1] is the TS and IRC overview, and [2] comes from the optimized TS and endpoints.
* **When `all` stops before IRC**: the TS files stay in `segments/seg_NN/ts/`, and the later segments are not post-processed.
* **Endpoint optimizations**: if one endpoint optimization does not converge, the result is `partial` and `segments/seg_NN/endpoint_opt/` is kept for inspection. If an endpoint optimization fails with an error, the error is written to `segments/seg_NN/endpoint_opt/failure.json`, and that segment stops before the frequency and DFT stages; the TS and IRC structures are kept.
* **Thermochemistry file**: with `--thermo`, `thermoanalysis.yaml` is kept even under `--no-dump`, because `all` reads the thermochemistry from it.
* **Extraction radius**: `-r 0` disables radius-based expansion, so the model starts from the residues selected by `-c` and `--selected-resn`. Structural safeguards can still add a disulfide partner or the backbone of an adjacent residue. A zero radius is evaluated internally as 0.001 Å.
* **Without `-c`**: extraction is skipped, and the full input structures are used. The ML region comes from the B-factors of the input (`--detect-layer`, on by default) or from `--model-pdb`; with `--no-detect-layer` and no `--model-pdb`, the run stops with an error. One structure still needs `-s` or `--tsopt`.
* **AmberTools**: without `--parm7`, `all` stops with an error when AmberTools is not found.
* **Input formats**: `all` reads PDB and mmCIF; XYZ input needs `--ref-pdb`, a PDB with the same atoms. All structures of one run must have the same atoms in the same order.
* **Charge and multiplicity**: `-q` and `-m` describe the ML region, not the whole enzyme. With `-c`, the ML-region charge is the sum over the extracted model of the first input: built-in values for amino acids, ions, and water, `-l` for the other residues, and 0 for residues not listed in `-l`. Without `-c`, the same sum is taken over the ML region from the B-factors or `--model-pdb`. `-q` overrides the derived value with a warning; when no value can be derived, `calc.model_charge` in the YAML file is used. The multiplicity is `-m`, otherwise `calc.model_mult` in the YAML file, otherwise 1. See [Common options and selectors](cli-conventions.md).
* **Frozen atoms and rigid-body motions**: the vibrational analysis projects out only the rigid translations and rotations that leave the frozen atoms in place, so with a frozen MM layer usually none are removed; see [freq → Rigid modes with frozen boundaries](freq.md#rigid-modes-with-frozen-boundaries).
* **Separately prepared structures**: when the input structures were prepared independently, their differences outside the reaction coordinate enter the barrier. Compare the structures before reading the barrier. For two mechanisms of the same composition, use one common atom set and atom order for both paths.
* **`--resume-segment`**: it needs `--tsopt`, `--thermo`, or `--dft`, and cannot be combined with `--dry-run`. The run stops with an error when the saved inputs, ML region, topology, layered structures, or MEP do not match the command.

### Comparing a mutant with the wild type

Within one path, every structure has the same atoms in the same order. A mutant and the wild type (WT) differ in residues and often in atom count, so their total energies cannot be subtracted. Compare the barriers computed within each system instead:

`ΔΔG‡ = (G_TS − G_R)_mutant − (G_TS − G_R)_WT`

* Select the same ML-region residues and the same layer rules for both systems, so that the mutation is the only designed difference. Two independent radius-based selections can differ, because a boundary residue may enter one ML region and not the other; compare the two `ml_region.pdb` files.
* Use the same protonation rules, charge assignment, force field, backend and model, precision, restraints, and thermochemistry settings. If the mutation changes a protonation state or a formal charge in the ML region, the ML-region charges differ; do not force the same `-q` on both.

The two runs use the same options except for the input and the output directory. Give R and P of each system (Endpoint mode), so that R is the chemical reactant; `G_TS − G_R` is `post_segments[].gibbs_mlip.barrier_kcal`:

```bash
mlmm all -i wt_R.pdb wt_P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir ./result_wt
mlmm all -i mutant_R.pdb mutant_P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' --tsopt --thermo --out-dir ./result_mutant
```

---

## See also

* [extract](extract.md) — selection of the ML region
* [mm-parm](mm-parm.md) — Amber topology of the full system
* [define-layer](define-layer.md) — the three layers in the B-factors
* [scan](scan.md) — staged scans of distances, angles, and dihedrals
* [path-opt](path-opt.md) — one MEP between two structures (GSM / DMF)
* [path-search](path-search.md) — recursive MEP search that splits the path into steps
* [tsopt](tsopt.md) — TS optimization
* [irc](irc.md) — IRC from a TS
* [freq](freq.md) — vibrational analysis and thermochemistry
* [dft](dft.md) — DFT single points of the ML region
* [Refine an MLIP TS with DFT](dft-backend.md) — `-b dft` and `--dft`
* [Tips for studying reaction mechanisms](mechanism-tips.md) — splitting the reaction, checking the TS, and what to try when it fails
* [Troubleshooting](troubleshooting.md) — what to do when a run fails
* [Getting Started](getting-started.md) — the shortest run and what to read next
