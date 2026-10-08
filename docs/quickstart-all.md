# Quickstart: `mlmm all`

`mlmm all` builds a reaction path from the reactant (R) and product (P) in one run. It cuts out the ML region around the substrates, builds the Amber topology of the full system (`mm-parm`) and its three layers (`define-layer`), and searches the minimum energy path (MEP) between R and P with ML/MM. With `--tsopt --thermo --dft`, the same run continues to transition-state (TS) optimization, an intrinsic reaction coordinate (IRC) calculation, frequencies, and DFT single points.

The commands below use the bundled example of the GPP (geranyl pyrophosphate) C6-methyltransferase BezA in [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza): `1.R.pdb` is the reactant, `2.IM.pdb` an intermediate, and `3.P.pdb` the product. They are full enzyme structures with every hydrogen. Get them with `git clone https://github.com/t-0hmura/mlmm_toolkit && cd mlmm_toolkit/examples/beza`. For your own reaction, replace them with your full-system structures.

---

## What it is for

* **A first run of the whole workflow**: run every stage once on the bundled example.
* **The MEP between R and P**: get the path and its highest-energy image (HEI), the TS candidate.
* **TS, IRC, frequencies, and DFT in the same run**: add `--tsopt --thermo --dft` to check the TS candidate.

## Minimal command

Give R and P in reaction order, the residues the ML region is built around (`-c`), and the ligand charges (`-l`).

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --out-dir ./result_all
```

The run succeeded when the `====== Pipeline summary ======` block near the end of the console shows `Scientific status: success`; `summary.json` holds the same value in `scientific_status`.

### (Optional) Add post-processing in the same run

`--tsopt` adds TS optimization and IRC for each [reactive segment](glossary.md) (here `seg_01`), `--thermo` adds frequencies and thermochemistry, and `--dft` adds DFT single points of the ML region on R, TS, and P.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
 --tsopt --thermo --dft --out-dir ./result_all
```

## Before you run

The structures need every hydrogen atom, and R and P must list the same atoms in the same order; see [Before you run: the input structures](getting-started.md#before-you-run-the-input-structures). Two more points apply to the topology build:

* **Charges and hydrogens**: give each ligand in `-l` the charge that matches its hydrogens in the file. In the bundled example SAM has 23 hydrogens, so it is `SAM:1`; with 22 it would be `SAM:0`. When they do not match, `mm-parm` stops with an electron-count error before it runs `antechamber`.
* **AmberTools**: `all` builds the topology with AmberTools (`tleap`, `antechamber`, `parmchk2`) and stops with an error when they are missing. `--parm7` with an existing topology skips the build.

## Output files

The minimal command writes:

```text
result_all/
├── summary.log                  # Run summary
├── summary.json                 # Results, with scientific_status
├── mep_trj.pdb                  # MEP over all segments
├── energy_diagram_MEP.png       # MEP energy profile over all segments
├── ml_region.pdb                # ML region (reusable with --model-pdb)
├── mm_parm/                     # Amber topology 1.R.parm7 and 1.R.rst7 (reusable with --parm7)
├── layered/                     # 1.R_layered.pdb and 3.P_layered.pdb, with the three layers in the B-factors
└── _work/                       # Intermediate files, including the HEI (TS candidate); kept after the run
    └── path_opt/                # MEP search (path_search/ with --refine-path, the recursive MEP search)
        ├── hei_seg_01.{xyz,pdb} # Highest-energy image of segment 1
        └── summary.json         # MEP search results
```

The minimal command stops after the MEP search and does not create `segments/`. With `--tsopt`, a reactive segment adds `segments/seg_NN/` with the R/TS/P structures (`reactant.pdb`, `ts.pdb`, `product.pdb`), `ts/`, and `irc/`; `--thermo` also adds `freq/`, and `--dft` adds `dft/`. With mmCIF input, `all` also writes `.cif` files such as `mep_trj.cif`.

## Checking the result

1. **Completion**: `scientific_status` is `success` when every requested stage converged; otherwise it is `partial` or `failed`, with the [reasons](json-output.md#execution-and-requested-stage-completion) in `scientific_status_reasons`. With `--tsopt`, two checks are left for you: that the imaginary mode moves the bonds that form or break, and that the endpoints are the intended R and P.
2. **TS candidate**: open `_work/path_opt/hei_seg_01.pdb`, the HEI of the first segment. It is a full-system PDB with the layers in the B-factors. With `--tsopt`, also open the optimized TS, `segments/seg_01/ts.pdb`.
3. **Energy profile**: `energy_diagram_MEP.png` should show a clear barrier between R and P.
4. **TS (with `--tsopt`)**: a successful TS optimization gives one imaginary mode along the reaction coordinate. The console then prints `[microiter] Converged!` and then `[Imaginary modes] n=1 (...)`, with the imaginary wavenumber in cm⁻¹ in the brackets. Open `segments/seg_01/ts/vib/imag_*_trj.xyz` in a viewer and check that the mode moves the bonds that form or break.
5. **Endpoints (with `--tsopt`)**: open `segments/seg_01/irc/finished_irc_trj.xyz` and the optimized endpoints `segments/seg_01/reactant.pdb` and `product.pdb`, and check that they are the intended R and P. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

For how `all` judges each stage, see [Reading the run status](all.md#reading-the-run-status).

## Notes

* **Reusing the preparation**: pass `mm_parm/1.R.parm7` with `--parm7` and `ml_region.pdb` with `--model-pdb` to the next run or to the single-stage commands, so that they compute the same system without building the topology again.
* **DFT and GPU memory**: `--dft` needs the DFT extra (step 7 of {ref}`Step-by-step installation <step-by-step-installation>`); for GPU memory, see the Notes of [Refine an MLIP TS with DFT](dft-backend.md#notes).
* **Barriers in `summary.json`**: `segments[].barrier_kcal` is the barrier on the MEP, before TS optimization. With `--tsopt`, `post_segments[].mlip.barrier_kcal` is the ML/MM barrier from the optimized TS and endpoints; `--thermo` adds `post_segments[].gibbs_mlip.barrier_kcal` and `--dft` adds `post_segments[].dft.barrier_kcal`. `rate_limiting_step.barrier_kcal` is the highest barrier among the segments, compared at the highest level that every segment has (`DFT//MLIP/MM_Gibbs` > `DFT` > `MLIP_Gibbs` > `MLIP` > `MEP`); `rate_limiting_step.method` names that level.
* **Run time**: it depends on the system size, the size of the ML region, the GPU, and the stages you request.

## Next steps

- [Quickstart: scan](quickstart-scan.md): start from one structure when there is no product structure
- [Quickstart: TS-only mode](quickstart-tsopt.md): optimize and check a TS candidate you already have
- [Building the ML region and layers](model-setup.md): shrink the ML region, or extend it when residues are missing
- [Tips for studying reaction mechanisms](mechanism-tips.md): plan the calculations, and what to try when the TS search fails
- [Refine an MLIP TS with DFT](dft-backend.md): refine and check the TS with DFT/MM
- [`all`](all.md): full option reference (also `mlmm all --help-advanced`)
- [JSON Output Reference](json-output.md): the fields of `summary.json`
- [Troubleshooting](troubleshooting.md): find an error message or symptom and its fix
