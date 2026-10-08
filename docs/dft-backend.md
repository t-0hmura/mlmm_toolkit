# Refine an MLIP TS with DFT

Once MLIP/MM has found a reasonable pathway, mlmm-toolkit can take its TS straight into a DFT/MM TS optimization. It runs the TS optimization → IRC → endpoint optimization → frequency workflow with GPU-accelerated DFT through GPU4PySCF. DFT computes only the ML region; the protein around it stays in MM.

The MLIP/MM pathway search is the main tool; DFT/MM is an add-on that checks the TS candidate you found with MLIP/MM.

## What it is for

- **Refine the TS at the DFT/MM level**: run TS optimization → IRC → endpoint optimization → frequencies with DFT for the ML region (`-b dft`).
- **Add DFT energies to an MLIP/MM run**: run DFT single points on the R, TS, and P from the MLIP/MM run (`--dft`).
- **Run a DFT single point on one structure**: get the DFT/MM energy with population analysis (`mlmm dft`), or the energy and forces (`mlmm sp -b dft`).

## Workflow

1. **Explore with MLIP/MM**: generate pathways, try variants, and pick the most promising TS candidate.
2. **Refine with DFT/MM**: run [TS-only mode](#examples) on that TS with `-b dft`, reusing the MM topology and the ML region of the first run.
3. **Check**: as in an MLIP/MM run, a successful TS optimization gives one imaginary mode along the reaction coordinate (`[Imaginary modes] n=1` in the log). `Scientific status: success` under `====== Pipeline summary ======` shows that every requested stage converged. Even if the IRC does not converge, the result is usable when the endpoint optimizations reach the intended R and P.

## Examples

### 1. Search with MLIP/MM on a small ML region

Run the MLIP/MM search with an ML region small enough for DFT. Its output directory, `result_all/`, holds the TS, the MM topology, and the ML region that example 2 reuses.

```bash
mlmm all -i 1.R.pdb 3.P.pdb -c 'SAM,GPP,MG' -l 'SAM:1,GPP:-3' \
    -r 0 --selected-resn '44,63,186' --tsopt
```

`-r 0` sets the extraction radius to 0 Å, which stops adding nearby residues by distance; the ML region is built from the `-c` and `--selected-resn` residues. In the bundled example in [`examples/beza/`](https://github.com/t-0hmura/mlmm_toolkit/tree/main/examples/beza), residues 44, 63, and 186 are the three closest to SAM's methyl carbon (CS1). The bundled PDB has an empty chain field; for such PDB files, give residues by name or number. For your own system, pick the residues that take part in the reaction.

### 2. Refine the TS with DFT/MM

Pass the TS from example 1 as the only input, which selects TS-only mode. `ts.pdb` holds the full system, and `--parm7` and `--model-pdb` reuse the topology and ML region of example 1, so both runs describe the same system. This needs the DFT extra (see [Notes](#notes)).

```bash
mlmm all -i result_all/segments/seg_01/ts.pdb \
    --parm7 result_all/mm_parm/1.R.parm7 --model-pdb result_all/ml_region.pdb \
    -l 'SAM:1,GPP:-3' --tsopt --thermo -b dft -o ./result_dft
```

`seg_01` is the first reaction segment; if there are several, pick the one for the step you want to refine.

## Keep the ML region under about 300 atoms

With `-b dft`, DFT computes the whole ML region, so keep it to roughly 300 atoms, counting the link hydrogens.

For ways to shrink it, see [Shrink the ML region](model-setup.md#shrink-the-ml-region).

## `-b dft` and `--dft`

| Option | What DFT computes | When to use it |
|---|---|---|
| `-b dft` | The ML region in every calculation of the run (MEP search, TS optimization, IRC, endpoint optimization, frequencies) | Refine and check a TS candidate at the DFT/MM level |
| `--dft` | Single points on the R, TS, and P from the MLIP/MM run (`all` only) | Get DFT energies on MLIP/MM geometries |

`-b dft` works in 11 commands: `all`, `opt`, `tsopt`, `irc`, `freq`, `scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `sp`.

## Output files

With `-b dft`, the output has the same layout as an MLIP/MM run in the same mode; for TS-only mode, see [Quickstart: TS-only mode](quickstart-tsopt.md). `--dft` adds these files:

| File | Content |
|---|---|
| `segments/seg_NN/dft/{R,TS,P}/` | DFT single-point result for each state |
| `segments/seg_NN/energy_diagram_DFT.png` | DFT energy diagram on the MLIP/MM geometries |
| `segments/seg_NN/energy_diagram_G_DFT_plus_MLIP.png` | DFT energy plus the MLIP/MM thermal correction (with `--thermo`) |
| `energy_diagram_DFT_all.png`, `energy_diagram_G_DFT_plus_MLIP_all.png` | The same diagrams over all segments, at the top of the output directory |

## Main options

| Option | Description | Default |
|---|---|---|
| `-b, --backend dft` | Use DFT for the ML region (GPU4PySCF; CPU PySCF with `--dft-engine cpu`). | `uma` |
| `--parm7 FILE`, `--model-pdb FILE` | Reuse the MM topology and the ML region of an earlier run. Without them, `all` builds them again from the input. | — |
| `--func-basis TEXT` | Functional and basis as `FUNCTIONAL/BASIS`. Applies to both `-b dft` and `--dft`. | `wb97m-v/def2-svp` |
| `--embedcharge/--no-embedcharge`, `--embedcharge-cutoff FLOAT` | With `-b dft`, put the MM point charges within the cutoff (Å) of the ML region into the DFT Hamiltonian. | `--no-embedcharge`, `12.0` |
| `--dft/--no-dft` | Add DFT single points on R, TS, and P (`all` only). | `--no-dft` |

The other DFT options (`--dft-engine`, `--dft-low-memory/--no-dft-low-memory`, `--dft-nprocs`, `--dft-memory`, SCF checkpoints) are listed in the [`all` reference](reference/commands/all.md).

> **Note:** In YAML, the same settings go under `calc.dft`. `calc.dft.pyscf` passes attributes to PySCF objects by name, for example `mf: {level_shift: 0.2}` for a hard-to-converge SCF. See [YAML Reference](yaml-reference.md).

## Notes

- **DFT extra**: install it with `pip install "mlmm-toolkit[dft]"` for the `cu130` or `cu132` PyTorch wheel, or with `pip install "mlmm-toolkit[dft-cuda12]"` for `cu126`; see step 7 of {ref}`Step-by-step installation <step-by-step-installation>`. Without a GPU, add `--dft-engine cpu`.
- **Keep `--parm7` and `--model-pdb`**: without them, the DFT/MM run builds a new topology from the TS structure and takes the ML region from its B-factor layers. The parm7 is named after the first input of example 1 (`mm_parm/1.R.parm7` here); take the name from `result_all/mm_parm/`.
- **Charge**: a smaller ML region usually has a different charge. Check `Total active site model charge` in the console output of example 1 before the DFT/MM run.
- **300 atoms is a guide**: the code sets no limit. The first line of `ml_region_with_linkH.xyz` in the output directory is the ML-region atom count including the link hydrogens.
- **Combinations**: `-b dft` and `--dft` cannot be used together, and the run stops with an error at startup. To add DFT single points after a `-b dft` run, run `mlmm sp -b dft` or `mlmm dft` as a separate job. `--dft` and `--thermo` require `--tsopt`.
- **File and key names**: with `-b dft`, the diagrams keep the file names `energy_diagram_MLIP.png` and `energy_diagram_G_MLIP.png` (with `--thermo`), and the `summary.json` blocks keep the names `mlip` and `gibbs_mlip`. Both hold the DFT/MM values, and the plot title shows DFT/MM.
- **Memory and threads**: `--dft-memory` is the host RAM for PySCF, not the GPU VRAM. If GPU memory runs out, shrink the ML region; if `--dft` runs out of memory, drop `--dft` and run `mlmm dft` separately.

## See also

- [`dft`](dft.md): DFT single point with population analysis
- [`sp`](sp.md): single-point energy and forces with any backend
- [Quickstart: TS-only mode](quickstart-tsopt.md): check a TS candidate with `all --tsopt`
- [Building the ML region and layers](model-setup.md): build, shrink, and extend the ML region
- {ref}`Installation <step-by-step-installation>`: step 7 installs the DFT extra
- [MLIP Backends](backends.md): choosing a backend
- [Troubleshooting](troubleshooting.md): what to do when a run fails
