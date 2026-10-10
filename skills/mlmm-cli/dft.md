# `mlmm dft`

Runs a DFT single point on the ML region of one ML/MM structure with
GPU4PySCF (GPU) or PySCF (CPU), and combines it with the MM energies into the
ML(DFT)/MM total energy. Run
`mlmm dft -i <structure.pdb> --parm7 <real.parm7> -l '<RES:Q,...>'`. Success
is `energy.converged: true` in `result.yaml` and the console line
`E_total ML(dft)/MM (Hartree): …`.

## When to use

- Single-point energies and atomic charges at a higher level than the MLIP,
  on R, TS, and P from `tsopt` and `irc`, or as a standalone DFT driver on any
  ML/MM structure. `mlmm dft` computes energies only, no forces.
- `-b dft` instead makes DFT the ML-region level of every step of a run. It
  works in `sp`, `opt`, `tsopt`, `irc`, `freq`, `scan`, `scan2d`, `scan3d`,
  `path-opt`, `path-search`, and `all`; iterative runs start each SCF from the
  last converged density. Keep the ML region to roughly 300 atoms, counting
  the link H. To refine an MLIP/MM TS, run [TS-only mode](all-ts-only.md)
  with `-b dft` on `segments/seg_NN/ts.pdb`, reusing the `--parm7` and
  `--model-pdb` of the first run; its frequencies and thermochemistry stay a
  PHVA of the ML and movable MM atoms, as in the MLIP/MM run. If you take the
  structure to another QM code instead (for example with
  [oniom-export](oniom.md)), first check that its frequency analysis also
  leaves out the frozen atoms; otherwise its n_imag and Gibbs corrections are
  not comparable.
- `all --dft` runs these single points on R, TS, and P after an MLIP/MM run
  and needs `--tsopt`; `-b dft` and `--dft` cannot be combined. How to choose:
  [DFT backend](../../docs/dft-backend.md).

## Minimal run

```bash
mlmm dft -i seg_01/ts.pdb --parm7 real.parm7 -l 'SAM:1,GPP:-3' \
    --func-basis 'wb97m-v/def2-svp' --dft-engine gpu -o result_dft_svp
```

CPU PySCF on an XYZ, with the full-system PDB as the template:

```bash
mlmm dft -i ts.xyz --ref-pdb real.pdb --parm7 real.parm7 -q 0 -m 1 \
    --func-basis 'wb97m-v/def2-svp' --dft-engine cpu -o result_dft_cpu
```

`-q` is the ML-region charge; `-l` derives it from the ligand charges and
the residues in the ML region. `--func-basis` (default `wb97m-v/def2-svp`)
takes PySCF names as `FUNC/BASIS`. `--embedcharge` puts the MM point charges
within `--embedcharge-cutoff` (12.0 Å) into the DFT Hamiltonian.
`--dft-nprocs` and `--dft-memory` (host RAM, not GPU memory) override the
values taken from the scheduler. The default output directory is
`./result_dft/`.

## Judge success

- The console prints `E_DFT (Hartree): …` and
  `E_total ML(dft)/MM (Hartree): …`; `result.yaml` has
  `energy.converged: true`, the MM energies, and the per-atom charges and
  spin densities.
- With `--out-json`, `result.json` (and `summary.json`, same content) has
  `converged`, `energy_hartree` (ML-region DFT energy),
  `total_dft_mm_energy_hartree`, `xc_functional`, `basis_set`, `engine`
  (`gpu4pyscf(rks_lowmem)`, `gpu4pyscf`, or `pyscf(cpu)`), `used_gpu`, and
  `used_lowmem`.
- The ML region is written as `ml_region_without_linkH.xyz` (the selection
  before link H) and `ml_region_with_linkH.xyz` (with the generated link H,
  as passed to PySCF; its first line is the atom count). PDB input with
  `--convert-files` also gives `ml_region_without_linkH.pdb`, which keeps
  the topology identifiers, and `ml_region_with_linkH.pdb`, with each link H
  as atom `HL` of residue `LKH`.

## Choosing the engine

- `gpu` (default): GPU4PySCF on x86_64 with a supported CUDA stack.
  `--dft-low-memory` (on by default) uses the low-memory RKS for closed-shell
  GPU runs, including embedding; `--no-dft-low-memory` uses density fitting.
- `cpu`: PySCF, no GPU needed. Use it on aarch64 or other machines without
  prebuilt GPU4PySCF wheels, or when no supported GPU stack exists.
- If GPU4PySCF cannot run, `dft` stops with
  `[gpu] GPU backend failed: … Set dft.engine: cpu …`; it does not switch to
  the CPU by itself. Rerun with `--dft-engine cpu`.
- Cost depends on elements, basis, functional, grid, engine, and software
  stack. Time one representative structure and size the job from its
  measured peak memory.

## Pitfalls and recovery

- `OSError: libcusolver.so.11: cannot open shared object file`:
  `LD_LIBRARY_PATH` shadows the CUDA
  libraries bundled with torch. Fix the order as in
  [Library-loading collisions](../mlmm-install/backends.md#library-loading-collisions).
- `cupy ... invalid device ordinal`: keep the scheduler's
  `CUDA_VISIBLE_DEVICES` and use a valid local ordinal (usually 0 on a
  one-GPU allocation).
- `RuntimeError: CUDA out of memory`: run the same method on CPU or a GPU
  with more memory. A smaller basis or grid is a different method; label and
  revalidate it. On a new GPU generation, check the GPU4PySCF and CuPy
  versions and the traceback first.
- GPU startup stalls or fails: capture the traceback, run `pip check`, and
  compare with the requirements of the installed GPU4PySCF version.
- SCF not converged: `WARNING: SCF did not converge.`; the results are still
  written with `converged: false`, and the exit code is 1. In low-memory
  mode, retry with `--no-dft-low-memory` when memory allows.
- A new run first removes `result.yaml`, `result.json`, `summary.json`, and
  the four `ml_region_*` files in its output directory; use a new `-o` to
  keep earlier results.
- A `def2` basis attaches the def2 ECP (`[dft] Using ECP: …`).

## DFT//MLIP/MM on the TS candidate

Evaluate the converged TS and the optimized, chemically identified IRC
endpoints with the same settings. For standalone runs, optimize the IRC
endpoints with `opt` first. The TS:

```bash
mlmm dft -i result_tsopt/final_geometry.pdb --parm7 real.parm7 \
    -l 'SAM:1,GPP:-3' --func-basis 'wb97m-v/def2-tzvpd' \
    --dft-engine gpu -o dft_TS
```

Repeat for both endpoints, then combine the energies with
[energy-diagram](utilities.md#energy-diagram). After `all --tsopt`, the
structures are `segments/seg_NN/reactant.pdb`, `ts.pdb`, and `product.pdb`.

## Next step

- [DFT (PySCF, GPU4PySCF)](../mlmm-install/backends.md#dft-pyscf-gpu4pyscf):
  install the `[dft]` extra (`[dft-cuda12]` for the `cu126` wheel).
- [tsopt.md](tsopt.md) and [irc.md](irc.md): make the geometries.
- `--show-config` prints the YAML given with `--config` and its top-level
  keys, then continues.
