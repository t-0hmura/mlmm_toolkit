# Troubleshooting

Find your symptom in the quick table, then read the fix in the section it points to.

(troubleshooting-quick-table)=
## Quick routing

| Symptom | Start here | Then read |
| --- | --- | --- |
| **Input & extraction** | | |
| Blank element columns stop `extract` (`Element symbols are missing in '...'`); `all` fills them itself and stops when some atoms cannot be assigned | Run `add-elem-info` on the original PDB | {ref}`Input / extraction <input--extraction>` |
| `[multi] Atom count mismatch` / `Coordinate shape mismatch` / `Element sequence mismatch` | Regenerate all PDBs with the same preparation tool and settings, build the topology again from the structure you compute with, and never reorder atoms after `mm-parm` | {ref}`Input / extraction <input--extraction>`, {ref}`AmberTools / mm-parm <ambertools--mm-parm>` |
| **Charge & spin** | | |
| `ML-region charge is unresolved` / `[all] ML-region charge could not be resolved` | Set `-q/--charge` or `-l/--ligand-charge` explicitly | {ref}`Charge / spin <charge--spin>` |
| Energies or states look wrong after a run | Re-check the charge and multiplicity of the ML region | {ref}`Charge / spin <charge--spin>` |
| **Calculation & convergence** | | |
| CUDA out of memory (`torch.cuda.OutOfMemoryError`) | Check the Frozen-MM layer, shrink the ML region (`--radius`), narrow the Hessian (`--hessian-cutoff`), return to the default `FiniteDifference` if you selected `Analytical`, or move to a larger GPU | {ref}`CUDA OOM <cuda-oom-torchcudaoutofmemoryerror>` |
| TS optimization does not converge (`TS optimization did not converge`), or n_imag is not 1 after it | Check the TS candidate first, then switch the optimizer (`tsopt --opt-mode` / `all --opt-mode-post`); for n_imag ≥ 2, add `--flatten` | {ref}`TS optimization <troubleshooting-ts>`, {ref}`When the TS search fails <ts-search-fails>` |
| IRC does not terminate | Check the optimized endpoints first, then reduce the step: `irc --step-size` or `all --irc-step-size` | {ref}`IRC <troubleshooting-irc>` |
| Optimizer stalls at a flat energy (possible MLIP noise floor) | Let `--max-cycles` bound the run, or opt in to `--stop-plateau`; tune `--stop-plateau-thresh` / `--stop-plateau-window` if it stops too early or too late | {ref}`Plateau stop <optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor>` |
| **Installation & environment** | | |
| UMA model 401 / 403 or gated-repo error (`huggingface_hub.errors.GatedRepoError`) | Run `hf auth login` and accept the UMA model license | {ref}`Installation / environment <installation-environment-problems>` |
| `orb-models is required for the ORB backend` (or the same for AIMNet2 / MACE) | Install the backend extra: `pip install "mlmm-toolkit[orb]"` or `"mlmm-toolkit[aimnet]"`; MACE goes in a separate environment | {ref}`Backend-specific <troubleshooting-backends>` |
| `mm-parm` cannot run (`AmberTools preflight failed`; `tleap` / `antechamber` / `parmchk2` missing) | Make AmberTools available first | {ref}`AmberTools / mm-parm <ambertools--mm-parm>` |
| `hessian_ff` build or import errors (`hessian_ff build attempts failed`) | Check the C++20 compiler, then rebuild the native extension | {ref}`hessian_ff build / import <hessian_ff-build--import>` |
| DMF mode import error (`DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2`) | Install `cyipopt` (conda-forge) and `pydmf[torch]>=1.2` (PyPI) | {ref}`Installation / environment <installation-environment-problems>` |
| CUDA / GPU runtime mismatch | Check the GPU, the PyTorch build, and the driver together | {ref}`Installation / environment <installation-environment-problems>` |
| Plot export fails | Run `plotly_get_chrome -y` to install headless Chrome | {ref}`Installation / environment <installation-environment-problems>` |

## Preflight checklist

Before a long run, check that:

- `mlmm -h` runs and shows the CLI help.
- A Hugging Face login is set up on this machine for the default UMA backend.
- Your input PDB/mmCIF structures contain hydrogens and element symbols.
- When you pass several PDBs, they share the same atoms in the same order.
- `tleap`, `antechamber`, and `parmchk2` are on `$PATH`.
- The `hessian_ff` C++ extension builds on first use; see {ref}`hessian_ff build / import <hessian_ff-build--import>`.

---

(input--extraction)=
## Input / extraction

### `Element symbols are missing in '...'`

- **Symptom**: `extract` stops with ``Element symbols are missing in '...'. For PDB input, run `mlmm add-elem-info -i ... --overwrite`, or write a fixed PDB with `-o` and pass that file to extract; ...``. `all` fills blank element columns itself before extraction, and stops with the same message when some atoms cannot be assigned.
- **Cause**: many PDBs leave the element column (columns 77–78) blank, and `extract` needs the elements to type the atoms. mmCIF input must provide `_atom_site.type_symbol`.
- **Fix**: fill the column with `add-elem-info` and rerun with the new file. For each atom listed after `[add-elem-info] WARNING: Could not confidently assign N atoms; left unchanged.`, type its element symbol into columns 77–78 by hand, right-aligned.

  ```bash
  mlmm add-elem-info -i input.pdb -o input_with_elem.pdb
  ```

### `[multi] Atom count mismatch` / `[multi] Atom order mismatch`

- **Symptom**: a run with several inputs stops with `[multi] Atom count mismatch between input #1 and input #2: ...` or `[multi] Atom order mismatch between input #1 and input #2.`
- **Cause**: the structures were prepared with different tools or settings, or the atom order changed after re-protonation.
- **Fix**: regenerate **all** structures with the same protonation tool and settings. For MD snapshots, take every frame from the same trajectory and topology. If matching several inputs is not practical, start from one PDB and build the path with [`--scan-lists`](quickstart-scan.md).

### The ML region is too small or misses catalytic residues

- **Symptom**: the extracted ML region is smaller than expected, or catalytic residues are missing.
- **Cause**: the radius (`-r/--radius`, default 2.6 Å) is too small for this site.
- **Fix**: raise `--radius` (for example 2.6 → 3.5 Å), or add the residue with `--selected-resn 'A:TYR:44'`, which adds it without starting a distance search from it; adding it to `-c` keeps it whole when `-r` is above 0 ([Make the model larger](model-setup.md#make-the-model-larger)). The accepted forms are in {ref}`Residue selectors <selected-resn-takes-ids>`; in a PDB with an empty chain column, use the name or the number alone, such as `'44'`. You can also select the ML atoms yourself and pass the PDB with `--model-pdb` ([Use a model you built yourself](model-setup.md#use-a-model-you-built-yourself)).

### Energies or barriers shift with the size of the ML region

Enlarge the ML region as in [The ML region is too small or misses catalytic residues](#the-ml-region-is-too-small-or-misses-catalytic-residues), and check how the result depends on the size of the ML region and the position of its boundary.

### A modified residue is not truncated

- **Symptom**: `extract` prints `[extract] WARNING: Residue ... may be an amino acid (has N, CA, C, O) but is not recognized as a standard residue name. Backbone truncation was not applied. ...`, and the residue keeps its full backbone.
- **Cause**: backbone truncation and link-hydrogen placement need an entry in the residue table. Residues such as SEP, TPO, and MLY are already in it.
- **Fix**: register the residue with its integer charge, for example `--modified-residue "HD1:0"` (also accepted by `mlmm all`); the rules for a bare name are in [extract](extract.md#5-non-standard-residues---modified-residue). If the backbone topology is unusual, build the ML region by hand and pass `--parm7` and `--model-pdb` to the downstream commands.

---

(charge--spin)=
## Charge / spin

`-q` is the charge of the ML region, not of the whole system. Check that each residue key in `-l/--ligand-charge` exists in the structure; the rules are in {ref}`Charge specification <charge-specification>`.

### `ML-region charge is unresolved` / `ML-region charge could not be resolved`

- **Symptom**: an individual command stops with `ML-region charge is unresolved. Provide -q/--charge or --ligand-charge.`, or `all` stops with `[all] ML-region charge could not be resolved. Provide -q/--charge, --ligand-charge, or calc.model_charge in YAML.`
- **Cause**: without `-q/--charge`, the charge is the total of the standard residues, ions, and `-l/--ligand-charge` values in the ML region, and then YAML `calc.model_charge`. None of them applied. With `--model-indices`, the charge cannot be derived from `-l`.
- **Fix**: give the charge and multiplicity, or a per-residue charge map with extraction. The console prints the derived charge as `Total active site model charge`.

  ```bash
  mlmm path-search -i R.pdb P.pdb --parm7 real.parm7 --model-pdb model.pdb -q 0 -m 1
  mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
  ```

---

(ambertools--mm-parm)=
## AmberTools / `mm-parm`

### `AmberTools preflight failed`

- **Symptom**: `mm-parm` stops with `AmberTools preflight failed. Missing required command(s): ... Required: tleap, antechamber, parmchk2`, or `all` stops with `[preflight] Missing required command(s) for mm_parm (AmberTools): ...`.
- **Fix**: install AmberTools with `conda install -c conda-forge ambertools=24.8 "numpy>=2,<2.5" -y`, load it with `module load ambertools` on HPC, or build it from source (<https://ambermd.org/AmberTools.php>). Verify with `which tleap antechamber parmchk2`. Without AmberTools, the individual commands still run when you pass a topology built elsewhere with `--parm7`.

### `antechamber` fails for a ligand

- **Symptom**: `mm-parm` stops with `[<RES>] antechamber failed (see log).`, or with `[<RES>] electron-count check failed before antechamber: ...` before antechamber runs.
- **Cause**: the charge or multiplicity does not match the hydrogens of the ligand, or its element symbols, connectivity, or TER records are wrong.
- **Fix**:
  - Check the element symbols, hydrogens, connectivity, and TER records of the ligand.
  - Give the formal charge with `-l 'LIG:-1'`, and the multiplicity of a non-singlet ligand with `--ligand-mult 'HEM:1,NO:2'` (`--auto-mm-ligand-mult` in `all`).
  - Rerun with `--keep-temp` (`--auto-mm-keep-temp` in `all`) to keep the working directory `parm7build_*`, and read `<resname>.antechamber.log` there.
  - Run antechamber by hand on the ligand: `antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc -3 -at gaff2`.
  - For RESP charges or other custom parameters, build the topology yourself with tleap and your own `frcmod` / `lib` files, and pass it with `--parm7` (see [mm-parm](mm-parm.md#notes)).

### `Coordinate shape mismatch for '...': got (N, 3), expected (M, 3)`

- **Symptom**: a calculation stops with this message.
- **Cause**: the structure does not have the same atoms as the `parm7` topology.
- **Fix**: build the topology again with `mm-parm` from the structure you compute with, or compute with the PDB that `mm-parm` writes with `-o` or `--add-h`. Never reorder PDB atoms after `mm-parm`.

### `oniom-export`: `Element sequence mismatch at atom index ...`

- **Fix**: give `-i` the same PDB that the `parm7` was built from. `--no-element-check` turns the check off (then verify the result by hand); a different atom count still stops the export (see [oniom-export](oniom-export.md#notes)).

---

(hessian_ff-build--import)=
## `hessian_ff` build / import

- **Symptom**: `hessian_ff build attempts failed: ...`, or a calculation stops with `native bonded extension is unavailable.`, `native nonbonded extension is unavailable; torch fallback is disabled.`, or `analytical Hessian native extension is unavailable.`
- **Cause**: the C++ extension is compiled on first use through `torch.utils.cpp_extension`. It needs a C++20 compiler (validated with GCC 13.3) and `ninja`, which is installed with `mlmm-toolkit`.
- **Fix**:
  - Check C++20 support with `g++ -std=c++20 -x c++ -fsyntax-only /dev/null`. To install a compiler, run `conda install -c conda-forge gxx_linux-64`, or on HPC load the site's C++20 compiler module.
  - Check that the PyTorch headers are found: `python -c "import torch; print(torch.utils.cmake_prefix_path)"`.
  - The build uses a local temporary directory by default, because a network-mounted directory (NFS/Lustre) can hang on the build lock of PyTorch; to choose another local path, set `TORCH_EXTENSIONS_DIR`.
  - Check that `hessian_ff` is importable from the Python environment you run `mlmm` in.
  - For a manual clean rebuild:

  ```bash
  cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
  ```

---

## B-factor layer assignment

The layers are stored in the B-factors: ML = 0.0, Movable-MM = 10.0, Frozen-MM = 20.0 (tolerance ±1.0). `--detect-layer` (on by default) reads them.

### Wrong layer assignments / ML region too small or too large

- Open the layered PDB in a viewer and color it by B-factor.
- Check that `--model-pdb` selects the intended atoms.
- Adjust the Movable-MM / Frozen-MM boundary with `define-layer --movable-cutoff` (default 8.0 Å).
- The atoms in the Hessian are set separately, with `--hessian-cutoff` or YAML `calc.hess_cutoff` / `calc.hess_mm_atoms`.

### B-factors are not recognized

- **Symptom**: `all` stops with `[all] ... does not contain a valid ML/MM B-factor partition (both ML and MM atoms are required). ...` or `[all] Automatic layer detection requires a valid 0/10/20 B-factor partition with both ML and MM atoms when extraction is skipped and --model-pdb is absent.`
- **Cause**: the B-factors are read as layers only when at least one atom is ML (0), at least one is MM (10 or 20), and at least 80% of the atoms carry one of these values.
- **Fix**: run `define-layer` again and use the PDB it writes. Do not hand-edit B-factors to arbitrary values.

### `--detect-layer` does not work as expected

- **Symptom**: the automatic layers split the system differently from what you intended, or `all` without `-c` stops with `[all] Skipping extraction (no -c/--center) with B-factor layer detection disabled requires --model-pdb. ...`
- **Fix**:
  - Give a PDB input, or an XYZ input with `--ref-pdb`.
  - Assign the layers with `define-layer` and use the PDB it writes.
  - A `--movable-cutoff` given to a calculation command turns off `--detect-layer`: the MM layers then come from the distance, not from the B-factors.

---

(installation-environment-problems)=
## Installation / environment

First confirm that the optional packages are installed in the active environment and that PyTorch sees the GPU. After a repair, check with `mlmm --version` and `python -c "import torch; print(torch.cuda.is_available())"`, then run the command once with `--dry-run` to check the options and input before the full run.

| Symptom | Cause | Fix |
| --- | --- | --- |
| UMA download fails (`huggingface_hub.errors.GatedRepoError`, `401`, `403`) | No Hugging Face login, or the UMA model license is not accepted | Run `hf auth login` once per environment and machine, and accept the UMA model license on its Hugging Face page. On HPC, make sure compute nodes can write to the Hugging Face cache directory |
| `torch.cuda.is_available()` returns `False`, or a CUDA runtime error at import | The PyTorch build does not match the GPU or driver of the node | Check the assigned GPU, the installed wheel, and the driver with `nvidia-smi`, `python -m torch.utils.collect_env`, and `python -m pip check`. The `CUDA Version` that `nvidia-smi` shows is the newest CUDA the driver supports; install a PyTorch wheel at or below it (`cu126`, `cu130`, or `cu132`) |
| `--mep-mode dmf` fails with `DMF mode (--mep-mode dmf) requires ase, cyipopt, and pydmf>=1.2` | `cyipopt` and `pydmf` are not installed with `mlmm-toolkit` (`ase` is) | Install them as in {ref}`step 3 of the installation <step-by-step-installation>`: `conda install -c conda-forge cyipopt -y`, then `pip install 'pydmf[torch]>=1.2'` (for `--dmf-backend cpu` only: `pip install 'pydmf>=1.2'`) |
| Plot export fails (Plotly / Chrome) | Headless Chrome is missing | Run `plotly_get_chrome -y` once; it downloads a Chromium binary and needs internet access |

### DMF is unusually slow inside IPOPT

If IPOPT/MUMPS uses multithreaded BLIS, nested threads can cause long waits.
Set `BLIS_NUM_THREADS=1` before starting Python or the CLI, for example in
the job script. Leave the outer OpenMP/MM thread count unchanged.
Manual `BLIS_JC_NT`, `BLIS_PC_NT`, `BLIS_IC_NT`, `BLIS_JR_NT`, or
`BLIS_IR_NT` settings override this limit; remove them from that job's
configuration. Restart an existing notebook kernel after changing the settings.
See [BLIS thread controls](https://github.com/flame/blis/blob/2.0/docs/Multithreading.md).

---

(calculation--convergence)=
## Calculation / convergence

First check the TS candidate: a successful TS optimization gives one imaginary mode along the reaction coordinate. A mode counts as imaginary when ν < −5.00 cm⁻¹ (set by YAML `freq.zero_cutoff_cm`). Inspect its displacement and the IRC endpoints.

(cuda-oom-torchcudaoutofmemoryerror)=
### CUDA OOM (`torch.cuda.OutOfMemoryError`)

Try in order:

1. **Check Frozen-MM**: `define-layer` should put distal atoms at B = 20.0. If the Frozen-MM region is too small, the Movable-MM region and its Hessian grow. A smaller `--movable-cutoff` enlarges Frozen-MM ([Thin the movable MM shell](model-setup.md#thin-the-movable-mm-shell)).
2. **Shrink the ML region**: a smaller `--radius` in `extract`, or a smaller ML region given with `--model-pdb` ([Shrink the ML region](model-setup.md#shrink-the-ml-region)).
3. **Narrow the Hessian**: `--hessian-cutoff` of `opt`, `tsopt`, `freq`, and `sp` keeps fewer Movable-MM atoms in the Hessian ([Narrow the Hessian](model-setup.md#narrow-the-hessian)).
4. **Compare Hessian modes**: finite differences often lower the ML autograd memory, but both modes form a dense Hessian over the active atoms; compare runtime and peak memory on the target system, and return to the default `FiniteDifference` if you selected `Analytical`.
5. **Use a bigger GPU**: try the same model, Hessian mode, and active region on the target device first.

(optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor)=
### Optimizer "stalls" with flat energy + forces just above threshold (MLIP force noise floor)

- **Symptom**: `opt` / `tsopt` keeps running while the energy stays flat and the forces stay just above the `gau` / `baker` thresholds.
- **Cause**: MLIPs have finite numerical precision. For large ML/MM systems the noise floor can exceed the force thresholds, so the forces never drop further.
- **Fix**:
  - Runs stop at `--max-cycles` (default 100000). To stop earlier, add `--stop-plateau` (`opt`, `tsopt`, and `all`): the run stops as `stalled`, which is not converged.
  - To tune it, see [YAML Reference](yaml-reference.md#opt).
  - The check skips GSM / DMF. It applies to the single-structure preoptimization in `path-opt` / `path-search`.

(troubleshooting-ts)=
### TS optimization does not converge / multiple imaginary modes remain

- **Symptom**: the TS optimization runs many cycles without converging (`TS optimization did not converge. Review the TS trajectory.` in `summary.log`), or after it converges n_imag is 2 or more (`TS imaginary-mode validation found n_imag=N.`) or 0 (`[tsopt] No imaginary mode detected. Try all --refine-path.`); see [Reading the TS result](tsopt.md#reading-the-ts-result).
- **Fix when the optimization does not converge**: inspect the stop reason and the mode displacements, then try the following in order.
  1. Switch the optimizer between RS-P-RFO (the default) and the Dimer method: `tsopt --opt-mode hess` / `dimer`, or `all --opt-mode-post hess` / `grad` (Dimer).
  2. Reduce the step size in YAML; see [YAML Reference](yaml-reference.md#ts-optimization-sections).
  3. Start from another candidate, such as a better HEI (highest-energy image) of the path; see {ref}`When the TS search fails <ts-search-fails>`.
- **Fix when n_imag ≥ 2 remains**: re-optimize with `--flatten`, or tighten the convergence preset from the default `baker` to `gau_tight` or `gau_vtight` (`tsopt --thresh` or `all --thresh-post`). Other moves, such as a wider Hessian or `--refine-path`, are in {ref}`When the TS search fails <ts-search-fails>`.

(troubleshooting-irc)=
### IRC does not terminate properly

An IRC that stops before it converges is still usable when the endpoint optimizations reach the intended R and P, so check the optimized endpoints first ([Judging the IRC](irc.md#judging-the-irc)).

- **Symptom**: the IRC stops before it reaches a clear minimum, or the energy oscillates and the gradient stays large.
- **Cause**: the step is too large for this surface, or the starting structure has more than one imaginary mode.
- **Fix**:
  - Standalone `irc`: `--step-size 0.05` (default 0.10 bohr); see [example 4 of irc](irc.md#4-retry-with-a-smaller-step).
  - `all`: `--irc-step-size 0.05`.
  - Confirm that the starting structure has n_imag = 1.
  - To ignore the physical stop criteria and trace to the cycle limit, use `irc --never-stop` or `all --irc-never-stop`, then inspect the trajectory and the endpoints.

(troubleshooting-mep)=
### MEP search (GSM / DMF) fails or misses bonds

- **Symptom**: the minimum energy path (MEP) search ends without a usable path (`MEP optimization did not converge. Review the MEP trajectory and convergence log.`), or misses an expected bond change.
- **Fix**:
  - Raise `--max-nodes` (default 20) to 30 or 40 for complex reactions.
  - Keep endpoint preoptimization on (the default); remove `--no-preopt` if you passed it.
  - Try the other method: `--mep-mode dmf` ↔ `gsm`.
  - Tune bond-change detection with YAML `bond.bond_factor` and `bond.delta_fraction`.

---

(troubleshooting-performance)=
## Performance / stability tips

- **Out of memory**: follow {ref}`CUDA OOM <cuda-oom-torchcudaoutofmemoryerror>`, or lower `--max-nodes`. In `opt` and `scan`, keep {ref}`--opt-mode grad <opt-mode-semantics>` (L-BFGS, no Hessian).
- **Analytical ML Hessian**: it can reduce evaluations, but its memory use depends on the backend and the system; compare it with `FiniteDifference` on a pilot.
- **MM Hessian**: the default `mm_fd: true` (finite differences) trades speed for memory; `mm_fd: false` is faster on small systems but uses more memory.
- **Multi-GPU**: ML uses one device (`ml_cuda_idx: 0`). The default `hessian_ff` MM backend runs on the CPU; to put MM on another GPU, select `mm_backend: openmm` together with `mm_device: cuda` and `mm_cuda_idx: 1`.
- **ML/MM parallelism**: ML (GPU) and MM (CPU) run in parallel by default; set the CPU threads with `mm_threads`.

(troubleshooting-backends)=
## Backend-specific

### A backend package is missing

- **Symptom**: `orb-models is required for the ORB backend. ...`, `aimnet is required for the AIMNet2 backend. ...`, or `mace-torch is required for the MACE backend. ...`
- **Fix**:
  - ORB: `pip install "mlmm-toolkit[orb]"`. AIMNet2: `pip install "mlmm-toolkit[aimnet]"`.
  - MACE: use a dedicated environment, `pip uninstall -y fairchem-core && pip install mace-torch`; `mace-torch` pins `e3nn==0.4.4`, while UMA (`fairchem-core`) needs `e3nn>=0.5`.
  - If ORB still fails to import after installing the extra, run `python -m pip check` and fix the package that the resolver or the import error names; do not add unrelated PyG packages.

### CUDA OOM while building a Hessian

Use `--hessian-cutoff` or the `FiniteDifference` Hessian mode. YAML `ml_device: cpu` avoids the GPU memory limit at a higher runtime cost.

---

## How to report an issue

Include the exact command, `summary.log` (or the console output), the smallest reproducing inputs, your environment (OS / Python / CUDA / PyTorch), and whether AmberTools and `hessian_ff` are installed and working.

## See also

- [Tips for studying reaction mechanisms](mechanism-tips.md) — what to try when the TS search fails
- [Installation](installation.md) — environment setup and optional backends
- [MLIP Backends](backends.md) — choosing a backend
- [Building the ML region and layers](model-setup.md) — check, trim, or enlarge the ML region and the layers
