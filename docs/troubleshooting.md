# Troubleshooting

Search this page for the error message you encounter. For a symptom-first index see [Common Error Recipes](recipes-common-errors.md).

## Preflight checklist

Before a long run, verify:

- `mlmm -h` runs and shows the CLI help.
- UMA weights can be downloaded (Hugging Face login is set up — see [Getting Started › Installation](getting-started.md#installation)).
- Your input PDB(s) contain hydrogens and element symbols.
- When you pass multiple PDBs, they share the same atoms in the same order.
- `tleap` is on `$PATH` (required by `mm-parm`).
- The `hessian_ff` C++ native extension built successfully (rebuild with `cd hessian_ff/native && make` if the auto-build failed).

---

(input--extraction)=
## Input / extraction

| Symptom | Cause | Fix |
|---|---|---|
| Missing element columns (cols 77–78) | Some design tools leave the column blank; `extract` needs element symbols for atom typing. `mlmm all` auto-runs `add-elem-info` as preflight. | `mlmm add-elem-info -i input.pdb -o input_with_elem.pdb`, then re-run. |
| `[multi] Atom count mismatch` / `[multi] Atom order mismatch` | Inputs prepared by different tools / settings. | Regenerate **all** structures with the same protonation tool + settings. For MD ensembles, extract frames from the same trajectory + topology. As a fallback, switch to the staged-scan workflow (one PDB + `--scan-lists`). |
| Pocket too small / catalytic residues missing | Default radius too small for the system. | Increase `--radius` (e.g. 2.6 → 3.5 Å); force-include residues with `--selected-resn 'A:123,B:456'`; or hand-craft an ML-region PDB in PyMOL and pass via `--model-pdb`. |
| Unreliable energies / barriers shifting with model size | Extracted pocket too small. | Increase `-r` (e.g. `mlmm extract -i complex.pdb -c 'SUB' -o pocket.pdb -r 4.0`). |
| Non-standard residues not truncated correctly (SEP / TPO / MLY / D-amino acids) | Backbone truncation + link-H placement only apply to known three-letter codes. | `--modified-residue "SEP,TPO,MLY"` (also accepted on `mlmm all`). If insufficient (unusual backbone topology), build the pocket manually and pass `--parm` + `--model-pdb` to downstream subcommands directly. |

---

(charge--spin)=
## Charge / spin

Per-stage calc subcommands require explicit `-q/--charge` (ML-region net charge) and `-m/--multiplicity`. In `mlmm all`, charge is resolved in order: `-q` override → extraction summary → `--ligand-charge` fallback (when extraction is skipped).

```bash
mlmm path-search -i R.pdb P.pdb --parm real.parm7 --model-pdb model.pdb -q 0 -m 1
mlmm all -i R.pdb P.pdb -c 'SAM,GPP' -l 'SAM:1,GPP:-3'
```

---

(ambertools--mm-parm)=
## AmberTools / `mm-parm`

| Symptom | Fix |
|---|---|
| `AmberTools preflight failed. Missing required command(s): tleap, antechamber, parmchk2` | `conda install -c conda-forge ambertools -y` (or `module load ambertools` on HPC, or build from source: <https://ambermd.org/AmberTools.php>). Verify with `which tleap`. Without AmberTools you can still run `opt` / `tsopt` / `path-search` if you supply `--parm` manually. |
| `antechamber` fails for a ligand | Check ligand element symbols + connectivity + TER records. Specify `-l 'LIG:-1'` and (for non-singlet) `--ligand-mult 'HEM:1,NO:2'`. Inspect `<resname>.antechamber.log` via `--keep-temp`. Try manually: `antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc -3 -at gaff2`. For higher-accuracy partial charges, generate RESP from HF/6-31G* and pass custom `frcmod` / `lib`. |
| `Atom count in parm7 does not match input PDB` / `parm7 topology does not match the input structure` / `Coordinate shape mismatch for... got (N,3), expected (M,3)` | Re-run `mm-parm` from the current PDB; use its output `<prefix>.pdb` for downstream subcommands (tleap may add / remove hydrogens). Never reorder PDB atoms after `mm-parm`. |
| `oniom-export` reports `Element sequence mismatch at atom index ...` | Use the same PDB for `-i` that was used to generate the parm7. As an escape hatch, `--no-element-check` disables the check (verify results manually). |

---

(hessian_ff-build--import)=
## `hessian_ff` build / import

```text
ImportError: cannot import name 'ForceFieldTorch' from 'hessian_ff'
RuntimeError: hessian_ff build attempts failed: ...
```

The C++ native extension is JIT-built on first use. If that fails: ensure `g++ >= 9` (`g++ --version`), that PyTorch headers are available (`python -c "import torch; print(torch.utils.cmake_prefix_path)"`), and that `ninja` is installed. On HPC: `module load gcc/11`. Then clean + rebuild:

```bash
conda install -c conda-forge ninja -y
cd $(python -c "import hessian_ff; print(hessian_ff.__path__[0])")/native && make clean && make
```

Also ensure `hessian_ff` is importable at all (it is if you installed mlmm-toolkit with `pip install -e .`).

---

## B-factor layer assignment

Encoding: ML = 0.0, Movable-MM = 10.0, Frozen-MM = 20.0 (tolerance ±1.0). Common symptoms:

- **Wrong layer assignments / ML region too small or too large** — verify `--model-pdb` selects the intended atoms; adjust `--radius-freeze` (default 8.0 Å) for the Movable / Frozen boundary; control Hessian-target MM separately via `hess_cutoff` / `hess_mm_atoms`. Inspect the layered PDB visually (color by B-factor).
- **B-factors not recognised** (calculator treats all atoms as one layer) — re-run `define-layer`; do not hand-edit B-factors to arbitrary values.
- **`--detect-layer` produces unexpected splits or fails without `--model-pdb`** — supply a PDB input (or XYZ + `--ref-pdb`); re-run `define-layer` explicitly; for distance-based control, set `hess_cutoff` / `movable_cutoff` and use `--no-detect-layer` (supplying `--movable-cutoff` already disables `--detect-layer`).

---

(installation--environment)=
## Installation / environment

| Symptom | Fix |
|---|---|
| MLIP model download fails / HF auth missing (`huggingface_hub.errors.GatedRepoError`, `401`, `403`) | `hf auth login` once per env / machine; accept the model licence on the HF page. On HPC, ensure HF cache dir is writable from compute nodes. |
| `torch.cuda.is_available()` returns `False` | Install PyTorch matching your cluster CUDA runtime; verify `nvidia-smi` and `python -c "import torch; print(torch.version.cuda, torch.cuda.is_available())"`. |
| `DMF mode requires ase, cyipopt, and pydmf to be installed.` | `conda install -c conda-forge ase cyipopt -y && pip install pydmf`. |
| Plot export fails (Plotly / Chrome) | `plotly_get_chrome -y`. |

(dmf-mode-fails-cyipopt--pydmf--ase-missing)=
### DMF mode fails (cyipopt / pydmf / ase missing)

See the `DMF mode requires ase, cyipopt, and pydmf to be installed.` row above.

(cuda--pytorch-mismatch)=
### CUDA / PyTorch mismatch

See the `torch.cuda.is_available()` row above.

(plot-export-fails-chrome-missing)=
### Plot export fails (Chrome missing)

See the `Plot export fails (Plotly / Chrome)` row above.

---

(calculation--convergence)=
## Calculation / convergence

(cuda-oom-torchcudaoutofmemoryerror)=
### CUDA OOM (`torch.cuda.OutOfMemoryError`)

ML/MM systems are larger than pure MLIP, so VRAM pressure is higher. Try in order:

1. **Verify Frozen-MM** — `define-layer` should put distal atoms at B=20.0. If the Frozen region is too small, the Movable-MM region (and its Hessian) inflates. Decrease `--radius-freeze` to expand Frozen.
2. **Shrink ML region** — smaller `--radius` in `extract`, or hand-craft a smaller `--model-pdb`.
3. **`--hessian-calc-mode FiniteDifference`** — slower but lower peak VRAM.
4. **Pre-define layers** with `define-layer` and `use_bfactor_layers: true` in YAML.
5. **Bigger GPU** — 24 GB+ for 500+ ML atoms; 48 GB+ for 1000+.

### TS optimisation does not converge / multiple imaginary modes remain

Try `--opt-mode grad` (Dimer) ↔ `--opt-mode hess` (RS-I-RFO); `--flatten` to flatten extra imaginary modes; `--max-cycles 20000`; tighter `--thresh baker` / `gau_tight`; expand Hessian-target atoms via `hess_cutoff`.

(optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor)=
### Optimizer "stalls" with flat energy + forces just above threshold (MLIP force noise floor)

MLIPs have a finite numerical precision. For large ML/MM systems the noise floor can exceed `gau` / `baker` gradient thresholds, so forces never drop further even though the geometry is stationary. **This is handled automatically** via `energy_plateau: true` (declares convergence when the 50-step energy range falls below 1.0e-4 au ≈ 0.06 kcal/mol). To tighten or disable:

```yaml
opt:
  energy_plateau_thresh: 1.0e-05   # stricter (au)
  energy_plateau_window: 100        # require a longer flat stretch
  # energy_plateau: false           # disable entirely (e.g. for benchmarking)
```

The plateau check is skipped automatically for chain-of-states optimisers (GS / DMF), so `path-opt` / `path-search` are unaffected.

### IRC does not terminate properly

Reduce `--step-size 0.05` (default 0.10); raise `--max-cycles 200`; verify the TS candidate has exactly one imaginary frequency before launching IRC.

### MEP search (GSM / DMF) fails or misses bonds

Raise `--max-nodes` (e.g. 15–20) for complex reactions; enable `--preopt`; try the alternate method (`--mep-mode dmf` ↔ `gsm`); tune bond-detection in YAML (`bond.bond_factor`, `bond.delta_fraction`).

---

## Performance / stability tips

- **OOM** — shrink ML region, shrink Hessian-target MM, lower `--max-nodes`, or use `--opt-mode grad`.
- **Analytical ML Hessian** is fastest when VRAM allows (24 GB+ recommended for 300+ ML atoms); else `FiniteDifference`.
- **MM Hessian** — default `mm_fd: true` (finite-difference) trades speed for memory; `mm_fd: false` is faster on small systems but heavier on memory. Cap MM atom count with `hess_cutoff`.
- **Large systems (2000+ atoms)** — make sure the Frozen layer is generous (`define-layer` with appropriate cutoffs) to keep the movable DOF count down.
- **Multi-GPU** — ML on one device (`ml_cuda_idx: 0`), MM on another (`mm_device: cuda`, `mm_cuda_idx: 1`).
- **ML/MM parallelism** — ML (GPU) and MM (CPU) run in parallel by default; tune CPU threads with `mm_threads`.

## Backend-specific

| Symptom | Fix |
|---|---|
| `ImportError: orb-models is required for the ORB backend` (or similar for AIMNet2 / MACE) | `pip install "mlmm-toolkit[orb]"` / `"[aimnet]"` / `pip install --no-deps mace-torch` (MACE in a dedicated env). |
| CUDA OOM on ORB / MACE / AIMNet2 | These backends use FD Hessians (more VRAM). Lower `hess_cutoff`, or set `ml_device: cpu` (slow but avoids the VRAM limit). |
| `XTBEmbedError: xTB command not found` | `conda install -c conda-forge xtb -y` and ensure `xtb` is on `$PATH`. Custom binary: set `xtb_cmd` in YAML. |

---

## How to report an issue

Include: the exact command line, `summary.log` (or console output), the smallest reproducing inputs, your env (OS / Python / CUDA / PyTorch), and whether AmberTools + `hessian_ff` are properly installed.
