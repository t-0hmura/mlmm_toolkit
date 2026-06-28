# Common Error Recipes

Use this page when you know the symptom but do not know which subcommand page to open first.
For full details, keep [Troubleshooting](troubleshooting.md) open in parallel.

## Quick routing

| Symptom | Start here | Then read |
| --- | --- | --- |
| **Input & extraction** | | |
| Missing element columns / extraction aborts | `add-elem-info` on the original PDB | [Input / extraction](troubleshooting.md#input--extraction) |
| `[multi] Atom count mismatch` / `Coordinate shape mismatch` / `Element sequence mismatch` | Regenerate all PDBs with the same prep tool + settings; re-run `mm-parm` from the current PDB; never reorder atoms after `mm-parm` | [Input / extraction](troubleshooting.md#input--extraction) |
| **Charge & spin** | | |
| "Charge is required" errors | Set `-q/--charge` and `-m/--multiplicity` explicitly | [Charge / spin](troubleshooting.md#charge--spin) |
| Energies/states look wrong after a run | Re-check charge/multiplicity policy in CLI conventions | [Charge / spin](troubleshooting.md#charge--spin) |
| **Installation & environment** | | |
| UMA model 401/403 / gated-repo error (`huggingface_hub.errors.GatedRepoError`) | `hf auth login` and accept the model license | [Installation / environment](troubleshooting.md#installation--environment) |
| `ImportError: orb-models is required` (or similar for AIMNet2 / MACE) | `pip install "mlmm-toolkit[orb]"` / `"[aimnet]"`; MACE installs into a separate env | [Installation / environment](troubleshooting.md#installation--environment) |
| `mm-parm` cannot run (`tleap`/`antechamber`/`parmchk2` missing) | Fix AmberTools availability first | [AmberTools / mm-parm](troubleshooting.md#ambertools--mm-parm) |
| `hessian_ff` import/build errors | Rebuild native extension (`hessian_ff/native`) | [hessian_ff build](troubleshooting.md#hessian_ff-build--import) |
| DMF mode import errors (`ase` / `cyipopt` / `pydmf`) | Install `ase` + `cyipopt` (conda-forge) and `pydmf>=1.2` (PyPI) | [DMF mode](troubleshooting.md#dmf-mode-fails-cyipopt--pydmf--ase-missing) |
| **GPU & CUDA** | | |
| CUDA out-of-memory at runtime (`torch.cuda.OutOfMemoryError`) | Shrink ML region (`--radius`), use `--hessian-calc-mode FiniteDifference`, or move to a larger GPU | [CUDA OOM](troubleshooting.md#cuda-oom-torchcudaoutofmemoryerror) |
| CUDA/GPU runtime mismatch | Verify `torch.cuda.is_available()` and CUDA build pairing | [CUDA / PyTorch](troubleshooting.md#cuda--pytorch-mismatch) |
| **Convergence** | | |
| TSOPT/IRC does not converge | Reduce step length (trust_radius for RFO/RS-I-RFO, max_step for L-BFGS) | [Convergence](troubleshooting.md#calculation--convergence) |
| TSOPT/IRC does not converge | Increase cycles | [Convergence](troubleshooting.md#calculation--convergence) |
| TSOPT/IRC does not converge | Validate TS quality first | [Convergence](troubleshooting.md#calculation--convergence) |
| Optimizer stalls at flat energy (MLIP noise floor) | Rely on the default `energy_plateau` fallback | [Plateau fallback](troubleshooting.md#optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor) |
| Optimizer stalls at flat energy (MLIP noise floor) | Tune `energy_plateau_thresh` / `energy_plateau_window` if the trigger fires too early or too late | [Plateau fallback](troubleshooting.md#optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor) |
| **Plotting** | | |
| Plot export failures | Install Chrome runtime for Plotly export | [Plot export](troubleshooting.md#plot-export-fails-chrome-missing) |

## Recipe 1: Extraction fails before MEP starts

**Signal:**

- Errors mention missing element symbols, atom-count mismatch, or empty pockets.

**First checks:**

- Confirm all inputs are prepared with the same workflow and atom ordering is consistent.
- Ensure element columns are present before running `extract` or `all`.

**Typical fix path:**

- Repair elements -> rerun extraction -> confirm pocket size and residue inclusion.

## Recipe 2: Charge/spin validation fails

**Signal:**

- If the log shows unexpected charges (e.g., protein charge is wrong, or total charge does not match expectations), review the charge resolution rules.

**First checks:**

- Ensure total charge and multiplicity are physically correct for the target state.
- If using residue maps, validate each residue key in `--ligand-charge`.
- Verify the resolution rules in [CLI Conventions](cli-conventions.md) when results look physically inconsistent.

**Typical fix path:**

- Prefer explicit `-q` and `-m` for critical runs, then retry scan/path/tsopt.

## Recipe 3: Build or environment blockers

**Signal:**

- `mm-parm` tooling not found, `hessian_ff` import failures, CUDA mismatch.

**First checks:**

- Confirm required executables and Python extension modules exist in the active env.
- Validate GPU visibility and PyTorch CUDA compatibility.

**Typical fix path:**

- Repair toolchain/build first, then rerun with `--help-advanced` to verify available options before full execution.

## Recipe 4: Convergence and post-processing failures

**Signal:**

- TSOPT stalls, IRC branches look unstable, or MEP refinement stops unexpectedly.

**First checks:**

- Confirm TS candidate quality with one dominant imaginary mode.
- Reduce step length (trust_radius / max_step) and increase cycle limits. Note that `trust_max` defaults to 0.10 bohr for both RFO and RS-I-RFO.
- Check whether the energy has already plateaued. If the last ~50 cycles show `|dE| < 1e-4` au (atomic units) while forces are flat, the cause is the noise floor of the machine-learned interatomic potential (MLIP) force rather than an optimization bug. In that case the default `energy_plateau` fallback declares convergence automatically (see [Troubleshooting](troubleshooting.md#optimizer-stalls-with-flat-energy--forces-just-above-threshold-mlip-force-noise-floor)).

**Typical fix path:**

- Run a smaller diagnostic case, tune thresholds/step sizes, then scale back up.
