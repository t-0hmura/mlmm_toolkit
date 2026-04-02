# Common Error Recipes

Use this page when you know the symptom but do not know which subcommand page to open first.
For full details, keep [Troubleshooting](troubleshooting.md) open in parallel.

## Quick routing

| Symptom | Start here | Then read |
| --- | --- | --- |
| Missing element columns / extraction aborts | `add-elem-info` on the original PDB | [Input / extraction](troubleshooting.md#input--extraction-problems) |
| "Charge is required" errors | Set `-q/--charge` and `-m/--multiplicity` explicitly | [Charge / spin](troubleshooting.md#charge--spin-problems) |
| Energies/states look wrong after a run | Re-check charge/multiplicity policy in CLI conventions | [Charge / spin](troubleshooting.md#charge--spin-problems) |
| `mm-parm` cannot run (`tleap`/`antechamber`/`parmchk2` missing) | Fix AmberTools availability first | [AmberTools / mm-parm](troubleshooting.md#ambertools--mm-parm-problems) |
| `hessian_ff` import/build errors | Rebuild native extension (`hessian_ff/native`) | [hessian_ff build](troubleshooting.md#hessian_ff-build-problems) |
| DMF mode import errors (`cyipopt`) | Install `cyipopt` in the active environment | [DMF mode](troubleshooting.md#dmf-mode-fails-cyipopt-missing) |
| TSOPT/IRC does not converge | Reduce step length (trust_radius for RFO/RS-I-RFO, max_step for L-BFGS), increase cycles, validate TS quality first | [Convergence](troubleshooting.md#calculation--convergence-problems) |
| CUDA/GPU runtime mismatch | Verify `torch.cuda.is_available()` and CUDA build pairing | [CUDA / PyTorch](troubleshooting.md#cuda--pytorch-mismatch) |
| Plot export failures | Install Chrome runtime for Plotly export | [Plot export](troubleshooting.md#plot-export-fails-chrome-missing) |

## Recipe 1: Extraction fails before MEP starts

Signal:
 - Errors mention missing element symbols, atom-count mismatch, or empty pockets.
First checks:
 - Confirm all inputs are prepared by the same workflow and atom ordering is consistent.
 - Ensure element columns are present before running `extract` or `all`.
Typical fix path:
 - Repair elements -> rerun extraction -> confirm pocket size and residue inclusion.

## Recipe 2: Charge/spin validation fails

Signal:
 - If the log shows unexpected charges (e.g., protein charge is wrong, or total charge does not match expectations), review the charge resolution rules.
First checks:
 - Ensure total charge and multiplicity are physically correct for the target state.
 - If using residue maps, validate each residue key in `--ligand-charge`.
 - Verify the resolution rules in [CLI Conventions](cli-conventions.md) when results look physically inconsistent.
Typical fix path:
 - Prefer explicit `-q` and `-m` for critical runs, then retry scan/path/tsopt.

## Recipe 3: Build or environment blockers

Signal:
 - `mm-parm` tooling not found, `hessian_ff` import failures, CUDA mismatch.
First checks:
 - Confirm required executables and Python extension modules exist in the active env.
 - Validate GPU visibility and PyTorch CUDA compatibility.
Typical fix path:
 - Repair toolchain/build first, then rerun with `--help-advanced` to verify available options before full execution.

## Recipe 4: Convergence and post-processing failures

Signal:
 - TSOPT stalls, IRC branches look unstable, or MEP refinement stops unexpectedly.
First checks:
 - Confirm TS candidate quality with one dominant imaginary mode.
 - Reduce step length (trust_radius / max_step) and increase cycle limits.
Typical fix path:
 - Run a smaller diagnostic case, tune thresholds/step sizes, then scale back up.
