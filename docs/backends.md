# MLIP Backends

This page explains how to choose the backend that computes the ML region and
lists, for each backend, the install command, model names, precision,
reproducibility settings, and Hessian evaluation mode. The default backend is
**UMA** (Meta's Universal Models for Atoms); `-b/--backend` also selects
**ORB**, **MACE**, and **AIMNet2**. All four are machine-learning interatomic
potentials (MLIPs). Whichever backend you choose, the MM region is computed
from the Amber topology (`--parm7`) and the two are combined by ONIOM.

## Per-backend characteristics

Select a backend with `-b/--backend` on any ML/MM calculation command, or set
`calc.backend` in YAML:

```bash
# UMA (default)
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0

# ORB
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b orb

# MACE
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b mace

# AIMNet2
mlmm opt -i complex.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -b aimnet2
```

| backend | install | model identifier | precision option |
|---------|---------|------------------|------------------|
| `uma` | `pip install mlmm-toolkit` (`fairchem-core` is a core dependency) + [Hugging Face login](installation.md) | `uma-s-1p2` (default) / `uma-m-1p1` | `uma_precision="fp32" \| "fp64"` |
| `orb` | `pip install "mlmm-toolkit[orb]"` | `orb_v3_conservative_omol` | `orb_precision="float32-high" \| "float32-highest" \| "float64"` (`fp32` / `float32` are accepted as other names) |
| `mace` | dedicated env: `pip uninstall -y fairchem-core && pip install mace-torch` (`mace-torch` pins an `e3nn` version that conflicts with UMA, so UMA does not run in this env) | `MACE-OMOL-0` | `mace_dtype="float32" \| "float64"` |
| `aimnet2` | `pip install "mlmm-toolkit[aimnet]"` | `aimnet2` | n/a |

`--backend-model NAME` overrides the model variant for the selected `--backend`
(e.g. `--backend uma --backend-model uma-m-1p1`). `-b dft` computes the ML region with DFT
([DFT/MM backend](#dftmm-backend)), and `--calc-file` plugs in any ASE
calculator ({ref}`Custom backend <backends-custom-calculator>`).

The run prints the backend and model it loads, for example
`[backend] Preparing MLIP model (UMA / UMA-S-1.2 (OMol))...`, and the JSON
output records them as `mlip_backend`, `mlip_model`, and `mlip_precision`
([JSON output](json-output.md#common-envelope)).

### Precision

`--precision fp32|fp64` sets the floating-point precision of MLIP inference.
When `--precision` is not given, each backend takes its own default:

| backend | without `--precision` | `--precision fp64` |
|---------|-----------------------|--------------------|
| `uma` | fp32 | accepted |
| `orb` | fp64 | accepted |
| `mace` | fp64 | accepted |
| `aimnet2` | fp32 (no precision setting; `--precision fp32` changes nothing) | error |
| `custom` (`--calc-file`) | the calculator's own setting | error (`--precision fp32` is also an error) |

Which value to choose depends on the purpose:

| Purpose | Recommended | Why |
| --- | --- | --- |
| Routine run | Leave unset | Keeps the defaults above: UMA/AIMNet2 fp32, ORB/MACE fp64. |
| Speed screening | `--precision fp32` only when needed | This lowers ORB/MACE precision (see [Notes](#notes)). |
| Final TS/Hessian | Leave unset; with UMA, compare `--precision fp64` when n_imag ≥ 2 ([tsopt](tsopt.md#wrong-imaginary-mode-count-after-optimization)) | Whatever the precision, check n_imag from the final Hessian of `tsopt` and confirm with IRC and the endpoint optimizations that the TS connects the intended R and P. |

Enable fp64 with:

```bash
mlmm tsopt -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
mlmm freq -i opt.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
mlmm irc -i ts.pdb --parm7 real.parm7 --model-pdb ml_region.pdb -q 0 -m 1 --precision fp64
```

Or via YAML config:

```yaml
calc:
 precision: fp64
```

## Determinism and reproducibility

`--deterministic` makes repeated runs with the same input give the same result
on the same software and GPU. Without it, two GPU runs with identical inputs can
differ in the last digits.

```bash
mlmm opt -i complex.pdb --parm7 enzyme.parm7 --model-pdb ml_region.pdb -q 0 --deterministic
mlmm all -i r_complex.pdb p_complex.pdb -c PRE -q -1 --deterministic
```

| ML backend | `--deterministic` |
|---|---|
| `uma` | Supported; check that two runs match for the installed model version |
| `orb` / `mace` | PyTorch's deterministic mode is turned on; check that two runs match for the installed backend version |
| `aimnet2` | **Not supported**: the run stops with an error (see [Notes](#notes)) |
| `custom` (`--calc-file`) | **Not supported**: the run stops with an error, because the supplied calculator is outside mlmm-toolkit's control |

## Workers and Hessian mode

`--uma-workers N` (default 1) runs N parallel UMA predictors (this needs
`fairchem-core[extras]`), and `--uma-workers-per-node` (default 1) sets how
many of them run on each node. Both flags exist on `opt`, `tsopt`, `freq`,
`irc`, `sp`, `all`, `path-opt`, `path-search`, `scan`, `scan2d`, and `scan3d`.
The other backends ignore them with a warning. HPC job templates are in
[Device Configuration & HPC Setup](device-hpc.md).

### Hessian evaluation mode

`--hessian-calc-mode` chooses how the Hessian of the ML region is computed; it
exists on `freq`, `irc`, `tsopt`, `sp`, and `all`, and YAML uses
`calc.hessian_calc_mode`. `FiniteDifference` (default) takes central differences of the forces;
`Analytical` uses the autograd or native Hessian of the backend. UMA (with one
worker), ORB, MACE, and AIMNet2 compute analytical Hessians when the installed
version provides them, and the DFT backend does so without `--embedcharge`. A custom calculator supports only
`FiniteDifference`. When the analytical Hessian is not available for the
selected backend, the run stops with an error. The [MM part of the Hessian](mlmm-calc.md) is
set separately.

Choose one of these two settings with UMA:

```bash
--uma-workers 1 --hessian-calc-mode Analytical       # analytical Hessian
--uma-workers 4 --hessian-calc-mode FiniteDifference # parallel UMA predictor + FD
```

Model precision and Hessian precision are separate settings. The Hessian is
assembled in float64 by default (`calc.H_double: true`); `H_double: false`
returns it in float32. With `--precision fp64`, the Hessian is always float64:
a `H_double: false` in the config is overridden with a warning.

## xTB electrostatic correction

For MLIP/MM workflows, `--embedcharge/--no-embedcharge` (default off) adds
`E_xTB(ML + MM point charges) - E_xTB(ML)` and the corresponding force and
Hessian difference, so the ML region feels the MM point charges.

With `-b dft`, `--embedcharge` instead places the MM point charges directly in
the PySCF Hamiltonian, without the xTB correction.

## DFT/MM backend

The 11 ML/MM calculation commands (`all`, `opt`, `tsopt`, `irc`, `freq`,
`scan`, `scan2d`, `scan3d`, `path-opt`, `path-search`, and `sp`) accept
`-b dft --func-basis FUNCTIONAL/BASIS --dft-engine gpu|cpu` (defaults
`wb97m-v/def2-svp` and `gpu`), which computes the ML region with
PySCF/GPU4PySCF while the MM region stays on the Amber force field. The
separate `mlmm dft` command gives single points with population analysis.
Low-memory mode is described in [`dft`](dft.md#how-it-works); CPU threads,
host RAM, and SCF checkpoints are listed in the
[`all` reference](reference/commands/all.md).

(backends-custom-calculator)=
## Custom backend — bring your own ASE Calculator (`--calc-file`)

Beyond the built-in MLIP backends, the **ML region** can be computed by any
[ASE](https://wiki.fysik.dtu.dk/ase/) Calculator supplied at run time with
`--calc-file`, without modifying mlmm-toolkit. This couples the ML side of the
ML/MM ONIOM scheme to GFN-xTB (via `tblite` / `xtb-python`), DFTB+, ORCA,
Psi4, or any other engine with an ASE calculator — the boundary is the standard ASE
Calculator interface (energy in eV, forces in eV/Å).

Write a Python file exposing a `get_calculator` factory that returns an ASE
Calculator:

```python
# my_calc.py  (minimal illustrative example)
from ase.calculators.emt import EMT

def get_calculator(charge=0, spin=1, device="auto", **kwargs):
    return EMT()
```

Swap `EMT()` for the engine you want — e.g. `tblite.ase.TBLite(...)` for
GFN-xTB, the DFTB+ ASE calculator, or `ase.calculators.orca.ORCA(...)`. Then
pass the file to a stage or to `all`; it selects the `custom` ML backend and
overrides `--backend`:

```bash
mlmm sp    -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm opt   -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm freq  -i complex.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
mlmm all   -i R.pdb P.pdb --parm7 system.parm7 --model-pdb ml_region.pdb --calc-file my_calc.py -q 0 -m 1
```

- The factory receives `charge`, `spin` (multiplicity; also offered as `mult` /
  `multiplicity`), and `device` when its signature accepts them, or
  unconditionally if it declares `**kwargs`, so engines that need the total
  charge (e.g. xTB) can be configured. Use a different factory name with
  `--calc-file-func-name NAME`; a Calculator instance assigned to that name is
  also accepted.
- The custom calculator computes the **ML region only**; the MM side keeps its
  usual `hessian_ff` / OpenMM engine and the ONIOM coupling is unchanged.
  Hessians use finite differences, so `freq` and `tsopt --opt-mode hess` work
  with any engine. Frozen atoms are honored as usual.
- Available on `all` and every standalone ML/MM calculation command.
  `all` forwards the same factory to every stage that uses a calculator. For a
  permanent, installable backend with its own `--backend` name, see
  [For developers](#for-developers).

## Python API

In Python, select the backend with the `backend` argument of `MLMMCore` or of
the pysisyphus calculator `mlmm`. The classes, their parameters, and a runnable
example are in [ML/MM Calculator › Python API](mlmm-calc.md#python-api).

## For developers

### Backend dispatcher pattern

`MLMMCore` hands the ML region to the adapter of the selected backend and
keeps the MM calculation and the ONIOM coupling to itself. An unknown backend name raises
`ValueError`; mlmm-toolkit has no `auto` backend, and the workflows pass the
backend chosen on the command line.

### File map

| file | role |
|------|------|
| `mlmm/backends/__init__.py` | Turns `--precision`, `--backend-model`, `--calc-file`, and `--uma-workers` into the calculator settings of the selected backend |
| `mlmm/backends/mlmm_calc.py` | `MLMMCore` (ML/MM ONIOM coupling), `MLMMASECalculator` (ASE), `mlmm` (pysisyphus Calculator), the per-backend adapters, finite-difference Hessian assembly, and unit conversion |
| `mlmm/backends/pyscf_dft.py` | PySCF/GPU4PySCF DFT backend for the ML region, with electrostatic embedding; reuses the SCF state between steps |

To add a built-in backend with its own `--backend` name, follow recipe 3.2
"Add an MLIP backend" in
[CONTRIBUTING](https://github.com/t-0hmura/mlmm_toolkit/blob/main/CONTRIBUTING.md).

### GPU memory during an ML/MM stage

During an ML/MM stage, GPU memory is used by the selected ML backend and its
Hessian intermediates; topology handling and the analytical MM force field run
on the CPU. Standalone DFT is a separate stage. The finite-difference Hessian
loop in `mlmm/backends/mlmm_calc.py` evaluates one displacement direction at a
time to bound the number of simultaneous evaluations; re-run the GPU smoke
suite and check peak VRAM before changing it to a batched implementation. Stage
runners release calculators between stages, so later stages do not keep
earlier models in memory.

### ONIOM coupling vs raw MLIP

The MLIP adapters in `mlmm/backends/mlmm_calc.py` evaluate the **ML region
only**. The subtractive ONIOM energy formula (`# CHEMISTRY-RULE:1`), the
link-atom Hessian projection (`# CHEMISTRY-RULE:2`), and the three-layer partial
Hessian assembly (`# CHEMISTRY-RULE:8`) live in the same file. A new MLIP
backend does not need to know the ONIOM coupling; it only needs to return the
ML-region energy, forces, and Hessian in the correct units.

## Notes

- `--precision fp32` on ORB or MACE is for screening only; check n_imag before
  you use the result.
- For ORB, `--precision fp32` selects the reduced `float32-high` mode.
- AIMNet2 supports neither `--precision fp64` nor `--deterministic`; both stop
  the run with an error. AIMNet2 casts its model inputs to float32, and it
  computes forces with its own CUDA code outside PyTorch's deterministic mode.
  When you need repeatable runs, use UMA, ORB, or MACE with `--deterministic`
  and run twice in the same environment to compare.
- `--calc-file` accepts neither `--precision` (`fp32` or `fp64`) nor
  `--deterministic`; both stop the run with an error. Set the precision inside
  your calculator.
- `--deterministic` turns on PyTorch's deterministic algorithms (`torch.use_deterministic_algorithms`) and replaces one PyTorch operation that has no deterministic GPU version.
- `--deterministic` applies to the whole process: set on `all`, it covers every stage that `all` runs, so you do not pass it per stage.
- `--deterministic` can be slower: deterministic GPU operations may use different, slower implementations. Use it only when you need repeated runs to match.
- `--deterministic` stops with an error when PyTorch has no deterministic version of an operation in the run, instead of silently giving non-reproducible output.
- The environment variable `MLMM_STRICT_DETERMINISTIC=1` turns on the same mode for CI jobs or the Python API; with it set, `--no-deterministic` does not turn the mode off.
- `--deterministic` alone does not guarantee bit-identical results across
  machines or software versions; compare two runs on the target setup.
- With UMA, `--uma-workers` above 1 cannot be combined with
  `--hessian-calc-mode Analytical`: the run stops with an error because the
  parallel predictor has no autograd model. Use `--uma-workers 1` for an
  analytical Hessian, or `FiniteDifference` with several workers.
- MACE cannot share an environment with UMA; install it in its own conda env.
- A backend whose package is not installed stops with an error such as ``orb-models is required for the ORB backend. Install with `pip install orb-models`.``
- `--embedcharge` runs xTB with and without the MM point charges at every
  energy, force, and Hessian evaluation; keep the ML region to about 200–300
  atoms and benchmark the actual system first.

## See Also

- [ML/MM Calculator](mlmm-calc.md) — ONIOM coupling, MM Hessian, and the Python API (`MLMMCore`, `MLMMASECalculator`, `mlmm`).
- [Architecture](architecture.md) — directory map and dependency direction.
- [Device Configuration & HPC Setup](device-hpc.md) — GPU/CPU placement and job templates.
- [Refine an MLIP TS with DFT](dft-backend.md) — DFT settings, memory, and checkpoints.
- [Troubleshooting](troubleshooting.md) — detailed troubleshooting guide.
