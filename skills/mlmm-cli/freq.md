# `mlmm freq`

## Purpose

Vibrational analysis: build the Hessian, diagonalize for normal-mode
frequencies, write per-mode geometry displacements, and compute
QRRHO thermochemistry. Default temperature 298.15 K, 1 atm.
Partial-Hessian variant (PHVA) activates automatically when
`freeze_atoms` is non-empty.

## Synopsis

```bash
mlmm freq -i geom.{pdb,xyz} --parm real.parm7 \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--temperature 298.15] [--pressure 1.0] [--symmetry-number 1] \
    [-b uma|orb|mace|aimnet2] [-o ./result_freq/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Derive layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); used only when `--model-pdb` is omitted (`--model-pdb` takes precedence) |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | Unavailable in v0.3.3; use `--no-embedcharge` |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--temperature` | float | 298.15 | K, for thermochemistry |
| `--pressure` | float | 1.0 | atm, for thermochemistry |
| `--symmetry-number` | int ≥ 1 | 1 | External rotational symmetry number. Point-group symmetry is not inferred. |
| `--hessian-calc-mode` | str | `FiniteDifference` | `Analytical` / `FiniteDifference`; check `FREQ_KW` / `MLMM_CALC_KW` |
| `--tr-projection` | str | `constrained` | Frozen-boundary TR treatment. `legacy-active` is deprecated comparison-only behavior; never use it for pass/HOSP transition-state certification. |
| `--precision` | str | backend-specific | UMA/AIMNet2 fp32; ORB/MACE fp64; AIMNet2 rejects fp64 |
| `--workers` | int | 1 | UMA predictor workers; `>1` requires `fairchem-core[extras]` and is incompatible with `Analytical` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default (298.15 K, 1 atm)

```bash
mlmm freq -i ts.xyz --parm real.parm7 --ref-pdb full_enzyme.pdb \
    -q 0 -m 1 -b uma -o result_freq
```

### Higher temperature for activation enthalpy

```bash
mlmm freq -i ts.xyz --parm real.parm7 --ref-pdb full_enzyme.pdb \
    -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

## Output

```
result_freq/
├── result.json                          # when --out-json
├── frequencies_cm-1.txt                 # all modes, sorted, cm⁻¹
├── thermoanalysis.yaml                  # when --dump (ZPE, S, H, G)
└── mode_NNNN_<±freq>cm-1_trj.xyz / .pdb # per-mode displacement (visualize in PyMOL)
```

`--dump-hess result_freq/hessian.npz` writes the Hessian at that exact path;
a relative path is resolved from the current working directory, not relocated
under `--out-dir`. The schema-2 NPZ stores atom order, Cartesian geometry,
active-DOF basis, PHVA metadata, model charge, and multiplicity with the
Hessian. `mlmm irc --read-hess` accepts it only for the matching geometry,
layer/Hessian settings, and electronic state. Schema-1 files lack charge/spin
identity and require IRC's explicit `--allow-unverified-hess-state` opt-in;
unidentified legacy NPZ files are rejected.

`result.json` keys:

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"])                 # 0 for minimum, 1 for TS
print(d["frequencies_cm"][:5])          # first five frequencies
print(d["thermochemistry"]["zpe_ha"])
print(d["thermochemistry"]["thermal_correction_energy_ha"])
print(d["thermochemistry"]["S_cal_per_mol_K"])
t = d["thermochemistry"]
print(t["electronic_energy_ha"], "+",
      t["thermal_correction_free_energy_ha"], "=",
      t["sum_EE_and_thermal_free_energy_ha"])  # E + G_corr = G
print(d["thermochemistry"]["symmetry_number"],
      d["thermochemistry"]["symmetry_number_source"])
print(d["rigid_projection"]["treatment"], d["rigid_projection"]["effective_rank"])
```

## QRRHO thermochemistry

Default thermochemistry uses the QRRHO (Grimme) treatment with a
100 cm⁻¹ rotor cutoff:

- low-frequency vibrations (< 100 cm⁻¹) are interpolated toward the
  free-rotor entropy limit,
- high-frequency vibrations use the standard harmonic-oscillator
  partition function.

The QRRHO/rotor cutoff (100 cm⁻¹) is fixed by the vendored
thermoanalysis default and is not user-tunable via `THERMO_KW`.
`mlmm.core.defaults.THERMO_KW` exposes `temperature`, `pressure_atm`,
`symmetry_number`, and `dump`. Supply the external rotational symmetry number
when it is not 1; the workflow does not infer a point group.

## Partial-Hessian Vibrational Analysis (PHVA)

When the input has frozen atoms (PDB B-factor or `freeze_atoms`
set), `freq` automatically computes the **partial Hessian**: only the
mobile-atom block is built and diagonalized; frozen atoms are projected
out. This is much cheaper for large clusters.

Frozen atoms are written by `extract` for link-H parents. To override,
use `--config` YAML and set `freeze_atoms`.

The default `constrained` TR treatment removes only full-system rigid motions
that leave every frozen anchor fixed. Generic effective ranks are 6/3/1/0 for
zero/one/two/at least three non-collinear anchors; realistic ML/MM boundaries
normally have rank 0. All-frozen input is an explicit error. `legacy-active`
is an isolated-active comparison treatment using the current common kernel;
bitwise identity is not guaranteed for rank-degenerate cases. It is deprecated
and must not be used for pass/HOSP transition-state certification. `result.json` and dumped
`thermoanalysis.yaml` record the treatment, effective rank, Hessian source,
and Hessian shape under `rigid_projection`.

## Caveats

- A minimum should have **0 imaginary frequencies**, a TS should have
  **exactly 1**.
- A small-magnitude imaginary frequency may be numerical or a real shallow
  mode. Inspect its displacement and repeat the Hessian at suitable precision;
  the QRRHO cutoff does not validate a stationary point.
- `--hessian-calc-mode FiniteDifference` often lowers peak model/autograd
  memory. Runtime depends on backend, model, system, precision, and hardware;
  benchmark the actual calculation, and remember that both paths materialize a
  dense active-space Hessian.
- An explicit analytical Hessian with `workers > 1` is a hard error. Use one
  worker for analytical curvature or select `FiniteDifference` before enabling
  the UMA parallel predictor.
- Thermochemistry depends on charge / spin — make sure `-q`/`-m` are
  correct or ZPE will be off.

## See also

- `tsopt.md`, `irc.md` — usual upstream stages.
- `mlmm-install-backends/uma.md` — `--hessian-calc-mode` knob.
- Defaults: `import mlmm.core.defaults as d; print(d.FREQ_KW, d.THERMO_KW, d.MLMM_CALC_KW)`
