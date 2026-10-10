# `mlmm freq`

Builds and diagonalizes the Hessian of the ML/MM model, writes the frequencies
and per-mode displacements, and computes QRRHO thermochemistry. Run
`mlmm freq -i <structure> --parm7 <real.parm7> -q <charge> --out-json`.
Success is the n_imag you expect: 0 for a minimum, 1 for a TS candidate.

## When to use

- Check that an optimized structure is a minimum (n_imag = 0) or a TS
  candidate (n_imag = 1).
- Free energies at 298.15 K and 1 atm by default, or at your own temperature
  and pressure.
- With frozen atoms, the partial-Hessian analysis (PHVA) is used
  automatically.

## Minimal run

```bash
mlmm freq -i ts.xyz --parm7 real.parm7 --ref-pdb full_enzyme.pdb \
    -q 0 -m 1 -b uma --out-json -o result_freq
```

At a higher temperature for the activation enthalpy:

```bash
mlmm freq -i ts.xyz --parm7 real.parm7 --ref-pdb full_enzyme.pdb \
    -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

## Judge success

The console thermochemistry summary prints `Number of Imaginary Freq = N`, and
`result.json` records it as `n_imaginary`. A minimum has n_imag = 0. A
successful TS optimization gives one imaginary mode along the reaction
coordinate, so a TS candidate has exactly 1. Residual imaginary modes in R or P
do not block thermochemistry. `freq` does not judge n_imag itself:
`scientific_status` is `success` whatever n_imag is.

`freq` retains every signed physical mode. The default imaginary criterion is
ν < −5.00 cm⁻¹, and YAML `freq.zero_cutoff_cm` sets another cutoff magnitude.
Positive modes between 0 and 5 cm⁻¹ remain in thermochemistry.
`n_negative_modes` counts every negative value. Raw negative counts are diagnostic and do not add a failure gate.

```
result_freq/
├── frequencies_cm-1.txt                  # all modes, sorted, cm⁻¹
├── mode_NNNN_<±freq>cm-1_trj.xyz / .pdb  # per-mode displacement
├── thermoanalysis.yaml                   # with --dump
└── result.json                           # with --out-json
```

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"])                 # minimum: 0; TS: 1
print(d["frequencies_cm"][:5])          # first five frequencies
t = d["thermochemistry"]
print(t["zpe_ha"], t["thermal_correction_energy_ha"], t["S_cal_per_mol_K"])
print(t["electronic_energy_ha"], "+",
      t["thermal_correction_free_energy_ha"], "=",
      t["sum_EE_and_thermal_free_energy_ha"])  # E + G_corr = G
print(t["symmetry_number"], t["symmetry_number_source"])
```

`--dump-hess result_freq/hessian.npy` writes the Hessian at that exact path;
a relative path is resolved from the current working directory, not under
`--out-dir`. The file is one `numpy.save` array: the Cartesian Hessian in
Hartree/bohr², not mass-weighted, atoms in input order, 3N×3N or only the
atoms selected by `--active-dof-mode`. `--read-hess` checks only size,
symmetry, and finiteness, so pass a Hessian computed for the same geometry,
charge, multiplicity, layers, and calculator.

## Thermochemistry

The QRRHO (Grimme) treatment uses a fixed 100 cm⁻¹ rotor cutoff: vibrations
below it are interpolated toward the free-rotor entropy limit, and higher ones
use the harmonic oscillator. The free energy is E + G_corr = G.

What you choose:

- `--temperature` (default 298.15 K) and `--pressure` (default 1.0 atm).
- The point group and the external rotational symmetry number are detected
  from each structure, and the 1/σ correction is always included. YAML
  `thermo.symmetry_number` is an advanced override.

## PHVA

When the input has frozen atoms, `freq` builds and diagonalizes only the
mobile-atom block. This is much cheaper for large systems. With three or more
frozen atoms not on one line, the usual ML/MM case, no rigid motion is removed
and every vibration of the mobile atoms is kept.

What you choose: the frozen atoms, assigned by `define-layer` or explicitly
through `geom.freeze_atoms`. `extract` only writes the capped pocket
structure.

## Pitfalls and recovery

- A small-magnitude imaginary frequency may be numerical or a real shallow
  mode. Inspect its displacement and repeat the Hessian at suitable precision;
  the QRRHO cutoff does not validate a stationary point.
- `--hessian-calc-mode FiniteDifference` often lowers peak model/autograd
  memory. Runtime depends on backend, model, system, precision, and hardware;
  benchmark the actual calculation, and remember that both modes build a dense
  active-space Hessian.
- An explicit analytical Hessian with `--uma-workers` above 1 is a hard error.
  Use one worker for analytical curvature, or select `FiniteDifference` before
  enabling the UMA parallel predictor.
- Thermochemistry depends on charge and spin; make sure `-q`/`-m` are correct
  or the ZPE will be off.

## Next step

- Usual upstream stages: [tsopt.md](tsopt.md), [irc.md](irc.md).
- The `--hessian-calc-mode` setting for UMA:
  [backends.md](../mlmm-install/backends.md#uma).
