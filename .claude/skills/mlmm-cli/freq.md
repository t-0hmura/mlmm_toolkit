# `mlmm freq`

## Purpose

Vibrational analysis: build the Hessian, diagonalize for normal-mode
frequencies, write per-mode geometry displacements, and compute
QRRHO thermochemistry. Default temperature 298.15 K, 1 atm.
Partial-Hessian variant (PHVA) activates automatically when
`freeze_atoms` is non-empty.

## Synopsis

```bash
mlmm freq -i geom.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--temperature 298.15] [--pressure 1.0] \
    [-b uma|orb|mace|aimnet2] [-o ./result_freq/]
```

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Input geometry |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--temperature` | float | 298.15 | K, for thermochemistry |
| `--pressure` | float | 1.0 | atm, for thermochemistry |
| `--hessian-calc-mode` | str | (live default) | `Analytical` / `FiniteDifference`; check `UMA_CALC_KW` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_freq/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default (298.15 K, 1 atm)

```bash
mlmm freq -i ts.xyz -q 0 -m 1 -b uma -o result_freq
```

### Higher temperature for activation enthalpy

```bash
mlmm freq -i ts.xyz -l 'SAM:1' \
    --temperature 310.15 --pressure 1.0 \
    -b uma -o result_freq_310K
```

## Output

```
result_freq/
├── result.json
├── frequencies_cm-1.txt           # all modes, sorted, cm⁻¹
├── thermoanalysis.yaml            # ZPE, thermal corrections, S, H, G
├── vib/
│   └── mode_NNN.{xyz,pdb}         # per-mode displacement (visualize in PyMOL)
└── freq.log
```

`result.json` keys:

```python
import json
d = json.load(open("result_freq/result.json"))
print(d["n_imaginary"])                 # 0 for minimum, 1 for TS
print(d["frequencies_cm"][:5])          # first five frequencies
print(d["thermochemistry"]["zpe_ha"])
print(d["thermochemistry"]["thermal_correction_ha"])
print(d["thermochemistry"]["entropy_kcal_per_K"])
print(d["thermochemistry"]["G_hartree"])
```

## QRRHO thermochemistry

Default thermochemistry uses the QRRHO (Grimme) treatment with a
100 cm⁻¹ rotor cutoff:

- low-frequency vibrations (< 100 cm⁻¹) are interpolated toward the
  free-rotor entropy formula,
- high-frequency vibrations use the standard harmonic-oscillator
  partition function.

Inspect the active QRRHO knob via `mlmm.defaults.THERMO_KW`.

## Partial-Hessian Vibrational Analysis (PHVA)

When the input has frozen atoms (PDB B-factor or `freeze_atoms`
set), `freq` automatically computes the **partial Hessian**: only the
mobile-atom block is built and diagonalized; frozen atoms are projected
out. This is much cheaper for large clusters.

Frozen atoms are written by `extract` for link-H parents. To override,
use `--config` YAML and set `freeze_atoms`.

## Caveats

- A minimum should have **0 imaginary frequencies**, a TS should have
  **exactly 1**.
- Imaginary frequencies < ~50 cm⁻¹ are often numerical noise, not real
  modes. The QRRHO cutoff (100 cm⁻¹) is one safeguard.
- `--hessian-calc-mode FiniteDifference` is more memory-friendly but
  ~3× slower than `Analytical`.
- Thermochemistry depends on charge / spin — make sure `-q`/`-m` are
  correct or ZPE will be off.

## See also

- `tsopt.md`, `irc.md` — usual upstream stages.
- `mlmm-install-backends/uma.md` — `--hessian-calc-mode` knob.
- Defaults: `import mlmm.defaults as d; print(d.FREQ_KW, d.THERMO_KW, d.FREQ_CALC_KW)`
