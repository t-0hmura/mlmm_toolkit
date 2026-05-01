# `mlmm tsopt`

## Purpose

Transition-state optimization. Two algorithms: Hessian-Guided Dimer
(default, `--opt-mode grad`/`dimer`) and RS-I-RFO (full-Hessian,
`--opt-mode hess`/`rsirfo`). Use after `path-search` or `scan` to refine
a HEI to a true first-order saddle, or as a standalone validator on an
externally-generated TS guess.

## Synopsis

```bash
mlmm tsopt -i ts_guess.{pdb,xyz,gjf} \
    [-q 0 -m 1] [-l 'RES:Q,...'] \
    [--opt-mode grad|hess|light|heavy|dimer|rsirfo] \
    [--max-cycles 10000] \
    [-b uma|orb|mace|aimnet2] [-o ./result_tsopt/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

Beyond the cluster-style flags below (inherited from `pdb2reaction`),
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); overrides `--model-pdb` |
| `--ref-pdb FILE` | Full-enzyme PDB used as topology reference for XYZ inputs |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | TS candidate; `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `hess` | `grad`/`dimer` (Hessian-Guided Dimer) or `hess`/`rsirfo` (RS-I-RFO); the mlmm-only `light` / `heavy` shortcuts are also accepted (light = Dimer, heavy = full-Hessian RS-I-RFO) |
| `--max-cycles` | int | 10000 | Optimization step cap |
| `--hessian-init` | str | (live default) | `'analytical'` / `'finite-diff'` / `'guess'`; check `RSIRFO_KW` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_tsopt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default RS-I-RFO

```bash
mlmm tsopt -i hei.xyz -q 0 -m 1 -b uma -o result_tsopt
```

### Dimer mode (lighter, no full Hessian)

```bash
mlmm tsopt -i hei.xyz -q 0 -m 1 \
    --opt-mode dimer -b uma -o result_tsopt_dimer
```

### Tighter convergence on an ill-conditioned saddle

```bash
mlmm tsopt -i hei.xyz -l 'SAM:1,GPP:-3' \
    --opt-mode rsirfo --max-cycles 200 -b mace \
    -o result_tsopt_rsirfo
```

## Output

```
result_tsopt/
├── result.json
├── final_geometry.{xyz,pdb}        # converged TS
├── tsopt_trj.xyz                   # full optimization trajectory
├── vib/                            # imaginary-mode vibrations
│   └── imag_*.{pdb,xyz}            # mode displacement visualization
└── tsopt.log
```

`result.json` keys:

```python
import json
d = json.load(open("result_tsopt/result.json"))
print(d["status"])                      # "converged" / "not_converged"
print(d["energy_hartree"])
print(d["n_imaginary_modes"])           # should be 1 for a real TS
print(d["imaginary_frequencies_cm"])    # list of cm⁻¹
print(d["structure_path"])              # final_geometry.xyz
```

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `hess` / `rsirfo` (default) | RS-I-RFO with full Hessian | Robust for tricky / multi-imaginary-mode candidates; slower per cycle but converges in fewer cycles |
| `grad` / `dimer` | Hessian-Guided Dimer | Cheaper per cycle; useful when full Hessian is too expensive (large clusters > 600 atoms with UMA-m) |

Switch from `dimer` → `rsirfo` if Dimer fails to converge after ~50
cycles.

## Validation: imaginary modes

A real TS has exactly one imaginary frequency that corresponds to the
reaction coordinate.

```python
import json
d = json.load(open("result_tsopt/result.json"))
if d["n_imaginary_modes"] == 1:
    print("OK: single imaginary mode at", d["imaginary_frequencies_cm"][0], "cm-1")
elif d["n_imaginary_modes"] == 0:
    print("BAD: collapsed to a minimum during refinement")
elif d["n_imaginary_modes"] > 1:
    print("AMBIGUOUS: multiple imaginary modes; inspect vib/imag_*.pdb")
```

For multi-imaginary cases, visualize the modes (`pymol vib/imag_*.pdb`)
to decide whether the extra modes are spurious (translation/rotation
of frozen residues) or real chemical second-order saddle points.

## Caveats

- A converged `tsopt` is **not** a complete validation; always follow
  with `irc.md` to confirm the TS connects the expected R and P.
- `--max-cycles` rarely needs to be increased above 200 in practice;
  failure usually means the TS guess is too far off.
- Backend choice matters here more than for minima: UMA / MACE are
  usually safer than Orb for TS curvature.

## See also

- `path-search.md` — produces TS candidates for `tsopt`.
- `irc.md`, `freq.md` — downstream validation.
- `mlmm-install-backends/uma.md` / `mace.md` — TS-accurate
  backends.
- Defaults: `import mlmm.defaults as d; print(d.RSIRFO_KW, d.DIMER_KW, d.HESSIAN_DIMER_KW)`
