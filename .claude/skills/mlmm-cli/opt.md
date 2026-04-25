# `mlmm opt`

## Purpose

Single-structure geometry optimization with LBFGS or RFO.
Use this to relax a starting geometry to its nearest local minimum
before feeding it to `path-search` / `path-opt`, or as a post-IRC
endpoint refinement.

## Synopsis

```bash
mlmm opt -i input.pdb [-q 0 -m 1] \
    [--opt-mode grad|hess|light|heavy|lbfgs|rfo] \
    [-b uma|orb|mace|aimnet2] [-o ./result_opt/]
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
| `-i, --input` | path | required | `.pdb` / `.xyz` / `.gjf` |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--opt-mode` | str | `grad` | `grad` (LBFGS) or `hess` (RFO); aliases `lbfgs` / `rfo` and the mlmm-only `light` / `heavy` shortcuts (light = LBFGS, heavy = RFO with full Hessian) are accepted |
| `--max-cycles` | int | (live default) | Stop after N cycles; check `OPT_BASE_KW` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `--solvent` | str | none | xTB-ALPB solvent |
| `-o, --out-dir` | path | `./result_opt/` | Output directory |
| `--config` / `--show-config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Default LBFGS

```bash
mlmm opt -i my.pdb -l 'SAM:1' -b uma -o result_opt
```

### RFO for stiffer convergence

```bash
mlmm opt -i my.xyz -q -1 -m 1 --opt-mode rfo -b mace -o result_opt_rfo
```

### Pre-relax endpoints before path-opt

```bash
mlmm opt -i 1.R.pdb -l '...' -o /tmp/relax_R
mlmm opt -i 3.P.pdb -l '...' -o /tmp/relax_P
mlmm path-opt -i /tmp/relax_R/final_geometry.xyz /tmp/relax_P/final_geometry.xyz ...
```

## Output

```
result_opt/
├── result.json
├── final_geometry.{xyz,pdb}   # converged geometry
├── opt_trj.xyz                # full optimization trajectory
└── opt.log                    # per-cycle gradient norm + energy
```

`result.json` keys: `status` (converged / not_converged), `n_cycles`,
`final_energy_hartree`, `gradient_max`, `structure_path`.

## `--opt-mode` choice

| Mode | Algorithm | When |
|---|---|---|
| `grad` / `lbfgs` | LBFGS | Default, fast, robust for most well-conditioned minima |
| `hess` / `rfo` | RFO with Hessian updates | Stiffer convergence; useful when LBFGS oscillates |

## Caveats

- Not a TS optimizer — for TS use `tsopt.md`.
- LBFGS occasionally walks past a saddle on shallow surfaces; if the
  resulting geometry has imaginary frequencies (run `freq` to check),
  re-run with `--opt-mode rfo`.
- `--config` YAML is the way to override less-common settings (step
  limits, trust radius, etc.); inspect `OPT_BASE_KW` and `LBFGS_KW`
  in `mlmm.defaults`.

## See also

- `tsopt.md` — TS analog.
- `freq.md` — verify the optimized minimum (zero imaginary modes).
- Defaults: `import mlmm.defaults as d; print(d.OPT_BASE_KW, d.LBFGS_KW, d.RFO_KW)`