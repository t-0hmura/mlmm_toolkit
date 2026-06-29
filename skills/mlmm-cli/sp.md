# `mlmm sp`

## Purpose

Single-point ML/MM ONIOM energy + forces (and optionally the full
Hessian). The cheapest stage: useful for spot-checking a geometry,
extracting forces, or producing a Hessian without running an
optimization.

## Synopsis

```bash
mlmm sp -i structure.pdb --parm real.parm7 [-q 0 -m 1] \
    [--hess] [--hessian-calc-mode analytical|finitedifference] \
    [-b uma|orb|mace|aimnet2] [-o ./result_sp/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** (`--real-parm7` alias) |
| `--model-pdb FILE` | PDB defining the ML-region atoms (optional with `--detect-layer`) |
| `--detect-layer / --no-detect-layer` | Pick layer assignment from PDB B-factor (0.0=ML, 10.0=movable-MM, 20.0=frozen). Default on. |
| `--model-indices` | Comma-separated atom indices for ML region (e.g. `'1-50,75,100-110'`); used only when `--model-pdb` is omitted (`--model-pdb` takes precedence) |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `--embedcharge / --no-embedcharge` | xTB point-charge embedding for MM→ML environment (default off) |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Layered `.pdb` (or `.xyz`) |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--hess` / `--no-hess` | flag | `no-hess` | Also compute the full ONIOM Hessian and save to `hessian.npy` |
| `--hessian-calc-mode` | str | (auto) | `analytical` (UMA only) or `finitedifference`; used only with `--hess` |
| `--mm-backend` | str | `hessian_ff` | MM backend: `hessian_ff` or `openmm` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_sp/` | Output directory |
| `--out-json / --no-out-json` | flag | `no-out-json` | Write machine-readable `result.json` to out-dir |
| `--config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Energy + forces

```bash
mlmm sp -i my.pdb -l 'SAM:1' -b uma -o result_sp
```

### Energy + full Hessian

```bash
mlmm sp -i my.pdb -q -1 -m 1 --hess -o result_sp_hess
```

## Output

```
result_sp/
├── result.json      # when --out-json
├── summary.json     # mirrored payload, written alongside result.json
├── forces.npy       # ML/MM forces (a.u./bohr), shape (N, 3)
└── hessian.npy      # full ONIOM Hessian (a.u.) — only when --hess
```

`result.json` reports `stage`, `status`, `backend`, `charge`, `spin`,
`energy_au`, `forces_path`, and (when `--hess`) `hessian_path`.

## Caveats

- Single-point only — no geometry change. For relaxation use `opt.md`;
  for TS search use `tsopt.md`.
- `--hessian-calc-mode analytical` requires the UMA backend; other
  backends fall back to FiniteDifference automatically.
- `--config` YAML overrides less-common settings; inspect
  `MLMM_CALC_KW` and `GEOM_KW_DEFAULT` in `mlmm.core.defaults`.

## See also

- `opt.md` — relax the geometry to a minimum.
- `freq.md` — vibrational analysis from the Hessian.
- `dft.md` — replace the ML high level with a DFT single point.
