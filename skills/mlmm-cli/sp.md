# `mlmm sp`

## Purpose

Single-point ML/MM ONIOM energy + forces (and optionally the full
Hessian). The cheapest stage: useful for spot-checking a geometry,
extracting forces, or producing a Hessian without running an
optimization.

## Synopsis

```bash
mlmm sp -i structure.pdb --parm real.parm7 [-q 0 -m 1] \
    [--hess] [--hessian-calc-mode Analytical|FiniteDifference] \
    [-b uma|orb|mace|aimnet2] [-o ./result_sp/]
```


## ML/MM-aware flags (mlmm-toolkit specific)

In addition to the common flags below,
**`mlmm-toolkit` requires an Amber topology** and supports layer-aware
selection. Most subcommands accept:

| flag | purpose |
|---|---|
| `--parm FILE` | Amber `parm7` topology of the whole enzyme — **required** (`--real-parm7` alias) |
| `--model-pdb FILE` | Explicit ML-region PDB; takes precedence over B-factor ML membership |
| `--detect-layer` | Automatically read B-factor layers; explicit ML membership retains valid movable/frozen MM layers. Enabled by default. |
| `--model-indices` | Explicit ML atom indices used when `--model-pdb` is omitted; takes precedence over B-factor ML membership |
| `--link-atom-method [scaled\|fixed]` | g-factor (default) or fixed 1.09/1.01 Å |
| `-q, --charge` | **ML-region** charge (not whole-system) |
| `-l, --ligand-charge` | Per-residue charge mapping for ML region |

Inspect via `mlmm <subcommand> --help` and `mlmm <subcommand> --help-advanced`.

## Key flags

| flag | type | default | description |
|---|---|---|---|
| `-i, --input` | path | required | Layered `.pdb`, or `.xyz` with `--ref-pdb` |
| `--ref-pdb` | path | required for XYZ | Atom-order-identical full-system PDB/mmCIF topology and layer metadata |
| `-q` / `-l` / `-m` | — | — | Charge / spin (common conventions) |
| `--hess` / `--no-hess` | flag | `no-hess` | Also compute the active-coordinate ONIOM Hessian block and save it to `hessian.npy` |
| `--hessian-calc-mode` | str | (auto) | `Analytical` (UMA/ORB/MACE/AIMNet2) or `FiniteDifference`; used only with `--hess` |
| `--mm-backend` | str | `hessian_ff` | MM backend: `hessian_ff` or `openmm` |
| `-b, --backend` | str | `uma` | MLIP backend |
| `-o, --out-dir` | path | `./result_sp/` | Output directory |
| `--out-json / --no-out-json` | flag | `no-out-json` | Write machine-readable `result.json` to out-dir |
| `--config` / `--dry-run` / `--help-advanced` | — | — | Standard |

## Examples

### Energy + forces

For XYZ coordinates, supply the matching full-system topology:

```bash
mlmm sp -i structure.xyz --ref-pdb structure.pdb --parm real.parm7 \
  -q 0 -m 1 -o result_sp
```

```bash
mlmm sp -i my.pdb --parm real.parm7 -l 'SAM:1' -b uma -o result_sp
```

### Energy + active-coordinate Hessian

```bash
mlmm sp -i my.pdb --parm real.parm7 -q -1 -m 1 --hess -o result_sp_hess
```

## Output

```
result_sp/
├── result.json      # when --out-json
├── summary.json     # mirrored payload, written alongside result.json
├── forces.npy       # ML/MM forces (a.u./bohr), shape (N, 3)
└── hessian.npy      # active-coordinate ONIOM Hessian block (a.u.) — only with --hess
```

`result.json` reports `stage`, `status`, `mlip_backend`, `mlip_model`,
`mlip_precision`, `mm_backend`, `link_atom_method`, `use_cmap`, `charge`,
`spin`, `energy_au`, `forces_path`, and `hessian_path` (null without
`--hess`). For `--calc-file`, `mlip_backend` is `custom`, `mlip_model` is
`filename:factory`, and `mlip_precision` is null.

## Caveats

- Single-point only — no geometry change. For relaxation use `opt.md`;
  for TS search use `tsopt.md`.
- An explicit `--hessian-calc-mode Analytical` requires the selected
  backend's analytical API. If it is unavailable, the command raises an error;
  it does not silently switch to `FiniteDifference`.
- `--config` YAML overrides less-common settings; inspect
  `MLMM_CALC_KW` and `GEOM_KW_DEFAULT` in `mlmm.core.defaults`.

## See also

- `opt.md` — relax the geometry to a minimum.
- `freq.md` — vibrational analysis from the Hessian.
- `dft.md` — replace the ML high level with a DFT single point.
