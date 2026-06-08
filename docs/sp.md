# `sp`

`mlmm sp` evaluates the ML/MM ONIOM energy + atomic forces (optionally the full ONIOM Hessian) at a single geometry.

## When to use

- Fast inspection of a layered structure before running an optimization.
- Comparing backends head-to-head on the same ONIOM partition.
- Generating reference Hessians outside the optimizer loop.

## Quick examples

```bash
# energy + forces on a layered PDB (B-factor encodes ML / movable-MM / frozen-MM)
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1
```

```bash
# also compute the full ONIOM Hessian (Analytical when --backend uma)
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1 --hess
```

## Inputs

Command form:

```bash
mlmm sp -i INPUT --parm PARM7 -q CHARGE [options]
```

| Input | Required | Notes |
|---|---|---|
| `-i, --input FILE` | yes | layered PDB (or XYZ) defining the ML / movable-MM / frozen-MM partition |
| `--parm FILE` | yes | Amber `parm7` topology of the full enzyme (`--real-parm7` retained as alias) |
| `-q, --charge INT` | yes | ML region total charge |
| `-l, --ligand-charge TEXT` | no | per-ligand charge mapping (e.g. `SAM:1,GPP:-3`) |
| `-m, --multiplicity INT` | no | ML region spin multiplicity, 2S+1 (default `1`) |

### ML region selection

Either embed the partition in the input PDB's B-factor (ML=0.0, movable-MM=10.0, frozen-MM=20.0) with `--detect-layer` (the default), or pass it explicitly:

| flag | meaning |
|---|---|
| `--detect-layer / --no-detect-layer` | use B-factor encoding (default `on`) |
| `--model-pdb FILE` | alternative PDB defining ML atoms |
| `--model-indices TEXT` | comma-separated 1-based atom indices (e.g. `1-50,75,100-110`) |

### Hessian backend

When `--hess` is set, the backend choice picks the Hessian computation strategy:

- `--backend uma` (default) → `Analytical` Hessian for the ML region via the UMA torch autograd path; MM region uses `hessian_ff` analytical Hessian
- `--backend orb` / `mace` / `aimnet2` → falls back to `FiniteDifference` for the ML region

`--hessian-calc-mode` lets you override per-call.

## Outputs

`sp` writes outputs under `result_sp/` by default. The ONIOM energy is also printed to stdout; the JSON files (identical payload mirrored to both names) are emitted only when `--out-json` is passed.

| file | contents | written |
|---|---|---|
| `forces.npy` | `(N, 3)` array of ONIOM forces in atomic units (Hartree / Bohr) | always |
| `hessian.npy` | `(3N, 3N)` mass-unweighted ONIOM Hessian (Hartree / Bohr²) | only with `--hess` |
| `result.json` / `summary.json` | ONIOM energy (a.u.), backend, charge/spin, paths to npy outputs, elapsed time | only with `--out-json` |

`sp` does not write a `summary.log`.

## CLI options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| flag | default | meaning |
|---|---|---|
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | MLIP backend for the ML region |
| `--hess / --no-hess` | `--no-hess` | also compute and write `hessian.npy` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | auto | force a specific Hessian mode (only with `--hess`) |
| `--embedcharge / --no-embedcharge` | off | xTB point-charge embedding correction for MM→ML coupling |
| `--link-atom-method [scaled\|fixed]` | `scaled` | link-atom positioning |
| `--mm-backend [hessian_ff\|openmm]` | `hessian_ff` | MM backend (analytical vs finite-difference Hessian) |
| `-o, --out-dir PATH` | `./result_sp/` | output directory |
| `--precision [fp32\|fp64]` | `fp32` | numeric precision passed to the backend |
| `--config PATH` | — | YAML config providing `calc.*`, `geom.*` defaults |
| `--show-config / --dry-run` | off | print effective merged config / validate without running |

Run `mlmm sp --help-advanced` for the full list (hess-cutoff override, MCP-style result.json, etc.).

## See Also

- [`opt`](opt.md) — optimize the layered structure (microiteration)
- [`tsopt`](tsopt.md) — refine a TS candidate (ML/MM ONIOM)
- [`freq`](freq.md) — ONIOM vibrational analysis + QRRHO thermochemistry
- [`dft`](dft.md) — single-point DFT counterpart on the ML region
