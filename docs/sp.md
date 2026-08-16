# `sp`

`mlmm sp` evaluates the ML/MM ONIOM energy + atomic forces (optionally the active-coordinate ONIOM Hessian) at a single geometry. Use it for fast inspection of a layered structure before running an optimization, for comparing backends directly on the same ONIOM partition, or for generating reference Hessians outside the optimizer loop.

## Examples

Energy + forces on a layered PDB (B-factor encodes ML / movable-MM / frozen):

```bash
# energy + forces on a layered PDB (B-factor encodes ML / movable-MM / frozen)
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1
```

Also compute the active-coordinate ONIOM Hessian:

```bash
# finite differences are used by default; select Analytical only for a backend that supports it
mlmm sp -i layered.pdb --parm real.parm7 -q 0 -m 1 --hess
```

## Outputs

`sp` writes outputs under `result_sp/` by default. The ONIOM energy is also printed to stdout; the JSON files (written to both `result.json` and `summary.json` with identical content) are emitted only when `--out-json` is passed.

| file | contents | written |
|---|---|---|
| `forces.npy` | `(N, 3)` array of ONIOM forces in atomic units (Hartree / Bohr) | always |
| `hessian.npy` | Mass-unweighted ONIOM Hessian for the calculator's active coordinates (Hartree / Bohr²); inspect the saved array shape | only with `--hess` |
| `result.json` / `summary.json` | ONIOM energy (a.u.), backend, charge/spin, paths to npy outputs, elapsed time | only with `--out-json` |

`sp` does not write a `summary.log`.

## CLI options

Command form:

```bash
mlmm sp -i INPUT --parm PARM7 -q CHARGE [options]
```

| Input | Required | Notes |
|---|---|---|
| `-i, --input FILE` | yes | layered PDB/mmCIF, or XYZ coordinates accompanied by `--ref-pdb` |
| `--ref-pdb FILE` | for XYZ | atom-order-identical full-system PDB/mmCIF supplying topology and layer metadata |
| `--parm FILE` | yes | Amber `parm7` topology of the full enzyme (`--real-parm7` retained as alias) |
| `-q, --charge INT` | yes (unless `-l` is given) | ML region total charge |
| `-l, --ligand-charge TEXT` | no | per-ligand charge mapping (e.g. `SAM:1,GPP:-3`); derives the net charge when `-q` is omitted |
| `-m, --multiplicity INT` | no | ML region spin multiplicity, 2S+1 (default `1`) |

### ML region selection

Either embed the partition in the input PDB's B-factor (ML=0.0, movable-MM=10.0, frozen=20.0) with `--detect-layer` (the default), or pass it explicitly:

| flag | meaning |
|---|---|
| `--detect-layer` | automatic B-factor layer detection (enabled by default) |
| `--model-pdb FILE` | alternative PDB defining ML atoms |
| `--model-indices TEXT` | comma-separated 1-based atom indices (e.g. `1-50,75,100-110`) |

### Hessian backend

When `--hess` is set, `--hessian-calc-mode Analytical` uses the selected
backend's analytical/native Hessian path (UMA, ORB, MACE, or AIMNet2), while
`FiniteDifference` uses central differences of forces. The MM backend defaults
to `hessian_ff`, but MM Hessians use finite differences by default. Set
`calc.mm_fd: false` for the `hessian_ff` analytical MM Hessian. An unavailable
requested path is an error.

### Other options

The full flag list is in the generated [command reference](reference/commands/index.md); the table below covers the options that need explanation.

| flag | default | meaning |
|---|---|---|
| `-b, --backend [uma\|orb\|mace\|aimnet2]` | `uma` | MLIP backend for the ML region |
| `--hess / --no-hess` | `--no-hess` | also compute and write `hessian.npy` |
| `--hessian-calc-mode [Analytical\|FiniteDifference]` | `FiniteDifference` | Hessian mode when `--hess` is set; `Analytical` uses the backend's native path |
| `--link-atom-method [scaled\|fixed]` | `scaled` | link-atom positioning |
| `--mm-backend [hessian_ff\|openmm]` | `hessian_ff` | MM backend; Hessian method is controlled separately by `calc.mm_fd` |
| `-o, --out-dir PATH` | `./result_sp/` | output directory |
| `--precision [fp32\|fp64]` | backend-specific | numeric precision passed to the backend (unset: UMA/AIMNet2 fp32, ORB/MACE fp64) |
| `--config PATH` | — | YAML config providing `calc.*`, `geom.*` defaults |
| `--show-config / --dry-run` | off | print effective merged config / validate without running |

Run `mlmm sp --help-advanced` for the full list (hess-cutoff override, MCP-style result.json, etc.).

## See Also

- [`opt`](opt.md) — optimize the layered structure (microiteration)
- [`tsopt`](tsopt.md) — refine a TS candidate (ML/MM ONIOM)
- [`freq`](freq.md) — ONIOM vibrational analysis + QRRHO thermochemistry
- [`dft`](dft.md) — single-point DFT counterpart on the ML region
