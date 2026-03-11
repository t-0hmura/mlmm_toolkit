# Tests Layout

`mlmm_toolkit` uses `tests/` as the single test root.

- Unit tests: `tests/*.py`
- Smoke fixtures/scripts: `tests/smoke/`

## Unit tests (CPU-only)

Run with pytest on any machine (no GPU required):

```bash
pytest tests/ -v --tb=short -x
```

Key test files:

| File | Coverage |
|------|----------|
| `test_cli_help.py` | `--help` / `--help-advanced` for all subcommands |
| `test_default_group.py` | DefaultGroup behavior (default subcmd, lazy import, bool compat) |
| `test_bool_compat_cli.py` | Bool toggle/value syntax across all subcommands |

CI runs these via `.github/workflows/pytest.yml` on every push/PR.

## Smoke tests (GPU)

GPU smoke tests live in `tests/smoke/`. They require a CUDA-capable GPU
and `conda activate mlmm`.

```bash
cd tests/smoke
bash run.sh        # local
qsub run.sh        # PBS
```

### Speed optimizations

All iterative tests use `--thresh gau_loose` and reduced `--max-cycles` (3-5)
to keep runtime short. Scan tests use `--preopt False --endopt False`, and
path commands use `--preopt False --climb False --max-nodes 5`.

### Test coverage (34 tests)

| # | Category | Command | Notes |
|---|----------|---------|-------|
| 1 | Preparation | `extract` | pocket extraction |
| 2 | Preparation | `define-layer` | B-factor layer assignment |
| 3 | Preparation | `mm-parm` | AMBER parm7 generation |
| 4 | Subcommand | `opt` (grad/lbfgs) | `--dump` for trj2fig |
| 5 | Subcommand | `opt` (hess/rfo) | |
| 6 | Subcommand | `opt` (hess, `--microiter`) | microiteration |
| 7 | Subcommand | `tsopt` (grad/dimer) | |
| 8 | Subcommand | `tsopt` (hess/rsirfo) | |
| 9 | Subcommand | `freq` | |
| 10 | Subcommand | `irc` | max-cycles 3 |
| 11 | Subcommand | `dft` (hf/sto-3g) | |
| 12 | Subcommand | `scan` (1D) | preopt/endopt False |
| 13 | Subcommand | `scan2d` | |
| 14 | Subcommand | `scan3d` | |
| 15 | Subcommand | `path-opt` (gsm) | max-nodes 5 |
| 16 | Subcommand | `path-opt` (dmf) | |
| 17 | Subcommand | `path-search` | |
| 18 | all command | `all` (minimal, no tsopt/thermo/dft) | |
| 19 | all command | `all` (full: tsopt+thermo+dft) | |
| 20 | all command | `all` (`--ref-pdb --parm --model-pdb`) | XYZ input, tsopt-only from test19 artifacts |
| 21 | tsopt variant | `tsopt` (`--radius-hessian 0.0`) | |
| 22 | tsopt variant | `tsopt` (`--radius-hessian 3.6`) | |
| 23-29 | Dry-run | opt/tsopt/freq/scan/dft/path-search/irc | `--dry-run` |
| 30 | Utility | `add-elem-info` | |
| 31 | Utility | `trj2fig` | conditional on test4 output |
| 32 | Utility | `energy-diagram` | |
| 33 | Utility | `oniom-export` (g16) | |
| 34 | xTB | `opt` (`--embedcharge`) | requires xTB binary |

## `oniom-export` Tips (Gaussian ONIOM)

### Default route line

`oniom-export` generates the following route section by default:

```
#p oniom(wB97XD/def2-TZVPD:amber=softonly)
scf=(xqc,intrep,maxconventionalcyc=80)
nosymm iop(2/15=3) geom=connectivity Amber=(FirstEquiv)
```

Users should add `opt=` or `opt(ts,...)` keywords manually before submission.

### Keyword reference

| Keyword | Meaning |
|---------|---------|
| `scf=(xqc,intrep,maxconventionalcyc=80)` | `xqc`: fall back to quadratically-convergent SCF if DIIS fails. `intrep`: use internal representation for SCF. `maxconventionalcyc=80`: allow up to 80 conventional SCF cycles before switching to QC-SCF. |
| `iop(2/15=3)` | Read MM atom types from input as-is (do not let Gaussian reassign them). Required for AMBER parm7-derived atom types (`CT`, `HC`, `N`, `OS`, etc.). |
| `nosymm` | Disable symmetry (required for ONIOM with frozen atoms). |
| `geom=connectivity` | Read connectivity table from input (appended after coordinates). |
| `Amber=(FirstEquiv)` | Use AMBER force field with FirstEquiv atom type matching. |

### Optimization keywords

For geometry optimization, add one of these to the route line:

```
opt=(quadmac,calcfc,maxstep=10)
```

For TS optimization:

```
opt=(calcfc,quadmac,ts,noeigentest,recalcfc=200)
```

A smaller basis set (e.g., `wB97XD/def2-SVP`) is recommended for optimization,
with single-point energy at `wB97XD/def2-TZVPD` afterwards.

### EmbedCharge

When the mlmm calculation used `--embedcharge` (QM/MM electrostatic embedding),
add `=EmbedCharge` to the ONIOM keyword:

```
#p oniom(wB97XD/def2-TZVPD:amber=softonly)=EmbedCharge
```
