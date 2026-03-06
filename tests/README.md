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
| `test_cli_smoke_low_coverage.py` | Utility subcommands (add-elem-info, fix-altloc, trj2fig, energy-diagram) |

CI runs these via `.github/workflows/pytest.yml` on every push/PR.

## Smoke tests (GPU)

GPU smoke tests live in `tests/smoke/`. They require a CUDA-capable GPU
and `conda activate mlmm`.

```bash
cd tests/smoke
bash run.sh        # local
qsub run.sh        # PBS
```

The script covers all major subcommands (extract, define-layer, path-opt,
path-search, opt, tsopt, freq, irc, dft, scan, scan2d, scan3d, all)
plus --dry-run validation, --microiter, trj2fig, add-elem-info, and
energy-diagram.
