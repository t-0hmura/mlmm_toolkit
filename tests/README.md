# Tests Layout

`mlmm_toolkit` uses `tests/` as the single test root.

- Unit tests: `tests/*.py`
- Smoke fixtures/scripts: `tests/smoke/`

## Smoke tests

GPU smoke tests live in `tests/smoke/`. They require a CUDA-capable GPU
and `conda activate mlmm`.

```bash
cd tests/smoke
bash run.sh        # local
qsub run.sh        # PBS
```

The script runs 21+ numbered test cases covering all major subcommands
(extract, define-layer, path-opt, path-search, opt, tsopt, freq, irc,
dft, scan, scan2d, scan3d, all) plus --dry-run validation, --config,
--microiter, trj2fig, add-elem-info, and energy-diagram.
