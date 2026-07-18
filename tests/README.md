# Tests Layout

`mlmm-toolkit` uses `tests/` as the primary test root (the bundled `hessian_ff/tests/` is a second root configured in `pyproject.toml` `testpaths`).

## Unit and CI Tests

CPU-only tests live directly under `tests/*.py`.

```bash
pytest tests/ -v --tb=short -x
```

These tests cover CLI parsing, helper contracts, logging summaries, geometry
regressions, frequency active-DOF behavior. They are suitable for CI.

## Manual Smoke Tests

GPU smoke fixtures and commands live in `tests/smoke/`.

```bash
cp -a tests/smoke /path/to/writable-scratch/mlmm-smoke
cd /path/to/writable-scratch/mlmm-smoke
bash run.sh
```

`tests/smoke/run.sh` assumes the caller has already configured the Python
environment, CUDA runtime, and AmberTools. It verifies that distribution
metadata, the imported module, and the module CLI report one consistent
version, without pinning a release number. Every CLI case is executed through
`python -m mlmm` from the same interpreter. Scheduler wrappers and environment
activation stay out of the repository.

The default strict lane numerically compares UMA and ORB analytical Hessians
with finite differences. Run the same required checker in the isolated MACE
and AIMNet2 environments:

```bash
bash run_backend_hessian.sh mace
bash run_backend_hessian.sh aimnet2
```
