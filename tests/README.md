# Tests Layout

`mlmm-toolkit` uses `tests/` as the single test root.

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
cd tests/smoke
bash run.sh
```

`tests/smoke/run.sh` assumes the caller has already configured the Python
environment, CUDA runtime, and AmberTools when needed. Scheduler wrappers and
environment activation must stay out of the repository.

