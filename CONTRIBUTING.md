# Contributing to mlmm_toolkit

Thank you for your interest in contributing to **mlmm_toolkit**.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/t-0hmura/mlmm_toolkit.git
cd mlmm_toolkit

# Create a conda environment
conda create -n mlmm python=3.11
conda activate mlmm

# Install in editable mode
pip install -e .
```

## Running Tests

```bash
# Run the full test suite
pytest tests/ -v --tb=short

# Run with coverage
pytest tests/ -v --tb=short --cov=mlmm_toolkit --cov-report=term-missing
```

## Code Style

- Use `from __future__ import annotations` in all modules.
- Use `logging.getLogger(__name__)` for structured logging; avoid bare `print()`.
- CLI options should include `show_default=True` where applicable.
- Warning messages via `click.echo()` should use `err=True` to write to stderr.

## Submitting Changes

1. Fork the repository and create a feature branch.
2. Make your changes and add tests where appropriate.
3. Ensure all tests pass locally.
4. Open a pull request with a clear description of the changes.

## License

By contributing, you agree that your contributions will be licensed under the
[GPL-3.0 License](LICENSE).
