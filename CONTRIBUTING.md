# Contributing to StaDE

Thanks for your interest! This project is a course assignment that has been
released as open source, so contributions are welcome but the review bandwidth
is limited.

## Development setup

```bash
git clone https://github.com/your-account/stade.git
cd stade
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate
pip install -e ".[dev]"
python -m pytest -q
```

## What makes a good pull request

- **One change per PR**, with a clear description of the statistical or
  engineering rationale.
- **Tests**: new statistics need a test with a known-value or property-based
  check. The suite must stay green on Windows and Linux (see CI).
- **No silent behaviour changes**: if a default changes, update `CHANGELOG.md`
  and `docs/method.md`.
- **Reproducibility**: any new random procedure must accept a `seed` argument.

## Reporting bugs

Please include:

1. the command you ran (or a minimal Python snippet),
2. the full traceback,
3. the output of `python -c "import stade; print(stade.__version__)"`,
4. the shape of your input matrix and the head of your metadata table
   (no personal or patient-identifiable data).

## Code style

- Python >= 3.10, type hints encouraged, `from __future__ import annotations`.
- Pure numpy/pandas/scipy in the statistics layer; plotting stays in `plotting.py`.
- Keep functions small and documented with a one-paragraph docstring that states
  the statistical definition and the assumptions.

## License

By contributing you agree that your contribution is licensed under the MIT
License (see `LICENSE`).
