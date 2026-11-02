# Contributing

SOuPE is a small Stevens operator utility for QuTiP.
Follow PEP8 style guidelines. 

```bash
uv sync
uv run ruff check src examples
uv run ruff format --check src examples
uv run pycodestyle src examples
uv build
```

Polynomial definitions follow Ballou's operator table and the publications
listed in [TRACEABILITY.md](docs/TRACEABILITY.md). Keep expressions explicit
for LaTeX export. NumPy and QuTiP are the runtime dependencies.


## Development and review

Generative coding agents were used during development to assist with
implementation, refactoring, documentation and software tests. The
LaTeX exporter in particular, was largly AI generated. Scientific
definitions and conventions are grounded in the references listed in
[TRACEABILITY.md](docs/TRACEABILITY.md). Agent-generated material is subject
to author review.

Contributions to SOuPE's original code, documentation and examples are made
under the [BSD 3-Clause License](LICENSE).
