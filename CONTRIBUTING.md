# Contributing to SeqForge

Thanks for your interest in improving SeqForge! This project follows a
strict **zero runtime dependency** policy: everything must run on the
Python 3.8+ standard library alone.

## Ground rules

1. **No third-party dependencies.** Standard library only. If you need a
   feature that seems to require a dependency, open an issue first.
2. **Python 3.8 compatibility.** No `match` statements, no PEP 604 unions
   in runtime annotations, no `str.removeprefix`.
3. **Deterministic behaviour.** No network access, no randomness without a
   seed, no wall-clock dependence in library code.
4. **Tests are mandatory.** Every bug fix and every new feature must come
   with unit tests in `tests/`.

## Development workflow

```bash
# 1. clone your fork
git clone https://github.com/gitstq/seqforge.git
cd seqforge

# 2. create a branch (Angular commit convention)
git checkout -b feat/my-feature

# 3. run the full test suite
python -m unittest discover -s tests -v

# 4. syntax check every file (same gate as CI)
python -m compileall -q seqforge tests
```

## Commit message convention

We use the Angular convention:

- `feat:` a new feature
- `fix:` a bug fix
- `docs:` documentation only
- `test:` tests only
- `refactor:` code change that neither fixes a bug nor adds a feature
- `chore:` tooling / CI / packaging

Example: `feat: add --min-len filter to fqtrim`

## Pull request checklist

- [ ] `python -m unittest discover -s tests` passes locally
- [ ] New behaviour is covered by tests
- [ ] Public functions have docstrings
- [ ] README updated when the CLI surface changes
- [ ] No new dependencies in `pyproject.toml`

## Reporting issues

Please include:

- SeqForge version (`seqforge --version`) and Python version
- The exact command you ran and a **minimal** input file that reproduces
  the problem (a few bases is usually enough)
- Expected vs actual output

Security-sensitive reports (e.g. a parser that crashes on malformed
input) can be sent privately via GitHub Security Advisories.
