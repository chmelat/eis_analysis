# Version Management

## Single Source of Truth

Project version is defined in `eis_analysis/version.py`:

```python
__version__ = '0.9.2'
__version_info__ = (0, 9, 2)
__release_date__ = '2026-01-05'
__breaking_changes__ = []  # For major versions
```

All other files import from here or reference it.

## Files with Versions

**Automatically synchronized (import/dynamic):**
- `eis_analysis/__init__.py` - imports `__version__`
- `eis.py` - imports `get_version_string()`
- `pyproject.toml` - dynamic: `version = {attr = "eis_analysis.version.__version__"}`

**Manually maintained:**
- `README.md` line 3: `**Version:** vX.Y.Z (YYYY-MM-DD)`
- `doc/PYTHON_API.md` line 3: `**Current version:** vX.Y.Z`
- `CHANGELOG.md` - newest entry: `## Version X.Y.Z (YYYY-MM-DD)`

## Version Change Workflow

1. **Update `eis_analysis/version.py`** - change `__version__`, `__version_info__`, `__release_date__`
2. **Update `CHANGELOG.md`** - add new entry at top with changes
3. **Update `README.md` line 3** - version + date
4. **Update `doc/PYTHON_API.md` line 3** - version only
5. **Verify**: `python3 -c "from eis_analysis import __version__; print(__version__)"`

## Semantic Versioning

Following [semver.org](https://semver.org/):
- **MAJOR** (X.0.0) - breaking changes
- **MINOR** (x.Y.0) - new features (backwards compatible)
- **PATCH** (x.y.Z) - bug fixes

## Consistency Check

```bash
grep -rn "v[0-9]\+\.[0-9]\+\.[0-9]\+" --include="*.py" --include="*.md" | grep -v CHANGELOG
```

## Important Notes

- **Always** start from `version.py` - never edit multiple files at once
- `pyproject.toml` imports automatically - do not edit manually
- `__init__.py` imports version - no hardcoded version in docstring
- Major version bump requires `__breaking_changes__` documentation
- CHANGELOG.md must have entry for every version
