# Version Management

## Single Source of Truth

The project version is defined in **one single place**: `eis_analysis/version.py`

All other files import the version from there, ensuring consistency across the project.

## Structure

### `eis_analysis/version.py`

```python
__version__ = '0.9.0'              # Main version (semantic versioning)
__version_info__ = (0, 9, 0)       # Tuple for programmatic comparison
__release_date__ = '2025-12-30'    # Release date
__breaking_changes__ = [...]       # List of breaking changes (for major versions)
```

### Functions

- `get_version_string()` -> `"v0.9.0 (2025-12-30)"` - formatted string for UI

## Usage

### Python modules

```python
# In modules
from eis_analysis import __version__, get_version_string

print(f"Version: {__version__}")  # 0.9.0
print(get_version_string())       # v0.9.0 (2025-12-30)
```

### CLI

CLI `eis.py` automatically imports and displays the version:

```bash
python eis.py data.DTA
# INFO: EIS Analysis (v0.9.0 (2025-12-30))
```

### Documentation

Documentation files should reference the version as:

```markdown
**Version:** v0.9.0 (2025-12-30)
```

Or dynamically import from `version.py` if generated.

## Files with Versions

### Automatically synchronized (import)

+ `eis_analysis/__init__.py` - imports `__version__`
+ `eis.py` - imports and displays `get_version_string()`

### Manually maintained

- `README.md` - line 3: `**Version:** vX.Y.Z`
- `doc/PYTHON_API.md` - line 3: `**Current version:** vX.Y.Z`
- `CHANGELOG.md` - newest entry: `## Version X.Y.Z`

### Legacy/ignore

- `doc/*.md` - feature-specific docs (have their own versions)

## Workflow for Version Change

### 1. Update `eis_analysis/version.py`

```python
__version__ = '1.0.0'              # New version
__version_info__ = (1, 0, 0)
__release_date__ = '2026-01-15'    # New date
```

### 2. Update `CHANGELOG.md`

Add new entry at the beginning:

```markdown
## Version 1.0.0 (2026-01-15)

### Added
...
```

### 3. Update README.md

Line 3:
```markdown
**Version:** v1.0.0 (2026-01-15)
```

### 4. Update PYTHON_API.md

Line 3:
```markdown
**Current version:** v1.0.0
```

### 5. Update `eis_analysis/__init__.py` docstring

Line 2:
```python
"""
EIS Analysis Toolkit v1.0.0
============================
...
```

### 6. Verify

```bash
# Check that import works
python3 -c "from eis_analysis import __version__, get_version_string; print(__version__); print(get_version_string())"

# Check CLI
python3 eis.py --no-fit -v | head -5
```

## Semantic Versioning

The project uses [Semantic Versioning 2.0.0](https://semver.org/):

- **MAJOR** (X.0.0) - Breaking changes (incompatible API changes)
- **MINOR** (x.Y.0) - New features (backwards compatible)
- **PATCH** (x.y.Z) - Bug fixes (backwards compatible)

### Examples

- `0.9.0` -> `0.9.1` - Bug fix (patch)
- `0.9.1` -> `0.10.0` - New feature (minor)
- `0.10.0` -> `1.0.0` - Breaking change or stable release (major)

**Breaking change example:**
- Removal of a public function
- Change of return format

## Consistency Check

To verify that all files have the correct version:

```bash
# Find all versions in project
grep -rn "v[0-9]\+\.[0-9]\+\.[0-9]\+" --include="*.py" --include="*.md" | grep -v CHANGELOG

# Should return consistent version (current: v0.9.0)
```

## Notes

- **Never** change version in multiple files at once - always start from `version.py`
- **CHANGELOG.md** must have an entry for each version
- **Major version bump** requires breaking change documentation
- CLI automatically displays version at startup (INFO level)
