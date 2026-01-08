# CLAUDE.md

This file provides guidance for Claude Code (claude.ai/code) when working on this project.

## Project Philosophy

### Goals

**EIS Analysis Toolkit** is a modular tool for electrochemical impedance spectroscopy analysis with emphasis on:

1. **Scientific accuracy** - Mathematically correct implementation of methods (DRT, GCV, GMM)
2. **Automation** - Minimizing manual tuning (auto-lambda, auto-circuit, GMM peaks)
3. **Reproducibility** - Objective, data-driven parameter choices
4. **Usability** - Easy to use for researchers (CLI + Python API)
5. **Modular design** - Clean separation of concerns, testability

### Design Principles

**Key principles:**
- **Modularity** - Each module has clearly defined responsibility
- **Single Responsibility** - One module = one functionality
- **DRY (Don't Repeat Yourself)** - Code is not duplicated
- **Testability** - Modules testable independently
- **Documentation** - Every function has docstring, every module has documentation

**Anti-patterns to avoid:**
- Monolithic files (>1000 lines)
- God objects (classes doing too much)
- Duplicate code between modules
- Hardcoded constants without explanation
- Magic numbers without documentation

## Code Organization

**Modular architecture:**
- **Core package** - Main package with analytical functions
  - I/O modules - Data loading, synthetic data
  - Validation - Data quality (KK validation)
  - Algorithms - Core scientific computations (DRT, GCV, fitting)
  - Visualization - Plots and diagnostics
  - Domain logic - Specific analyses (e.g., oxide layers)
- **CLI script** - User interface, workflow orchestration
- **Utility scripts** - Standalone tools for data generation

**Separation rules:**
- CLI code separated from core logic (enables use as library)
- I/O separated from computations (testing without files)
- Visualization separated from algorithms (headless use)
- Domain logic isolated (extensibility for other domains)

## Documentation Organization

### Document Types

**Main documentation (README.md):**
- Structured by chapters according to functionalities:
  - Introduction and quick start
  - Kramers-Kronig validation (theory + usage)
  - DRT analysis (main chapter with theory, implementation, interpretation)
  - Circuit fitting (theory, auto-suggestion, usage)
  - Oxide analysis (domain-specific)
  - CLI reference, synthetic data, troubleshooting
  - Installation (at the end)
- Each chapter: brief theory -> implementation -> usage -> examples
- **WITHOUT** change history (-> CHANGELOG.md)
- **WITHOUT** Python API (-> PYTHON_API.md)
- **WITHOUT** emoji (due to PDF export)

**Change history (CHANGELOG.md):**
- Complete history of all versions
- Breaking changes, new features, bugfixes
- Migration guides between versions
- Chronologically from newest

**Python API documentation (PYTHON_API.md):**
- Complete Python API reference
- Usage as library (import)
- Documentation of all functions and parameters
- Examples of programmatic usage
- **WITHOUT** CLI documentation (-> README.md)

**Specialized guides (`./doc/` directory):**
- In-depth documentation of one function/module
- Mathematical background, algorithms
- Usage examples, interpretation
- Literature references

**Development guide (CLAUDE.md - this file):**
- Project philosophy
- Design principles
- Code style guidelines
- Development workflow
- **WITHOUT** implementation details (-> README.md)
- **WITHOUT** specific versions (general rules)

### Documentation Separation Rules

**Single source of truth principle:**
- Each piece of information has one primary location
- Other documents reference it, don't duplicate it
- Change in one place = change everywhere

**What belongs where:**
- Version history -> CHANGELOG.md (not README)
- Implementation details of current version -> README.md
- Mathematical theory of specific method -> specialized MD
- API reference -> docstrings in code
- Philosophy and design -> CLAUDE.md (this file)

**Consistency across documents:**
- Unified terminology
- Consistent formatting (especially for PDF export)
- Cross-references between documents (markdown links)

### Characters for PDF Export

**CRITICAL RULE: NO EMOJI in documents intended for PDF!**

DejaVu Serif font (used by pandoc for PDF) doesn't support emoji. Every emoji causes PDF warning.

**USE (PDF-safe):**
- Standard ASCII characters and Markdown: `#`, `##`, `**bold**`, `*italic*`, `` `code` ``
- Tables: `|`, `-`
- Lists: `-`, `+`, `*`, `1.`
- Code blocks: ` ```python `
- Math in code: `` `Z = R + jωL` ``
- Basic Unicode: ×, ±, ≈, ≤, ≥, →, °, ∞
- Greek letters: α, β, γ, λ, τ, ω, σ, φ, μ (or in words: alpha, beta, lambda, tau, omega)

**NEVER (causes PDF warnings):**
- Emoji: all emoji (sparkles, charts, warnings, checkmarks, tools, folders, etc.)
- Box drawing characters: ┌ ─ └ │ ┐ ┘ ├ ┤
- Fancy quotes: „" '' (use standard " ')

**Rules for README.md:**
- README.md intended for PDF -> no emoji
- For emphasis use `**bold text**` or `## Headings`
- Advantage lists: use `+` instead of emoji checkmarks
- Section headings: without emoji (they're unnecessary, heading is descriptive enough)

**Examples:**
```markdown
# GOOD (PDF compatible)
## Key Features

**Advantages:**
+ Objective, reproducible results
+ Eliminates manual experimentation

## Parameters
- Parameter λ (lambda) for regularization
- Impedance Z = R + jωL [Ω]
- Phase φ = -45° for Warburg

# BAD (PDF warnings)
## 📋 Key Features

**Advantages:**
- ✅ Objective results
- ✅ Eliminates manual work

## 🔧 Parameters
- Parameter λ 🔧 for regularization
```

## Version Management

**Single Source of Truth:** Version is defined in `eis_analysis/version.py` - all other files import from there or synchronize manually.

**When changing version:** Always start from `version.py`, then update CHANGELOG.md and documentation.

**Details:** See `VERSION_MANAGEMENT.md`

## Development Workflow

### Adding a New Feature

1. **Planning**
   - Is it backward compatible?
   - Does it need a new module or extension of existing one?
   - What tests are needed?

2. **Implementation**
   - Choose the correct module in core package
   - Follow existing API conventions
   - Add docstrings (NumPy style)
   - Add type hints

3. **Testing**
   - Test with synthetic data
   - Test with real data (if available)
   - Test edge cases

4. **Documentation**
   - Update README.md (if user-facing)
   - Update CHANGELOG.md
   - Create specialized MD (if complex)

5. **Integration**
   - Update CLI script if needed
   - Update package exports (`__init__.py`)
   - Test entire workflow end-to-end

### Git Workflow

**NEVER:**
- Commit to main without testing
- Force push to main
- Commit without message
- Commit large binary files

**ALWAYS:**
- Descriptive commit messages
- Test before committing
- Update CHANGELOG.md
- One commit = one logical change

### Code Style

- **Python:** PEP 8 (but not strictly - longer lines OK if more readable)
- **Docstrings:** NumPy style
- **Type hints:** Yes, but not dogmatically
- **Comments:** Explain "why", not "what"
- **Naming:** Descriptive names (avoid abbreviations except common ones: freq, Z, R, C)
- **Linting:** `ruff check` for static analysis (unused imports, PEP 8, potential bugs)
- **Type checking:** `mypy` for type annotation checking

## Types of Development Tasks

### Adding a New Feature

**Consider:**
- Which architecture layer is appropriate? (I/O, algorithms, visualization, domain logic)
- Is a new module needed or extension of existing one?
- How does the feature fit into CLI workflow?
- Does it need specialized documentation?

**Proceed:**
1. Implement in appropriate module
2. Add to package exports if public API
3. Integrate into CLI if user-facing
4. Test with synthetic and real data
5. Document (docstring + README/specialized MD)
6. Update CHANGELOG.md

### Bug Fix

**Proceed:**
1. Identify root cause
2. Write test reproducing the bug (regression test)
3. Fix in appropriate module
4. Verify all tests pass
5. Update CHANGELOG.md (Bugfixes section)
6. Commit message: `fix(module): brief description`

### Adding Utility Script (data generators, etc.)

**Proceed:**
1. Standalone script at root level (not in package)
2. Argparse CLI interface
3. Export to standard formats
4. Optional `--plot` flag for quick preview
5. Create specialized MD documentation
6. Reference from README.md
7. Update CHANGELOG.md

## Notes for AI Assistant

### Environment

- Use `python3` command (not `python`) - system has no `python` symlink

### When User Requests a Change

1. **Determine scope:**
   - Is it a bugfix, feature, or refactoring?
   - Which modules will it affect?
   - Is it a breaking change?

2. **Propose solution:**
   - Where exactly to change code
   - What tests to add
   - What to update in docs

3. **Implement gradually:**
   - First core logic
   - Then integration
   - Then documentation
   - Then tests

4. **Verify:**
   - Run script on synthetic data
   - Check backward compatibility
   - Verify docs are up to date

### Prioritization

**High priority:**
- Bugfixes affecting results
- Missing/incorrect documentation
- Breaking API changes (requires major version bump)

**Medium priority:**
- Performance improvements
- New features (if user request)
- Code refactoring

**Low priority:**
- Cosmetic changes
- Alternative implementations (if current one works)

### Anti-patterns to Avoid

**Code:**
- Adding features "just in case" -> YAGNI (You Ain't Gonna Need It)
- Premature optimization -> first correct, then fast
- Ignoring edge cases -> always think about min/max/empty/NaN

**Documentation:**
- Copying docstring to README -> reference the module
- Duplicating information -> single source of truth
- Writing documentation-for-documentation-sake -> every section must have purpose

**Testing:**
- "It looks correct" -> always quantify (error < X%)
- Test only happy path -> test edge cases
- Manual testing -> prefer automated

---

**This document defines development philosophy and principles. For implementation details see README.md and other documentation.**
