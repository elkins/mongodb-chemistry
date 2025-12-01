# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-12-01

### Major Modernization Release

This release represents a complete modernization of the mongodb-chemistry project for 2025.

### Added

- Modern `pyproject.toml` configuration (PEP 517/518 compliant)
- GitHub Actions CI/CD workflow for automated testing
- Pre-commit hooks configuration for code quality
- Basic pytest test suite with coverage reporting
- `.gitignore` file for common Python/IDE artifacts
- Type checking setup with mypy
- Code formatting with Black
- Linting with Ruff
- Comprehensive updated README with modern usage examples

### Changed

- **Python Version**: Now requires Python 3.8+ (dropped Python 2.7 support)
- **Dependencies**: Updated all dependencies to current versions:
  - `click`: 3.3 → 8.0+
  - `numpy`: 1.8.1 → 1.21+
  - `scipy`: 0.14.0 → 1.7+
  - `matplotlib`: 1.4.2 → 3.5+
  - `pymongo`: 2.7.2 → 4.0+ (major API changes)
  - `psycopg2`: 2.5.4 → 2.9+ (now using psycopg2-binary)
  - Replaced deprecated `nose` with `pytest`
- **Build System**: Migrated from `setup.py` to modern `pyproject.toml`
- **Version**: 0.0.1 → 0.1.0

### Fixed

- Updated MongoDB API calls for pymongo 4.x compatibility:
  - `.count()` → `.count_documents({})`
  - `.insert()` → `.insert_one()`
  - `.ensure_index()` → `.create_index()`
  - `timeout=False` → `no_cursor_timeout=True`
  - `.aggregate()` now returns cursor, wrap with `list()`
- Removed deprecated `urllib2`, now using `urllib.request`
- Fixed file opening for gzip files in text mode
- Removed unnecessary `.encode()` calls for SMILES strings
- Updated string formatting from `%` to f-strings throughout

### Removed

- Python 2 compatibility code (`from __future__ import` statements)
- `# -*- coding: utf-8 -*-` markers (Python 3 default)
- Deprecated `setup.py` file (replaced by `pyproject.toml`)
- Old requirements pinning that prevented updates

### Security

- Updated all dependencies to versions without known vulnerabilities
- Added security scanning via pre-commit hooks

### Development

- Added modern development tooling configuration
- Standardized code style with Black (100 char line length)
- Improved code quality with Ruff linter
- Added type checking foundation with mypy
- CI/CD pipeline tests on Python 3.8-3.12 across Ubuntu and macOS

## [0.0.1] - 2014

### Initial Release

- Basic MongoDB chemistry database functionality
- Chemical similarity search implementation
- Morgan fingerprint generation
- Performance profiling and analysis tools
- Command-line interface
- PostgreSQL comparison tools

---

**Migration Notes for Users:**

If upgrading from the original version:

1. **Python 3.8+ Required**: Update your Python environment
2. **MongoDB 4.0+**: Ensure MongoDB is up to date
3. **Reinstall**: `pip install -e .` or `pip install -e ".[dev]"` for development
4. **API Changes**: If you've written custom code using this library, check for:
   - Updated pymongo API calls
   - Removed Python 2 compatibility imports
   - String handling changes (no more `.encode()` needed in most cases)

[0.1.0]: https://github.com/elkins/mongodb-chemistry/compare/v0.0.1...v0.1.0
[0.0.1]: https://github.com/elkins/mongodb-chemistry/releases/tag/v0.0.1
