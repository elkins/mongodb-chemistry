# Contributing to MongoDB Chemistry

Thank you for your interest in contributing to mongodb-chemistry! This document provides guidelines and instructions for contributing.

## Development Setup

### Prerequisites

- Python 3.8 or higher
- MongoDB 4.0 or higher
- RDKit (installed via conda or system package manager)
- Git

### Setting Up Your Environment

1. **Fork and clone the repository**

```bash
git clone https://github.com/YOUR_USERNAME/mongodb-chemistry.git
cd mongodb-chemistry
```

2. **Create a virtual environment**

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install in development mode**

```bash
pip install -e ".[dev]"
```

4. **Install pre-commit hooks**

```bash
pre-commit install
```

## Development Workflow

### Before Making Changes

1. Create a new branch for your feature or fix:

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/issue-description
```

2. Ensure all tests pass:

```bash
pytest
```

### Making Changes

1. **Write code** following the project style:
   - Use Python 3.8+ features
   - Follow PEP 8 (enforced by Black and Ruff)
   - Add type hints where appropriate
   - Use f-strings for string formatting
   - Keep functions focused and documented

2. **Write tests** for new functionality:
   - Add tests in the `tests/` directory
   - Aim for good test coverage
   - Test edge cases and error conditions

3. **Update documentation**:
   - Add docstrings to new functions/classes
   - Update README.md if adding new features
   - Update CHANGELOG.md with your changes

### Code Quality Checks

The project uses several tools to maintain code quality:

```bash
# Format code with Black
black .

# Check and fix linting issues
ruff check . --fix

# Type checking (optional but encouraged)
mypy mchem

# Run all pre-commit hooks manually
pre-commit run --all-files
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=mchem --cov-report=html

# Run specific test file
pytest tests/test_basic.py -v

# Run tests matching a pattern
pytest -k "test_fingerprint"
```

### Committing Changes

The pre-commit hooks will automatically run when you commit:

```bash
git add .
git commit -m "feat: add new fingerprint type support"
```

**Commit Message Format:**

Use conventional commit format:
- `feat:` - New feature
- `fix:` - Bug fix
- `docs:` - Documentation changes
- `style:` - Code style changes (formatting, etc.)
- `refactor:` - Code refactoring
- `test:` - Adding or updating tests
- `chore:` - Maintenance tasks

Example:
```
feat: add support for ECFP fingerprints

- Implement ECFP fingerprint generation
- Add CLI option for ECFP selection
- Add tests for ECFP functionality
```

### Submitting a Pull Request

1. **Push your branch** to your fork:

```bash
git push origin feature/your-feature-name
```

2. **Create a Pull Request** on GitHub:
   - Go to the original repository
   - Click "New Pull Request"
   - Select your branch
   - Fill in the PR template with:
     - Description of changes
     - Related issue numbers (if any)
     - Testing performed
     - Breaking changes (if any)

3. **Address review feedback**:
   - Make requested changes
   - Push additional commits
   - Respond to comments

## Testing Guidelines

### Writing Good Tests

- **Test one thing**: Each test should verify one specific behavior
- **Use descriptive names**: Test names should clearly indicate what they test
- **Arrange-Act-Assert**: Structure tests in three clear sections
- **Use fixtures**: Share common setup using pytest fixtures

Example:

```python
def test_morgan_fingerprint_generation(sample_smiles):
    """Test that Morgan fingerprints are generated correctly."""
    # Arrange
    fingerprinter = MorganFingerprinter(radius=2)
    mol = Chem.MolFromSmiles(sample_smiles[0])
    
    # Act
    fp = fingerprinter.generate(mol)
    
    # Assert
    assert isinstance(fp, list)
    assert len(fp) > 0
    assert all(isinstance(bit, int) for bit in fp)
```

### Test Coverage

- Aim for >80% coverage for new code
- Critical functions should have comprehensive tests
- Don't sacrifice test quality for coverage numbers

## Documentation

### Docstring Format

Use Google-style docstrings:

```python
def similarity_search(mol, fingerprinter, fp_collection, threshold=0.8, count_collection=None):
    """Perform a similarity search using aggregation framework.

    Args:
        mol: The query molecule (RDKit Mol object).
        fingerprinter: The Fingerprinter instance to use.
        fp_collection: MongoDB fingerprint collection to query.
        threshold: The Tanimoto similarity threshold (default: 0.8).
        count_collection: Optional MongoDB collection containing fingerprint bit frequencies.

    Returns:
        List of dictionaries containing matched molecule IDs and similarity scores.

    Raises:
        ValueError: If threshold is not between 0 and 1.
        
    Example:
        >>> mol = Chem.MolFromSmiles("CC(=O)O")
        >>> fingerprinter = MorganFingerprinter(radius=2)
        >>> results = similarity_search(mol, fingerprinter, fp_coll, threshold=0.7)
    """
```

## Code Style Guidelines

### General Principles

- **Readability counts**: Clear code is better than clever code
- **Explicit is better than implicit**: Don't hide functionality
- **Simple is better than complex**: Avoid unnecessary complexity

### Specific Guidelines

- **Line length**: Maximum 100 characters (enforced by Black)
- **Imports**: Organized by standard library, third-party, local
- **Naming conventions**:
  - Functions/variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_CASE`
- **Type hints**: Add to function signatures when practical
- **Error handling**: Use specific exceptions, include helpful messages

### Python 3.8+ Features

Feel free to use modern Python features:
- f-strings for formatting
- Type hints (from `typing` module)
- Dataclasses where appropriate
- Pathlib for file operations
- Assignment expressions (`:=`) when they improve readability

## Getting Help

- **Questions**: Open a GitHub Discussion or issue
- **Bugs**: Open a GitHub Issue with reproduction steps
- **Security**: Email maintainers privately for security issues

## Recognition

Contributors will be acknowledged in:
- CHANGELOG.md for their contributions
- GitHub contributors list
- Release notes

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Assume good intentions

Thank you for contributing to mongodb-chemistry! 🎉
