# Pre-Commit Checklist

Before committing your changes, verify:

## Code Quality
- [ ] All Python files pass syntax validation: `make smoke-test`
- [ ] Code is formatted: `make format`
- [ ] No linting errors: `make lint`
- [ ] Type checking passes (if applicable): `make type-check`

## Testing
- [ ] All tests pass: `make test`
- [ ] New features have tests
- [ ] Test coverage is maintained

## Documentation
- [ ] Updated README.md if adding features
- [ ] Updated CHANGELOG.md with changes
- [ ] Added docstrings to new functions/classes
- [ ] Updated type hints where applicable

## Git Hygiene
- [ ] Removed debug code and print statements
- [ ] No sensitive data (passwords, API keys, etc.)
- [ ] No large binary files
- [ ] Removed test/debug files

## Final Checks
- [ ] Commit message follows conventional format
- [ ] Branch is up to date with main/master
- [ ] No merge conflicts
- [ ] Pre-commit hooks pass: `make pre-commit`

## Quick Commands

```bash
# Run all quality checks
make all

# Clean and test
make clean test

# Format and verify
make format lint test
```

## Common Issues

**Import errors**: `make install-dev`  
**Test failures**: Check MongoDB is running  
**Linting errors**: `make format` to auto-fix  
**Pre-commit hook fails**: Fix the issue, then re-commit
