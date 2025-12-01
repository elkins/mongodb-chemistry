.PHONY: help install install-dev test lint format clean run-tests smoke-test

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install package
	pip install -e .

install-dev:  ## Install package with dev dependencies
	pip install -e ".[dev]"
	pre-commit install

test:  ## Run tests with pytest
	pytest -v

test-cov:  ## Run tests with coverage
	pytest -v --cov=mchem --cov-report=html --cov-report=term

lint:  ## Run linting checks
	ruff check .
	black --check .

format:  ## Format code with black and ruff
	black .
	ruff check . --fix

type-check:  ## Run type checking with mypy
	mypy mchem || true

smoke-test:  ## Run syntax validation on all Python files
	@python -m compileall -q mchem/ tests/
	@echo "✓ Smoke test passed - all Python files are valid"

clean:  ## Clean up build artifacts and cache
	rm -rf build/ dist/ *.egg-info .pytest_cache .coverage htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete

pre-commit:  ## Run pre-commit hooks on all files
	pre-commit run --all-files

run-tests: lint test  ## Run linting and tests

all: format lint test  ## Format, lint, and test

.DEFAULT_GOAL := help
