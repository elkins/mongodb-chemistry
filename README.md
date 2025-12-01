# MongoDB Chemistry

Chemical similarity search implementation in MongoDB, with performance analysis.

[![CI](https://github.com/elkins/mongodb-chemistry/workflows/CI/badge.svg)](https://github.com/elkins/mongodb-chemistry/actions)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Modernized for Python 3.8+ with updated dependencies and modern build system.**

See [this blog post](http://blog.matt-swain.com/post/87093745652/chemical-similarity-search-in-mongodb) for more information.

See also:

- [Presentation slides](https://speakerdeck.com/mcs07/chemical-structure-handling-in-mongodb)
- [Blog post by Rajarshi Guha](http://blog.rguha.net/?p=1261)
- [Blog post by Davy Suvee](http://datablend.be/?p=254) ([part 2](http://datablend.be/?p=256), [part 3](http://datablend.be/?p=265))

## Features

- Chemical similarity search using MongoDB's aggregation framework
- Support for Morgan fingerprints with customizable radius and folding
- Performance benchmarking and analysis tools
- Command-line interface for database operations
- PostgreSQL comparison tools

## Requirements

- Python 3.8 or higher
- MongoDB 4.0 or higher
- RDKit

## Installation

### Using pip (recommended)

```bash
# Install RDKit first (system package)
# macOS:
brew install rdkit mongodb-community

# Ubuntu/Debian:
sudo apt-get install python3-rdkit mongodb

# Install mongodb-chemistry
pip install git+https://github.com/elkins/mongodb-chemistry.git
```

### Development Installation

```bash
# Clone the repository
git clone https://github.com/elkins/mongodb-chemistry.git
cd mongodb-chemistry

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (optional)
pre-commit install
```

### Ubuntu Bootstrap

For a complete Ubuntu setup, see [scripts/bootstrap.sh](https://github.com/mcs07/mongodb-chemistry/blob/master/scripts/bootstrap.sh).

## Quick Start

```bash
# Start MongoDB
mongod

# Load molecules from SDF file
mchem load mymols.sdf

# Generate Morgan fingerprints (radius=2)
mchem addfp

# Pre-calculate bit frequencies for faster searching
mchem countfp

# Perform similarity search (Aspirin SMILES)
mchem similar "O=C(Oc1ccccc1C(=O)O)C" --threshold 0.8
```

## Usage Examples

### Basic Operations

```bash
# Load with custom collection name
mchem --db mydb load molecules.sdf --collection compounds

# Generate fingerprints with custom parameters
mchem addfp --fp morgan --radius 3 --length 2048

# Similarity search with custom threshold
mchem similar "CC(C)Cc1ccc(cc1)C(C)C(=O)O" --threshold 0.7
```

### Advanced Usage

```bash
# Create a random sample for benchmarking
mchem sample --size 1000 --output sample_ids.txt

# Run similarity search on sample molecules
mchem samplesim --sample sample_ids.txt --threshold 0.8

# Benchmark performance
mchem benchmark --sample sample_ids.txt

# Analyze screening methods
mchem analyse screening --sample sample_ids.txt
```

### Configuration

Set environment variables for default connection:

```bash
export MCHEM_MONGODB_URI="mongodb://localhost:27017"
export MCHEM_MONGODB_DB="mchem"
export MCHEM_MONGODB_COLL="mols"
```

## Development

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=mchem --cov-report=html

# Run specific test file
pytest tests/test_basic.py -v
```

### Code Quality

```bash
# Format code with black
black .

# Lint with ruff
ruff check .

# Type checking with mypy
mypy mchem
```

## Project Structure

```
mongodb-chemistry/
├── mchem/              # Main package
│   ├── build.py        # Database building functions
│   ├── fps.py          # Fingerprint generation
│   ├── similarity.py   # Similarity search algorithms
│   ├── cli.py          # Command-line interface
│   ├── screening.py    # Screening analysis
│   ├── profile.py      # Performance profiling
│   ├── plot.py         # Results visualization
│   └── postgres.py     # PostgreSQL comparison tools
├── tests/              # Test suite
├── scripts/            # Utility scripts
├── data/               # Sample data
└── pyproject.toml      # Modern Python project configuration
```

## Modernization Changes (2025)

This fork includes significant modernization:

- ✅ **Python 3.8+**: Removed Python 2 compatibility code
- ✅ **Modern Dependencies**: Updated all dependencies to current versions
- ✅ **Build System**: Migrated from `setup.py` to `pyproject.toml`
- ✅ **MongoDB 4.x**: Updated for pymongo 4.x API changes
- ✅ **Code Quality**: Added black, ruff, mypy, pre-commit hooks
- ✅ **CI/CD**: GitHub Actions workflow for automated testing
- ✅ **Testing**: Basic pytest setup with coverage reporting
- ✅ **Type Safety**: Prepared for gradual type hint adoption

## Performance

See [scripts/chembl.sh](https://github.com/mcs07/mongodb-chemistry/blob/master/scripts/chembl.sh) for benchmarking examples using the ChEMBL database.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Install pre-commit hooks: `pre-commit install`
4. Make your changes with tests
5. Ensure all tests pass: `pytest`
6. Submit a pull request

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```
Swain, M. (2014). Chemical Similarity Search in MongoDB.
https://github.com/mcs07/mongodb-chemistry
```

## Acknowledgments

- Original implementation by [Matt Swain](https://github.com/mcs07)
- Built with [RDKit](https://www.rdkit.org/) and [MongoDB](https://www.mongodb.com/)
- Inspired by work from Rajarshi Guha and Davy Suvee
