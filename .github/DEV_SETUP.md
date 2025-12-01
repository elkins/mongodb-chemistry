# MongoDB Chemistry Development Environment

This directory contains your local development environment for mongodb-chemistry.

## Quick Start

```bash
# Install dependencies
make install-dev

# Run tests
make test

# Format and lint code
make format
make lint
```

## Environment Variables

Create a `.env` file (not tracked by git) with your local settings:

```bash
MCHEM_MONGODB_URI=mongodb://localhost:27017
MCHEM_MONGODB_DB=mchem_dev
MCHEM_MONGODB_COLL=mols
```

## MongoDB Setup

```bash
# Start MongoDB (macOS with Homebrew)
brew services start mongodb-community

# Or run manually
mongod --dbpath ~/data/db
```

## Running Examples

```bash
# Load sample data
mchem load data/sample_chembl_1000.txt

# Generate fingerprints
mchem addfp

# Run similarity search
mchem similar "CC(=O)O" --threshold 0.8
```

## Troubleshooting

- **Import errors**: Run `make install-dev` to ensure all dependencies are installed
- **MongoDB connection**: Check that MongoDB is running on port 27017
- **RDKit issues**: Install via conda: `conda install -c conda-forge rdkit`

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed development guidelines.
