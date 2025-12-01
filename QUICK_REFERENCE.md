# Quick Reference - MongoDB Chemistry v0.1.0

## Installation

```bash
# Install with pip
pip install git+https://github.com/elkins/mongodb-chemistry.git

# Development install
git clone https://github.com/elkins/mongodb-chemistry.git
cd mongodb-chemistry
pip install -e ".[dev]"
pre-commit install
```

## CLI Quick Commands

```bash
# Basic workflow
mchem load molecules.sdf              # Load molecules
mchem addfp                           # Generate fingerprints
mchem countfp                         # Count bit frequencies
mchem similar "CC(=O)O" --threshold 0.8  # Search

# With options
mchem --db mydb --uri mongodb://localhost:27017 load data.sdf
mchem addfp --fp morgan --radius 3 --length 2048
mchem similar "SMILES" --threshold 0.7 --collection compounds

# Benchmarking
mchem sample --size 1000 --output sample.txt
mchem benchmark --sample sample.txt
mchem analyse screening --sample sample.txt
```

## Python API

### Basic Usage

```python
import pymongo
from rdkit import Chem
from mchem import fps, similarity

# Connect to MongoDB
client = pymongo.MongoClient('mongodb://localhost:27017')
db = client['mchem']

# Create fingerprinter
fingerprinter = fps.MorganFingerprinter(radius=2)

# Generate fingerprint
mol = Chem.MolFromSmiles('CC(=O)O')
fp = fingerprinter.generate(mol)

# Similarity search
fp_collection = db['mols.m2']
count_collection = db['mols.m2.counts']
results = similarity.similarity_search(
    mol, 
    fingerprinter, 
    fp_collection, 
    threshold=0.8,
    count_collection=count_collection
)
```

### Loading Data

```python
from mchem import build

# Load SDF files
sdf_paths = ['molecules.sdf', 'more_molecules.sdf.gz']
collection = db['molecules']
build.load_sdfs(sdf_paths, collection, idfield='CHEMBL_ID')
```

### Generating Fingerprints

```python
from mchem import fps

# Different fingerprint types
fp1 = fps.MorganFingerprinter(radius=2)          # Unfolded
fp2 = fps.MorganFingerprinter(radius=2, length=1024)  # Folded
fp3 = fps.MorganFingerprinter(radius=3)          # Different radius

# Generate for collection
mol_collection = db['molecules']
fp_collection = db['molecules.m2']
fps.generate(mol_collection, fp_collection, fp1)

# Count bit frequencies
count_collection = db['molecules.m2.counts']
fps.count(fp_collection, count_collection)
```

### Similarity Searching

```python
from mchem import similarity

# Search with molecule
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
results = similarity.similarity_search(
    mol, 
    fingerprinter, 
    fp_collection, 
    threshold=0.8,
    count_collection=count_collection
)

# Search with pre-computed fingerprint
qfp = fingerprinter.generate(mol)
results = similarity.similarity_search_fp(
    qfp, 
    fp_collection, 
    threshold=0.8,
    count_collection=count_collection
)

# Process results
for result in results:
    mol_id = result['_id']
    tanimoto = result['tanimoto']
    print(f"{mol_id}: {tanimoto:.3f}")
```

## MongoDB Queries

### Direct MongoDB Access

```python
# Count documents
count = collection.count_documents({})
count = collection.count_documents({'field': 'value'})

# Find molecules
mol_doc = collection.find_one({'_id': 'CHEMBL123'})
mol = Chem.Mol(mol_doc['rdmol'])

# Find with criteria
cursor = collection.find({'count': {'$gte': 50}})
for doc in cursor:
    process(doc)

# Aggregation
pipeline = [
    {'$match': {'count': {'$gte': 10}}},
    {'$limit': 100}
]
results = list(collection.aggregate(pipeline))
```

## Development Workflow

```bash
# Format code
black .

# Lint code
ruff check .
ruff check . --fix  # Auto-fix

# Type check
mypy mchem

# Run tests
pytest                    # All tests
pytest -v                 # Verbose
pytest --cov=mchem       # With coverage
pytest tests/test_basic.py -k test_version  # Specific test

# Pre-commit
pre-commit run --all-files
```

## Common Patterns

### Batch Processing

```python
from tqdm import tqdm  # Optional: progress bar

# Process molecules in batches
batch_size = 1000
cursor = collection.find(no_cursor_timeout=True)

for i, mol_doc in enumerate(tqdm(cursor)):
    mol = Chem.Mol(mol_doc['rdmol'])
    # Process molecule
    
    if i % batch_size == 0:
        # Commit/checkpoint every batch
        pass
```

### Error Handling

```python
from pymongo.errors import DuplicateKeyError, ConnectionFailure
from rdkit import Chem

try:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    collection.insert_one({'_id': mol_id, 'rdmol': mol.ToBinary()})
    
except DuplicateKeyError:
    print(f"Molecule {mol_id} already exists")
except ConnectionFailure:
    print("MongoDB connection failed")
except Exception as e:
    print(f"Error: {e}")
```

### Custom Fingerprints

```python
from mchem.fps import Fingerprinter

class CustomFingerprinter(Fingerprinter):
    def __init__(self, param1, param2):
        self.param1 = param1
        self.param2 = param2
    
    def generate(self, mol):
        # Your custom fingerprint logic
        # Must return sorted list of bit indices
        return sorted([1, 5, 10, 25, 50])
    
    @property
    def name(self):
        return f'custom_{self.param1}_{self.param2}'
```

## Environment Variables

```bash
# MongoDB connection
export MCHEM_MONGODB_URI="mongodb://localhost:27017"
export MCHEM_MONGODB_DB="mchem"
export MCHEM_MONGODB_COLL="mols"

# PostgreSQL (for pgchem)
export MCHEM_POSTGRES_DB="mchem"
export MCHEM_POSTGRES_USER="username"
export MCHEM_POSTGRES_PASSWORD="password"
```

## Database Schema

### Molecule Collection
```javascript
{
    "_id": "CHEMBL123",           // Molecule identifier
    "smiles": "CC(=O)O",          // SMILES string
    "rdmol": BinData(...)         // RDKit binary molecule
}
```

### Fingerprint Collection
```javascript
{
    "_id": "CHEMBL123",           // Matches molecule _id
    "bits": [1, 5, 10, 25, ...],  // Sorted list of set bits
    "count": 234                  // Number of set bits
}
```

### Count Collection
```javascript
{
    "_id": 42,                    // Bit index
    "count": 1523                 // Number of molecules with this bit
}
```

## Performance Tips

1. **Always create indexes** on `bits` and `count` fields
2. **Use count collection** for faster searches (pre-calculate with `mchem countfp`)
3. **Adjust threshold** - higher threshold = faster search
4. **Fold fingerprints** for memory efficiency (at cost of some accuracy)
5. **Use projection** to return only needed fields
6. **Batch operations** instead of individual inserts

## Troubleshooting

```bash
# Check installation
python -c "import mchem; print(mchem.__version__)"

# Check RDKit
python -c "from rdkit import Chem; print(Chem.__version__)"

# Check MongoDB connection
python -c "import pymongo; c = pymongo.MongoClient(); print(c.server_info())"

# Check database contents
mongo mchem --eval "db.mols.countDocuments({})"

# View logs with verbose mode
mchem --verbose load molecules.sdf
```

## Useful MongoDB Commands

```bash
# In mongo shell
use mchem

# List collections
show collections

# Count documents
db.mols.countDocuments({})

# View sample document
db.mols.findOne()

# Check indexes
db.mols.m2.getIndexes()

# Database stats
db.stats()

# Collection stats
db.mols.stats()
```

## Links

- **GitHub**: https://github.com/elkins/mongodb-chemistry
- **Issues**: https://github.com/elkins/mongodb-chemistry/issues
- **RDKit Docs**: https://www.rdkit.org/docs/
- **PyMongo Docs**: https://pymongo.readthedocs.io/
- **MongoDB Docs**: https://docs.mongodb.com/

---

**Version**: 0.1.0  
**Updated**: December 2025  
**Python**: 3.8+  
**MongoDB**: 4.0+
