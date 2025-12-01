# Upgrade Guide: v0.0.1 → v0.1.0

This guide helps you upgrade from the original mongodb-chemistry (Python 2.7, circa 2014) to the modernized version (Python 3.8+, 2025).

## Quick Start

If you're in a hurry:

```bash
# 1. Check Python version
python --version  # Should be 3.8 or higher

# 2. Reinstall
pip uninstall mongodb-chemistry
pip install git+https://github.com/elkins/mongodb-chemistry.git

# 3. Test
mchem --version
```

## Prerequisites

### 1. Python 3.8+

Check your Python version:

```bash
python --version
```

If you need to upgrade:

**macOS (using Homebrew):**
```bash
brew install python@3.11
```

**Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install python3.11 python3.11-venv
```

**Using pyenv (recommended for managing multiple versions):**
```bash
curl https://pyenv.run | bash
pyenv install 3.11.6
pyenv global 3.11.6
```

### 2. MongoDB 4.0+

Check your MongoDB version:

```bash
mongod --version
```

Upgrade if needed:

**macOS:**
```bash
brew upgrade mongodb-community
```

**Ubuntu:**
```bash
# Follow official MongoDB upgrade guide
# https://docs.mongodb.com/manual/tutorial/upgrade-revision/
```

### 3. RDKit

RDKit installation hasn't changed:

**macOS:**
```bash
brew install rdkit
```

**Ubuntu:**
```bash
sudo apt-get install python3-rdkit
```

**Conda (any platform):**
```bash
conda install -c conda-forge rdkit
```

## Step-by-Step Upgrade

### Step 1: Backup Your Data

**Critical:** Backup your MongoDB database before upgrading:

```bash
# Backup entire database
mongodump --db mchem --out /path/to/backup

# Or export specific collection
mongoexport --db mchem --collection mols --out mols_backup.json
```

### Step 2: Create New Virtual Environment

It's recommended to start fresh:

```bash
# Deactivate old environment if active
deactivate

# Create new environment
python3.11 -m venv venv-mchem-new
source venv-mchem-new/bin/activate

# Verify Python version in venv
python --version
```

### Step 3: Install Modernized Version

```bash
# Uninstall old version
pip uninstall mongodb-chemistry -y

# Install new version
pip install git+https://github.com/elkins/mongodb-chemistry.git

# Verify installation
mchem --version  # Should show v0.1.0
```

### Step 4: Test Basic Functionality

```bash
# Test MongoDB connection
mchem --help

# If you have existing data, try a simple query
mchem similar "CC(=O)O" --threshold 0.8
```

## Code Changes

### If You Have Custom Scripts

If you've written Python scripts that use mongodb-chemistry, you'll need to update them:

#### 1. Remove Python 2 Compatibility

**Before:**
```python
from __future__ import print_function, unicode_literals, division
```

**After:**
```python
# Simply remove these lines - Python 3 default behavior
```

#### 2. Update String Handling

**Before:**
```python
mol = Chem.MolFromSmiles(smiles.encode())
idfield = rdmol.GetProp(idfield.encode())
```

**After:**
```python
mol = Chem.MolFromSmiles(smiles)  # No .encode() needed
idfield = rdmol.GetProp(idfield)  # No .encode() needed
```

#### 3. Update MongoDB API Calls

**Before:**
```python
# Counting documents
count = collection.count()

# Inserting documents
collection.insert(document)

# Creating indexes
collection.ensure_index('field_name')

# Find with no timeout
cursor = collection.find(timeout=False)

# Aggregate results
results = collection.aggregate(pipeline)['result']
```

**After:**
```python
# Counting documents
count = collection.count_documents({})  # Empty filter for all

# Inserting documents
collection.insert_one(document)  # or insert_many([docs])

# Creating indexes
collection.create_index('field_name')

# Find with no timeout
cursor = collection.find(no_cursor_timeout=True)

# Aggregate results
results = list(collection.aggregate(pipeline))  # Returns cursor
```

#### 4. Update String Formatting

**Before:**
```python
log.info('Processing %s of %s items' % (current, total))
message = 'Result: %s' % result
```

**After:**
```python
log.info(f'Processing {current} of {total} items')
message = f'Result: {result}'
```

#### 5. Update File Operations

**Before:**
```python
if filename[-3:] == '.gz':
    f = gzip.open(filename)
```

**After:**
```python
if filename.endswith('.gz'):
    f = gzip.open(filename, 'rt')  # Text mode
```

## Testing Your Upgrade

### 1. Test CLI Commands

```bash
# Load test data
mchem load data/sample_chembl_1000.txt

# Generate fingerprints
mchem addfp --fp morgan --radius 2

# Count fingerprints
mchem countfp

# Test similarity search (Aspirin)
mchem similar "CC(=O)Oc1ccccc1C(=O)O" --threshold 0.7
```

### 2. Test Your Custom Scripts

Create a simple test script:

```python
# test_upgrade.py
import pymongo
from rdkit import Chem
from mchem import fps, similarity

# Test connection
client = pymongo.MongoClient('mongodb://localhost:27017')
db = client['mchem']

# Test fingerprint generation
fingerprinter = fps.MorganFingerprinter(radius=2)
mol = Chem.MolFromSmiles('CC(=O)O')
fp = fingerprinter.generate(mol)
print(f"Fingerprint generated: {len(fp)} bits")

# Test similarity search
fp_collection = db['mols.m2']
results = similarity.similarity_search(
    mol, 
    fingerprinter, 
    fp_collection, 
    threshold=0.8
)
print(f"Found {len(results)} similar molecules")
```

Run it:
```bash
python test_upgrade.py
```

## Common Issues and Solutions

### Issue 1: Import Errors

**Error:**
```
ImportError: cannot import name 'X' from 'mchem'
```

**Solution:**
Ensure you've uninstalled the old version and installed the new one:
```bash
pip uninstall mongodb-chemistry -y
pip install --force-reinstall git+https://github.com/elkins/mongodb-chemistry.git
```

### Issue 2: MongoDB Connection Errors

**Error:**
```
pymongo.errors.ServerSelectionTimeoutError
```

**Solution:**
1. Ensure MongoDB is running: `brew services start mongodb-community` (macOS)
2. Check MongoDB version: `mongod --version`
3. Try explicit URI: `mchem --uri mongodb://localhost:27017 --db mchem`

### Issue 3: RDKit Import Errors

**Error:**
```
ModuleNotFoundError: No module named 'rdkit'
```

**Solution:**
RDKit must be installed separately:
```bash
# Using conda (recommended)
conda install -c conda-forge rdkit

# Or system package
brew install rdkit  # macOS
sudo apt-get install python3-rdkit  # Ubuntu
```

### Issue 4: Aggregate Query Returns Empty

**Error:**
Empty results from similarity searches that worked before.

**Solution:**
The aggregate method now returns a cursor. Wrap it with `list()`:
```python
# In your custom code
results = list(collection.aggregate(pipeline))
```

This is already fixed in the modernized library.

### Issue 5: Slow Performance

**Issue:**
Queries seem slower after upgrade.

**Investigation:**
```bash
# Check indexes
mongo mchem --eval "db.mols.m2.getIndexes()"

# Rebuild indexes if needed
mchem addfp --fp morgan --radius 2  # Regenerates with new indexes
```

## Rollback Plan

If you encounter critical issues:

### 1. Restore Old Environment

```bash
# Keep old environment around temporarily
source venv-mchem-old/bin/activate
# Use old version
```

### 2. Restore Data

```bash
# Restore from backup
mongorestore --db mchem /path/to/backup/mchem
```

### 3. Report Issues

Open an issue on GitHub with:
- Python version
- MongoDB version
- Error messages
- Steps to reproduce

## Performance Tuning

After upgrading, optimize for the new version:

### 1. Rebuild Indexes

```bash
# Drop and recreate fingerprint collections
mongo mchem --eval "db.mols.m2.drop()"
mchem addfp --fp morgan --radius 2
mchem countfp
```

### 2. Update MongoDB Configuration

For MongoDB 4.0+, consider:
- Enabling WiredTiger compression
- Adjusting cache size
- Monitoring with MongoDB Compass

### 3. Profile Queries

```python
# Add timing to your scripts
import time

start = time.time()
results = similarity.similarity_search(...)
elapsed = time.time() - start
print(f"Query took {elapsed:.2f}s")
```

## Getting Help

- **Documentation**: See updated README.md
- **Issues**: https://github.com/elkins/mongodb-chemistry/issues
- **Discussions**: GitHub Discussions
- **Original docs**: Blog posts still relevant for concepts

## Verification Checklist

After upgrading, verify:

- [ ] Python version is 3.8+
- [ ] MongoDB version is 4.0+
- [ ] RDKit is installed and importable
- [ ] `mchem --version` shows v0.1.0
- [ ] Can connect to MongoDB
- [ ] Can load SDF files
- [ ] Can generate fingerprints
- [ ] Can perform similarity searches
- [ ] Custom scripts work (if applicable)
- [ ] Performance is acceptable

## Next Steps

Once upgraded successfully:

1. **Update documentation** of your own scripts
2. **Consider using new features** like:
   - Pre-commit hooks for code quality
   - pytest for testing
   - Type hints for better IDE support
3. **Contribute back** if you find issues or improvements
4. **Stay updated** by watching the repository

## Timeline

- **Immediate**: Can use old version in old environment
- **Month 1-2**: Test new version in parallel
- **Month 3**: Fully migrate to new version
- **Month 6+**: Old Python 2.7 version unsupported

---

**Need Help?** Open an issue on GitHub with the "upgrade" label.

**Updated:** December 1, 2025
