"""Test configuration for pytest."""

import pytest
import sys


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "rdkit: tests that require RDKit to be installed"
    )
    config.addinivalue_line(
        "markers", "mongodb: tests that require MongoDB connection"
    )
    config.addinivalue_line(
        "markers", "integration: integration tests requiring all dependencies"
    )
    config.addinivalue_line(
        "markers", "slow: slow-running tests"
    )


def pytest_collection_modifyitems(config, items):
    """Auto-skip tests based on available dependencies."""
    # Check for RDKit
    rdkit_available = False
    try:
        import rdkit
        rdkit_available = True
    except ImportError:
        pass
    
    # Check for pymongo
    pymongo_available = False
    try:
        import pymongo
        pymongo_available = True
    except ImportError:
        pass
    
    skip_rdkit = pytest.mark.skip(reason="RDKit not installed")
    skip_mongodb = pytest.mark.skip(reason="pymongo not installed")
    
    for item in items:
        if "rdkit" in item.keywords and not rdkit_available:
            item.add_marker(skip_rdkit)
        if "mongodb" in item.keywords and not pymongo_available:
            item.add_marker(skip_mongodb)


@pytest.fixture
def sample_smiles():
    """Return a list of sample SMILES strings for testing."""
    return [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    ]


@pytest.fixture
def sample_fingerprints():
    """Return sample fingerprint data for testing."""
    return {
        'aspirin': [10, 25, 42, 101, 256, 512, 1024],
        'caffeine': [15, 30, 45, 102, 257, 513, 1025],
        'ibuprofen': [20, 35, 50, 103, 258, 514, 1026],
    }


@pytest.fixture
def mock_mongodb_collection():
    """Create a mock MongoDB collection for testing."""
    from unittest.mock import Mock
    
    collection = Mock()
    collection.name = 'test_collection'
    collection.count_documents.return_value = 0
    collection.find.return_value = []
    collection.insert_one.return_value = Mock(inserted_id='test_id')
    collection.insert_many.return_value = Mock(inserted_ids=['id1', 'id2'])
    
    return collection
