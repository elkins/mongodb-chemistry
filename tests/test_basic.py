"""Basic tests for mchem package."""

import pytest
from mchem import __version__


def test_version():
    """Test that version is defined."""
    assert __version__ == "0.1.0"


def test_import():
    """Test that main modules can be imported."""
    from mchem import similarity
    
    # These require RDKit, so only test if available
    try:
        from mchem import build, fps, cli
        assert build is not None
        assert fps is not None
        assert cli is not None
    except ImportError:
        pytest.skip("RDKit not installed")
    
    assert similarity is not None
