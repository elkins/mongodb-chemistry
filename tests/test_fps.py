"""Unit tests for mchem.fps module."""

from unittest.mock import Mock, patch

import pytest

# Only import if RDKit is available
pytest.importorskip("rdkit")

from mchem import fps


class TestMorganFingerprinter:
    """Tests for MorganFingerprinter class."""

    def test_init_default(self):
        """Test MorganFingerprinter initialization with defaults."""
        fp = fps.MorganFingerprinter()
        assert fp.radius == 2
        assert fp.length is None

    def test_init_custom_radius(self):
        """Test MorganFingerprinter initialization with custom radius."""
        fp = fps.MorganFingerprinter(radius=3)
        assert fp.radius == 3
        assert fp.length is None

    def test_init_with_length(self):
        """Test MorganFingerprinter initialization with folding length."""
        fp = fps.MorganFingerprinter(radius=2, length=1024)
        assert fp.radius == 2
        assert fp.length == 1024

    def test_name_unfolded(self):
        """Test fingerprinter name for unfolded fingerprint."""
        fp = fps.MorganFingerprinter(radius=2)
        assert fp.name == "m2"

    def test_name_folded(self):
        """Test fingerprinter name for folded fingerprint."""
        fp = fps.MorganFingerprinter(radius=2, length=512)
        assert fp.name == "m2l512"

    def test_name_different_radius(self):
        """Test fingerprinter name with different radius."""
        fp = fps.MorganFingerprinter(radius=3, length=2048)
        assert fp.name == "m3l2048"

    @pytest.mark.rdkit
    def test_generate_unfolded(self):
        """Test fingerprint generation without folding."""
        from rdkit import Chem

        fp = fps.MorganFingerprinter(radius=2)
        mol = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
        result = fp.generate(mol)

        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(bit, int) for bit in result)
        assert result == sorted(result)  # Should be sorted

    @pytest.mark.rdkit
    def test_generate_folded(self):
        """Test fingerprint generation with folding."""
        from rdkit import Chem

        fp = fps.MorganFingerprinter(radius=2, length=1024)
        mol = Chem.MolFromSmiles("CC(=O)O")
        result = fp.generate(mol)

        assert isinstance(result, list)
        assert len(result) > 0
        assert all(0 <= bit < 1024 for bit in result)
        assert result == sorted(result)


class TestFingerprinterInterface:
    """Tests for Fingerprinter base class."""

    def test_fingerprinter_generate_not_implemented(self):
        """Test that Fingerprinter.generate raises NotImplementedError."""
        fp = fps.Fingerprinter()
        with pytest.raises(NotImplementedError):
            fp.generate(None)

    def test_fingerprinter_name_not_implemented(self):
        """Test that Fingerprinter.name raises NotImplementedError."""
        fp = fps.Fingerprinter()
        with pytest.raises(NotImplementedError):
            _ = fp.name


class TestGetFingerprinter:
    """Tests for get_fingerprinter function."""

    def test_get_morgan_fingerprinter(self):
        """Test getting Morgan fingerprinter."""
        fp = fps.get_fingerprinter("morgan", radius=2)
        assert isinstance(fp, fps.MorganFingerprinter)
        assert fp.radius == 2

    def test_get_morgan_with_length(self):
        """Test getting Morgan fingerprinter with length."""
        fp = fps.get_fingerprinter("morgan", radius=3, length=2048)
        assert isinstance(fp, fps.MorganFingerprinter)
        assert fp.radius == 3
        assert fp.length == 2048

    def test_get_invalid_fingerprinter(self):
        """Test that invalid fingerprinter name raises KeyError."""
        with pytest.raises(KeyError):
            fps.get_fingerprinter("invalid", radius=2)


@pytest.mark.mongodb
class TestGenerate:
    """Tests for generate function (requires MongoDB)."""

    def test_generate_creates_fingerprints(self):
        """Test that generate creates fingerprints in collection."""
        # Mock collections
        mock_mol_collection = Mock()
        mock_mol_collection.name = "test_mols"
        mock_mol_collection.find.return_value = [{"_id": "mol1", "rdmol": b"mock_rdmol_data"}]

        mock_fp_collection = Mock()
        mock_fp_collection.name = "test_fps"

        # Mock fingerprinter
        mock_fingerprinter = Mock()
        mock_fingerprinter.name = "m2"
        mock_fingerprinter.generate.return_value = [1, 5, 10, 25]

        # Mock Chem.Mol
        with patch("mchem.fps.Chem.Mol"):
            fps.generate(mock_mol_collection, mock_fp_collection, mock_fingerprinter)

        # Verify insert_one was called
        assert mock_fp_collection.insert_one.called
        call_args = mock_fp_collection.insert_one.call_args[0][0]
        assert call_args["_id"] == "mol1"
        assert call_args["bits"] == [1, 5, 10, 25]
        assert call_args["count"] == 4


@pytest.mark.mongodb
class TestCount:
    """Tests for count function (requires MongoDB)."""

    def test_count_builds_frequency_collection(self):
        """Test that count builds bit frequency collection."""
        # Mock fingerprint collection
        mock_fp_collection = Mock()
        mock_fp_collection.find.return_value = [
            {"_id": "mol1", "bits": [1, 2, 3]},
            {"_id": "mol2", "bits": [2, 3, 4]},
            {"_id": "mol3", "bits": [1, 3, 5]},
        ]

        # Mock count collection
        mock_count_collection = Mock()
        mock_count_collection.name = "test_counts"

        fps.count(mock_fp_collection, mock_count_collection)

        # Verify drop was called
        assert mock_count_collection.drop.called

        # Verify insert_one was called for each bit
        assert mock_count_collection.insert_one.call_count == 5

        # Check that counts are correct
        inserted_docs = [call[0][0] for call in mock_count_collection.insert_one.call_args_list]
        bit_counts = {doc["_id"]: doc["count"] for doc in inserted_docs}

        assert bit_counts[1] == 2  # Bit 1 appears in mol1 and mol3
        assert bit_counts[2] == 2  # Bit 2 appears in mol1 and mol2
        assert bit_counts[3] == 3  # Bit 3 appears in all three
        assert bit_counts[4] == 1  # Bit 4 appears in mol2
        assert bit_counts[5] == 1  # Bit 5 appears in mol3


@pytest.mark.rdkit
class TestIntegrationWithRDKit:
    """Integration tests with actual RDKit molecules."""

    def test_aspirin_fingerprint(self):
        """Test fingerprint generation for aspirin."""
        from rdkit import Chem

        fp = fps.MorganFingerprinter(radius=2)
        aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = fp.generate(aspirin)

        assert len(result) > 10  # Aspirin should have many bits set
        assert all(isinstance(bit, int) for bit in result)

    def test_different_molecules_different_fingerprints(self):
        """Test that different molecules produce different fingerprints."""
        from rdkit import Chem

        fp = fps.MorganFingerprinter(radius=2)
        aspirin = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        caffeine = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

        fp_aspirin = set(fp.generate(aspirin))
        fp_caffeine = set(fp.generate(caffeine))

        # Fingerprints should be different
        assert fp_aspirin != fp_caffeine
        # But may have some bits in common
        assert len(fp_aspirin & fp_caffeine) >= 0

    def test_same_molecule_same_fingerprint(self):
        """Test that the same molecule always produces the same fingerprint."""
        from rdkit import Chem

        fp = fps.MorganFingerprinter(radius=2)
        mol = Chem.MolFromSmiles("CC(=O)O")

        fp1 = fp.generate(mol)
        fp2 = fp.generate(mol)

        assert fp1 == fp2
