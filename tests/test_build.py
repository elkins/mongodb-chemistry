"""Unit tests for mchem.build module."""

from unittest.mock import Mock, mock_open, patch

import pytest

# Only import if RDKit is available
pytest.importorskip("rdkit")

from mchem import build


class TestLoadSdf:
    """Tests for load_sdf function."""

    @patch("mchem.build.Chem.ForwardSDMolSupplier")
    def test_load_sdf_success(self, mock_supplier):
        """Test successful SDF loading."""
        # Mock RDKit molecule
        mock_mol = Mock()
        mock_mol.GetProp.return_value = "CHEMBL123"
        mock_mol.ToBinary.return_value = b"mock_binary_data"

        # Mock SMILES generation
        with patch("mchem.build.Chem.MolToSmiles", return_value="CC(=O)O"):
            # Mock supplier to return one molecule
            mock_supplier.return_value = [mock_mol]

            # Mock collection
            mock_collection = Mock()
            mock_sdf = Mock()

            build.load_sdf(mock_sdf, mock_collection, idfield="CHEMBL_ID")

            # Verify insert_one was called
            assert mock_collection.insert_one.called
            call_args = mock_collection.insert_one.call_args[0][0]
            assert call_args["smiles"] == "CC(=O)O"
            assert call_args["rdmol"]

    @patch("mchem.build.Chem.ForwardSDMolSupplier")
    def test_load_sdf_skip_duplicate(self, mock_supplier):
        """Test that duplicate molecules are skipped."""
        from pymongo.errors import DuplicateKeyError

        mock_mol = Mock()
        mock_mol.GetProp.return_value = "CHEMBL123"
        mock_mol.ToBinary.return_value = b"mock_data"

        with patch("mchem.build.Chem.MolToSmiles", return_value="CC(=O)O"):
            mock_supplier.return_value = [mock_mol]

            mock_collection = Mock()
            mock_collection.insert_one.side_effect = DuplicateKeyError("Duplicate")

            mock_sdf = Mock()

            # Should not raise exception
            build.load_sdf(mock_sdf, mock_collection, idfield="CHEMBL_ID")

    @patch("mchem.build.Chem.ForwardSDMolSupplier")
    def test_load_sdf_skip_invalid_molecules(self, mock_supplier):
        """Test that invalid molecules (None) are skipped."""
        # Return None for invalid molecule
        mock_supplier.return_value = [None, None]

        mock_collection = Mock()
        mock_sdf = Mock()

        build.load_sdf(mock_sdf, mock_collection, idfield=None)

        # insert_one should not be called for None molecules
        assert not mock_collection.insert_one.called

    @patch("mchem.build.Chem.ForwardSDMolSupplier")
    def test_load_sdf_without_idfield(self, mock_supplier):
        """Test loading without explicit ID field."""
        mock_mol = Mock()
        mock_mol.ToBinary.return_value = b"mock_data"

        with patch("mchem.build.Chem.MolToSmiles", return_value="CC(=O)O"):
            mock_supplier.return_value = [mock_mol]

            mock_collection = Mock()
            mock_sdf = Mock()

            build.load_sdf(mock_sdf, mock_collection, idfield=None)

            # Verify molecule was inserted without explicit _id
            assert mock_collection.insert_one.called
            call_args = mock_collection.insert_one.call_args[0][0]
            assert "_id" not in call_args or call_args.get("_id") is None


class TestLoadSdfs:
    """Tests for load_sdfs function."""

    @patch("mchem.build.gzip.open")
    @patch("mchem.build.open")
    @patch("mchem.build.load_sdf")
    def test_load_sdfs_gzip_file(self, mock_load_sdf, mock_open_file, mock_gzip_open):
        """Test loading gzipped SDF file."""
        paths = ["molecules.sdf.gz"]
        mock_collection = Mock()

        build.load_sdfs(paths, mock_collection, idfield="CHEMBL_ID")

        # Verify gzip.open was called
        assert mock_gzip_open.called
        # Verify load_sdf was called
        assert mock_load_sdf.called

    @patch("mchem.build.open")
    @patch("mchem.build.load_sdf")
    def test_load_sdfs_regular_file(self, mock_load_sdf, mock_open_file):
        """Test loading regular SDF file."""
        paths = ["molecules.sdf"]
        mock_collection = Mock()

        build.load_sdfs(paths, mock_collection, idfield="CHEMBL_ID")

        # Verify open was called
        assert mock_open_file.called
        # Verify load_sdf was called
        assert mock_load_sdf.called

    @patch("mchem.build.open")
    @patch("mchem.build.load_sdf")
    def test_load_sdfs_multiple_files(self, mock_load_sdf, mock_open_file):
        """Test loading multiple SDF files."""
        paths = ["file1.sdf", "file2.sdf", "file3.sdf"]
        mock_collection = Mock()

        build.load_sdfs(paths, mock_collection, idfield=None)

        # Verify load_sdf was called for each file
        assert mock_load_sdf.call_count == 3


class TestDownloadChembl:
    """Tests for download_chembl function."""

    @patch("mchem.build.os.path.isfile")
    @patch("mchem.build.urllib.request.urlopen")
    def test_download_chembl_new_file(self, mock_urlopen, mock_isfile):
        """Test downloading ChEMBL when file doesn't exist."""
        mock_isfile.return_value = False

        mock_response = Mock()
        mock_response.read.return_value = b"mock_data"
        mock_urlopen.return_value = mock_response

        with patch("builtins.open", mock_open()) as mock_file:
            build.download_chembl()

            # Verify download was attempted
            assert mock_urlopen.called
            assert mock_file.called

    @patch("mchem.build.os.path.isfile")
    def test_download_chembl_existing_file(self, mock_isfile):
        """Test that existing file is not re-downloaded."""
        mock_isfile.return_value = True

        with patch("mchem.build.urllib.request.urlopen") as mock_urlopen:
            build.download_chembl()

            # Verify download was not attempted
            assert not mock_urlopen.called


@pytest.mark.rdkit
class TestBuildIntegration:
    """Integration tests with actual RDKit."""

    def test_process_valid_smiles(self):
        """Test processing a valid SMILES string."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles("CC(=O)O")
        assert mol is not None

        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        assert smiles == "CC(=O)O"

        binary = mol.ToBinary()
        assert isinstance(binary, bytes)
        assert len(binary) > 0

    def test_process_invalid_smiles(self):
        """Test handling of invalid SMILES."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles("INVALID_SMILES_XXX")
        assert mol is None

    def test_round_trip_binary_conversion(self):
        """Test molecule binary serialization round-trip."""
        from rdkit import Chem

        original = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        binary = original.ToBinary()
        restored = Chem.Mol(binary)

        original_smiles = Chem.MolToSmiles(original)
        restored_smiles = Chem.MolToSmiles(restored)

        assert original_smiles == restored_smiles
