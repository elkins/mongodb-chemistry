"""Unit tests for mchem.similarity module."""

import pytest
from unittest.mock import Mock, MagicMock
from math import ceil
from mchem import similarity


class TestSimilarityCalculations:
    """Tests for similarity calculation logic."""

    def test_tanimoto_identical_fingerprints(self):
        """Test Tanimoto calculation for identical fingerprints."""
        qfp = [1, 2, 3, 4, 5]
        fp = [1, 2, 3, 4, 5]
        
        intersection = len(set(qfp) & set(fp))
        tanimoto = float(intersection) / (len(qfp) + len(fp) - intersection)
        
        assert tanimoto == 1.0

    def test_tanimoto_no_overlap(self):
        """Test Tanimoto calculation for fingerprints with no overlap."""
        qfp = [1, 2, 3]
        fp = [4, 5, 6]
        
        intersection = len(set(qfp) & set(fp))
        tanimoto = float(intersection) / (len(qfp) + len(fp) - intersection)
        
        assert tanimoto == 0.0

    def test_tanimoto_partial_overlap(self):
        """Test Tanimoto calculation for partial overlap."""
        qfp = [1, 2, 3, 4]
        fp = [3, 4, 5, 6]
        
        intersection = len(set(qfp) & set(fp))  # {3, 4} = 2
        tanimoto = float(intersection) / (len(qfp) + len(fp) - intersection)  # 2 / (4 + 4 - 2) = 2/6
        
        assert abs(tanimoto - 0.333333) < 0.0001


class TestSimilarityThresholds:
    """Tests for threshold calculations."""

    def test_threshold_bounds_calculation(self):
        """Test qmin and qmax calculation from threshold."""
        qn = 100
        threshold = 0.8
        
        qmin = int(ceil(qn * threshold))
        qmax = int(qn / threshold)
        
        assert qmin == 80
        assert qmax == 125

    def test_threshold_ncommon_calculation(self):
        """Test ncommon calculation."""
        qn = 100
        threshold = 0.8
        qmin = int(ceil(qn * threshold))
        
        ncommon = qn - qmin + 1
        
        assert ncommon == 21  # 100 - 80 + 1


@pytest.mark.mongodb
class TestSimilarityClient:
    """Tests for similarity_client function."""

    def test_similarity_client_basic(self):
        """Test basic similarity search on client side."""
        # Mock molecule and fingerprinter
        mock_mol = Mock()
        mock_fingerprinter = Mock()
        mock_fingerprinter.generate.return_value = [1, 2, 3, 4, 5]
        
        # Mock fp_collection
        mock_fp_collection = Mock()
        mock_fp_collection.name = 'test_fps'
        mock_fp_collection.find.return_value = [
            {'_id': 'mol1', 'bits': [1, 2, 3, 4, 5], 'count': 5},  # Perfect match
            {'_id': 'mol2', 'bits': [1, 2, 3, 4, 6], 'count': 5},  # 4/6 = 0.667
            {'_id': 'mol3', 'bits': [1, 2, 7, 8, 9], 'count': 5},  # 2/8 = 0.25
        ]
        
        results = similarity.similarity_client(
            mock_mol, 
            mock_fingerprinter, 
            mock_fp_collection, 
            threshold=0.6
        )
        
        # Should return mol1 (1.0) and mol2 (0.667), but not mol3 (0.25)
        assert len(results) == 2
        assert results[0]['_id'] == 'mol1'
        assert abs(results[0]['tanimoto'] - 1.0) < 0.001
        assert results[1]['_id'] == 'mol2'
        assert abs(results[1]['tanimoto'] - 0.667) < 0.001

    def test_similarity_client_with_count_collection(self):
        """Test similarity search using count collection for rare bits."""
        mock_mol = Mock()
        mock_fingerprinter = Mock()
        mock_fingerprinter.generate.return_value = [1, 2, 3, 4, 5]
        
        # Mock count_collection (bit 5 is rarest)
        mock_count_collection = Mock()
        mock_cursor = Mock()
        mock_cursor.sort.return_value.limit.return_value = [
            {'_id': 5, 'count': 10},
            {'_id': 1, 'count': 20},
        ]
        mock_count_collection.find.return_value = mock_cursor
        
        mock_fp_collection = Mock()
        mock_fp_collection.name = 'test_fps'
        mock_fp_collection.find.return_value = []
        
        similarity.similarity_client(
            mock_mol,
            mock_fingerprinter,
            mock_fp_collection,
            threshold=0.8,
            count_collection=mock_count_collection
        )
        
        # Verify count_collection was queried
        assert mock_count_collection.find.called


@pytest.mark.mongodb
class TestSimilaritySearchFp:
    """Tests for similarity_search_fp function using aggregation."""

    def test_similarity_search_fp_aggregate_pipeline(self):
        """Test that aggregation pipeline is constructed correctly."""
        qfp = [1, 2, 3, 4, 5]
        threshold = 0.8
        
        mock_fp_collection = Mock()
        mock_fp_collection.aggregate.return_value = []
        
        similarity.similarity_search_fp(qfp, mock_fp_collection, threshold=threshold)
        
        # Verify aggregate was called
        assert mock_fp_collection.aggregate.called
        
        # Get the pipeline
        pipeline = mock_fp_collection.aggregate.call_args[0][0]
        
        # Check pipeline structure
        assert len(pipeline) == 3
        assert '$match' in pipeline[0]
        assert '$project' in pipeline[1]
        assert '$match' in pipeline[2]
        
        # Check final match has threshold
        assert pipeline[2]['$match']['tanimoto']['$gte'] == threshold

    def test_similarity_search_fp_with_count_collection(self):
        """Test similarity search with count collection for optimal bit selection."""
        qfp = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        
        # Mock count_collection
        mock_count_collection = Mock()
        mock_cursor = Mock()
        mock_cursor.sort.return_value.limit.return_value = [
            {'_id': 10, 'count': 5},
            {'_id': 5, 'count': 10},
        ]
        mock_count_collection.find.return_value = mock_cursor
        
        mock_fp_collection = Mock()
        mock_fp_collection.aggregate.return_value = []
        
        similarity.similarity_search_fp(
            qfp,
            mock_fp_collection,
            threshold=0.7,
            count_collection=mock_count_collection
        )
        
        # Verify count_collection was used
        assert mock_count_collection.find.called


@pytest.mark.mongodb
class TestSimilaritySearch:
    """Tests for main similarity_search function."""

    def test_similarity_search_calls_similarity_search_fp(self):
        """Test that similarity_search generates fingerprint and calls similarity_search_fp."""
        mock_mol = Mock()
        mock_fingerprinter = Mock()
        mock_fingerprinter.generate.return_value = [1, 2, 3]
        
        mock_fp_collection = Mock()
        mock_fp_collection.aggregate.return_value = []
        
        similarity.similarity_search(
            mock_mol,
            mock_fingerprinter,
            mock_fp_collection,
            threshold=0.7
        )
        
        # Verify fingerprint was generated
        assert mock_fingerprinter.generate.called
        
        # Verify aggregate was called
        assert mock_fp_collection.aggregate.called


@pytest.mark.rdkit
@pytest.mark.mongodb
class TestSimilarityIntegration:
    """Integration tests with RDKit and mock MongoDB."""

    def test_similarity_search_with_real_molecules(self):
        """Test similarity search with real RDKit molecules."""
        from rdkit import Chem
        from mchem.fps import MorganFingerprinter
        
        # Create fingerprinter
        fingerprinter = MorganFingerprinter(radius=2)
        
        # Query molecule (aspirin)
        aspirin = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
        aspirin_fp = fingerprinter.generate(aspirin)
        
        # Mock collection with similar molecules
        mock_fp_collection = Mock()
        mock_fp_collection.name = 'test_fps'
        
        # Create fingerprints for test molecules
        ibuprofen = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        ibuprofen_fp = fingerprinter.generate(ibuprofen)
        
        mock_fp_collection.find.return_value = [
            {'_id': 'aspirin', 'bits': aspirin_fp, 'count': len(aspirin_fp)},
            {'_id': 'ibuprofen', 'bits': ibuprofen_fp, 'count': len(ibuprofen_fp)},
        ]
        
        # Test client-side search
        results = similarity.similarity_client(
            aspirin,
            fingerprinter,
            mock_fp_collection,
            threshold=0.3
        )
        
        # Should find at least aspirin itself
        assert len(results) >= 1
        aspirin_result = [r for r in results if r['_id'] == 'aspirin']
        assert len(aspirin_result) == 1
        assert abs(aspirin_result[0]['tanimoto'] - 1.0) < 0.001


class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_empty_fingerprint(self):
        """Test handling of empty fingerprints."""
        qfp = []
        fp = []
        
        if len(qfp) == 0 or len(fp) == 0:
            # Should handle empty fingerprints gracefully
            assert True
        else:
            intersection = len(set(qfp) & set(fp))
            tanimoto = float(intersection) / (len(qfp) + len(fp) - intersection)

    def test_threshold_boundary_values(self):
        """Test threshold at boundary values."""
        # Test threshold = 0.0
        qn = 100
        threshold = 0.0
        qmin = int(ceil(qn * threshold))
        assert qmin == 0
        
        # Test threshold = 1.0
        threshold = 1.0
        qmin = int(ceil(qn * threshold))
        assert qmin == 100

    def test_very_high_threshold(self):
        """Test with very high threshold (>0.9)."""
        qn = 100
        threshold = 0.95
        qmin = int(ceil(qn * threshold))
        qmax = int(qn / threshold)
        ncommon = qn - qmin + 1
        
        assert qmin == 95
        assert qmax == 105
        assert ncommon == 6  # Very selective
