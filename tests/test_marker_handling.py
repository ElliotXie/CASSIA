"""
Unit tests for marker data handling functions in CASSIA.
"""

import unittest
import pandas as pd
import os
import tempfile
from CASSIA import split_markers, get_top_markers

class TestMarkerHandling(unittest.TestCase):
    """Test marker data handling functions."""
    
    def setUp(self):
        """Set up test data before each test."""
        # Create a sample Seurat-format marker dataframe
        self.seurat_data = pd.DataFrame({
            'cluster': [1, 1, 1, 2, 2, 2, 3, 3, 3],
            'gene': ['CD4', 'CD8A', 'IL7R', 'CD19', 'MS4A1', 'CD79A', 'CD14', 'LYZ', 'CSF1R'],
            'p_val_adj': [0.001, 0.002, 0.003, 0.001, 0.002, 0.003, 0.001, 0.002, 0.003],
            'avg_log2FC': [2.5, 2.0, 1.8, 3.0, 2.8, 2.6, 3.2, 3.0, 2.8],
            'pct.1': [0.8, 0.7, 0.6, 0.9, 0.85, 0.8, 0.95, 0.9, 0.85],
            'pct.2': [0.1, 0.15, 0.2, 0.1, 0.12, 0.15, 0.05, 0.1, 0.15]
        })
        
        # Create sample marker strings
        self.marker_strings = {
            'comma_separated': 'CD4,CD8A,IL7R',
            'comma_space_separated': 'CD4, CD8A, IL7R',
            'space_separated': 'CD4 CD8A IL7R',
            'mixed_separators': 'CD4, CD8A IL7R',
            'single_marker': 'CD4'
        }
    
    def test_split_markers_comma_separated(self):
        """Test splitting comma-separated markers."""
        markers = split_markers(self.marker_strings['comma_separated'])
        self.assertEqual(markers, ['CD4', 'CD8A', 'IL7R'])
    
    def test_split_markers_comma_space_separated(self):
        """Test splitting comma-and-space-separated markers."""
        markers = split_markers(self.marker_strings['comma_space_separated'])
        self.assertEqual(markers, ['CD4', 'CD8A', 'IL7R'])
    
    def test_split_markers_space_separated(self):
        """Test splitting space-separated markers."""
        markers = split_markers(self.marker_strings['space_separated'])
        self.assertEqual(markers, ['CD4', 'CD8A', 'IL7R'])
    
    def test_split_markers_mixed_separators(self):
        """Test splitting markers with mixed separators."""
        markers = split_markers(self.marker_strings['mixed_separators'])
        self.assertEqual(markers, ['CD4', 'CD8A', 'IL7R'])
    
    def test_split_markers_single_marker(self):
        """Test splitting a single marker."""
        markers = split_markers(self.marker_strings['single_marker'])
        self.assertEqual(markers, ['CD4'])
    
    def test_get_top_markers_seurat(self):
        """Test getting top markers from Seurat format data."""
        result = get_top_markers(self.seurat_data, n_genes=2)
        
        # Check that the result has the right columns
        self.assertIn('cluster', result.columns)
        self.assertIn('markers', result.columns)
        
        # Check that all clusters are represented
        self.assertSetEqual(set(result['cluster']), {1, 2, 3})
        
        # Check the marker content for a specific cluster
        cluster_1_markers = result[result['cluster'] == 1]['markers'].iloc[0]
        self.assertIn('CD4', cluster_1_markers)
        self.assertIn('CD8A', cluster_1_markers)
        
        # Check the marker content for another cluster
        cluster_2_markers = result[result['cluster'] == 2]['markers'].iloc[0]
        self.assertIn('CD19', cluster_2_markers)
        self.assertIn('MS4A1', cluster_2_markers)
    
    def test_get_top_markers_n_genes(self):
        """Test getting different numbers of top markers."""
        # Get just 1 top marker per cluster
        result_1 = get_top_markers(self.seurat_data, n_genes=1)
        
        # Check that each cluster's markers string has just 1 gene
        for _, row in result_1.iterrows():
            self.assertEqual(len(row['markers'].split(',')), 1)
        
        # Get all 3 markers per cluster
        result_3 = get_top_markers(self.seurat_data, n_genes=3)
        
        # Check that each cluster's markers string has 3 genes
        for _, row in result_3.iterrows():
            self.assertEqual(len(row['markers'].split(',')), 3)
    
    def test_get_top_markers_filter(self):
        """Test filtering of markers based on statistical criteria."""
        # Create data with some markers that should be filtered out
        data = self.seurat_data.copy()
        # Add a row with low p-value
        data.loc[len(data)] = [1, 'GENE1', 0.06, 2.5, 0.8, 0.1]  # p_val_adj > 0.05
        # Add a row with low log fold change
        data.loc[len(data)] = [1, 'GENE2', 0.01, 0.2, 0.8, 0.1]  # avg_log2FC < 0.25
        # Add a row with low percentage
        data.loc[len(data)] = [1, 'GENE3', 0.01, 2.5, 0.05, 0.05]  # pct.1 and pct.2 < 0.1
        
        result = get_top_markers(data, n_genes=10)
        
        # Get the markers for cluster 1
        cluster_1_markers = result[result['cluster'] == 1]['markers'].iloc[0].split(',')
        
        # Check that filtered genes are not in the result
        self.assertNotIn('GENE1', cluster_1_markers)
        self.assertNotIn('GENE2', cluster_1_markers)
        self.assertNotIn('GENE3', cluster_1_markers)
        
        # Check that the original genes are still there
        self.assertIn('CD4', cluster_1_markers)
        self.assertIn('CD8A', cluster_1_markers)
        self.assertIn('IL7R', cluster_1_markers)

    def test_marker_format_detection(self):
        """Test auto-detection of marker format (Seurat vs Scanpy)."""
        # This is a simplified test since we can't easily create a Scanpy format dataset
        # Make sure Seurat format is properly detected
        result = get_top_markers(self.seurat_data, format_type=None)  # Auto-detect
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn('cluster', result.columns)
        self.assertIn('markers', result.columns)


if __name__ == "__main__":
    unittest.main() 