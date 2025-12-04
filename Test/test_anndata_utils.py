"""
Unit tests for add_cassia_to_anndata function.

Tests the integration of CASSIA annotation results with Scanpy AnnData objects.
"""

import os
import sys
import tempfile
import unittest
import warnings

import numpy as np
import pandas as pd

# Add CASSIA to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'CASSIA_python'))

try:
    import anndata
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False


class TestAnndataUtils(unittest.TestCase):
    """Tests for add_cassia_to_anndata function."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        if not ANNDATA_AVAILABLE:
            return

        # Create sample CASSIA results DataFrame
        cls.cassia_df = pd.DataFrame({
            'Cluster ID': ['0', '1', '2', '3'],
            'Predicted General Cell Type': ['T cells', 'B cells', 'Macrophages', 'NK cells'],
            'Predicted Detailed Cell Type': [
                'CD4+ T cells, CD8+ T cells, Regulatory T cells',
                'Memory B cells, Naive B cells',
                'M1 Macrophages',
                'CD56bright NK cells, CD56dim NK cells'
            ],
            'Possible Mixed Cell Types': ['', 'B cells / Plasma cells', '', ''],
            'Score': [85, 90, 78, 92],
            'Merged_Grouping_1': ['Lymphocytes', 'Lymphocytes', 'Myeloid', 'Lymphocytes'],
            'Merged_Grouping_2': ['T cells', 'B cells', 'Macrophages', 'NK cells'],
            'Merged_Grouping_3': ['CD4+ T cells', 'Memory B cells', 'M1 Macrophages', 'CD56bright NK']
        })

        # Create sample AnnData object
        n_cells = 100
        n_genes = 50

        # Random expression matrix
        np.random.seed(42)
        X = np.random.rand(n_cells, n_genes)

        # Assign cells to clusters
        cluster_assignments = np.random.choice(['0', '1', '2', '3'], size=n_cells)

        cls.adata = anndata.AnnData(X=X)
        cls.adata.obs['leiden'] = pd.Categorical(cluster_assignments)
        cls.adata.obs['louvain'] = pd.Categorical(cluster_assignments)

    def setUp(self):
        """Reset adata before each test."""
        if not ANNDATA_AVAILABLE:
            self.skipTest("anndata not installed")

        # Create fresh copy for each test
        self.test_adata = self.adata.copy()

    def test_basic_functionality(self):
        """Test basic add_cassia_to_anndata with DataFrame input."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        # Check that columns were added
        self.assertIn('CASSIA_general_celltype', self.test_adata.obs.columns)
        self.assertIn('CASSIA_sub_celltype', self.test_adata.obs.columns)
        self.assertIn('CASSIA_score', self.test_adata.obs.columns)
        self.assertIn('CASSIA_combined_celltype', self.test_adata.obs.columns)

        # Check that uns was populated
        self.assertIn('CASSIA', self.test_adata.uns)

    def test_file_path_input(self):
        """Test loading from CSV file path."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        # Save to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            self.cassia_df.to_csv(f.name, index=False)
            temp_path = f.name

        try:
            add_cassia_to_anndata(self.test_adata, temp_path)
            self.assertIn('CASSIA_general_celltype', self.test_adata.obs.columns)
        finally:
            os.unlink(temp_path)

    def test_auto_detect_leiden(self):
        """Test auto-detection of leiden cluster column."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        # Remove louvain to ensure leiden is detected
        del self.test_adata.obs['louvain']

        add_cassia_to_anndata(self.test_adata, self.cassia_df)
        self.assertIn('CASSIA_general_celltype', self.test_adata.obs.columns)

    def test_auto_detect_louvain(self):
        """Test auto-detection of louvain when leiden not present."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        # Remove leiden to ensure louvain is detected
        del self.test_adata.obs['leiden']

        add_cassia_to_anndata(self.test_adata, self.cassia_df)
        self.assertIn('CASSIA_general_celltype', self.test_adata.obs.columns)

    def test_explicit_cluster_col(self):
        """Test explicit cluster column specification."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        self.test_adata.obs['my_clusters'] = self.test_adata.obs['leiden']

        add_cassia_to_anndata(
            self.test_adata,
            self.cassia_df,
            cluster_col='my_clusters'
        )
        self.assertIn('CASSIA_general_celltype', self.test_adata.obs.columns)

    def test_sub_celltype_splitting(self):
        """Test that sub-celltypes are properly split."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        # Check split columns exist
        self.assertIn('CASSIA_sub_celltype_all', self.test_adata.obs.columns)
        self.assertIn('CASSIA_sub_celltype_1', self.test_adata.obs.columns)
        self.assertIn('CASSIA_sub_celltype_2', self.test_adata.obs.columns)
        self.assertIn('CASSIA_sub_celltype_3', self.test_adata.obs.columns)

        # Check that sub_celltype equals sub_celltype_1
        pd.testing.assert_series_equal(
            self.test_adata.obs['CASSIA_sub_celltype'],
            self.test_adata.obs['CASSIA_sub_celltype_1'],
            check_names=False
        )

    def test_merged_groupings(self):
        """Test that merged grouping columns are added."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        self.assertIn('CASSIA_merged_grouping_1', self.test_adata.obs.columns)
        self.assertIn('CASSIA_merged_grouping_2', self.test_adata.obs.columns)
        self.assertIn('CASSIA_merged_grouping_3', self.test_adata.obs.columns)

    def test_columns_to_include_1(self):
        """Test columns_to_include=1 (merged groupings only)."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(
            self.test_adata,
            self.cassia_df,
            columns_to_include=1
        )

        # Should have merged groupings
        self.assertIn('CASSIA_merged_grouping_1', self.test_adata.obs.columns)

        # Should NOT have general_celltype when columns_to_include=1
        self.assertNotIn('CASSIA_general_celltype', self.test_adata.obs.columns)

    def test_combined_celltype_format(self):
        """Test combined celltype column format."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        # Check format contains "::"
        combined = self.test_adata.obs['CASSIA_combined_celltype']
        has_separator = combined.str.contains(' :: ', na=False)
        self.assertTrue(has_separator.any())

    def test_inplace_false(self):
        """Test inplace=False returns copy."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        result = add_cassia_to_anndata(
            self.test_adata,
            self.cassia_df,
            inplace=False
        )

        # Original should not be modified
        self.assertNotIn('CASSIA_general_celltype', self.test_adata.obs.columns)

        # Result should have columns
        self.assertIn('CASSIA_general_celltype', result.obs.columns)

    def test_replace_existing(self):
        """Test replace_existing parameter."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        # First call
        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        # Second call without replace should raise
        with self.assertRaises(ValueError):
            add_cassia_to_anndata(self.test_adata, self.cassia_df)

        # With replace=True should succeed
        add_cassia_to_anndata(
            self.test_adata,
            self.cassia_df,
            replace_existing=True
        )

    def test_custom_prefix(self):
        """Test custom prefix parameter."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(
            self.test_adata,
            self.cassia_df,
            prefix='MY_'
        )

        self.assertIn('MY_general_celltype', self.test_adata.obs.columns)
        self.assertNotIn('CASSIA_general_celltype', self.test_adata.obs.columns)

    def test_uns_storage(self):
        """Test cluster-level summary in adata.uns."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        add_cassia_to_anndata(self.test_adata, self.cassia_df)

        uns_df = self.test_adata.uns['CASSIA']
        self.assertIsInstance(uns_df, pd.DataFrame)
        self.assertEqual(len(uns_df), 4)  # 4 clusters

    def test_missing_cluster_col_error(self):
        """Test error when cluster column not found."""
        from CASSIA.core.anndata_utils import add_cassia_to_anndata

        # Remove both cluster columns
        del self.test_adata.obs['leiden']
        del self.test_adata.obs['louvain']

        with self.assertRaises(ValueError):
            add_cassia_to_anndata(self.test_adata, self.cassia_df)


class TestFuzzyMatching(unittest.TestCase):
    """Tests for cluster fuzzy matching."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        if not ANNDATA_AVAILABLE:
            return

    def setUp(self):
        """Set up for each test."""
        if not ANNDATA_AVAILABLE:
            self.skipTest("anndata not installed")

    def test_fuzzy_match_cluster_prefix(self):
        """Test fuzzy matching with 'cluster_' prefix."""
        from CASSIA.core.anndata_utils import _fuzzy_match_clusters

        adata_clusters = ['0', '1', '2']
        cassia_clusters = ['cluster_0', 'cluster_1', 'cluster_2']

        mapping = _fuzzy_match_clusters(adata_clusters, cassia_clusters)

        self.assertEqual(mapping['0'], 'cluster_0')
        self.assertEqual(mapping['1'], 'cluster_1')
        self.assertEqual(mapping['2'], 'cluster_2')

    def test_fuzzy_match_case_insensitive(self):
        """Test case-insensitive matching."""
        from CASSIA.core.anndata_utils import _fuzzy_match_clusters

        adata_clusters = ['T_cells', 'B_cells']
        cassia_clusters = ['T Cells', 'B Cells']

        mapping = _fuzzy_match_clusters(adata_clusters, cassia_clusters)

        # Normalized comparison should match
        self.assertIsNotNone(mapping['T_cells'])
        self.assertIsNotNone(mapping['B_cells'])

    def test_fuzzy_match_exact_priority(self):
        """Test that exact matches take priority."""
        from CASSIA.core.anndata_utils import _fuzzy_match_clusters

        adata_clusters = ['0', '1']
        cassia_clusters = ['0', '1', 'cluster_0', 'cluster_1']

        mapping = _fuzzy_match_clusters(adata_clusters, cassia_clusters)

        # Should prefer exact match
        self.assertEqual(mapping['0'], '0')
        self.assertEqual(mapping['1'], '1')


class TestHelperFunctions(unittest.TestCase):
    """Tests for helper functions."""

    def test_normalize_cluster_name(self):
        """Test cluster name normalization."""
        from CASSIA.core.anndata_utils import _normalize_cluster_name

        # Underscores are converted to spaces for equivalence
        self.assertEqual(_normalize_cluster_name('Cluster_0'), 'cluster 0')
        self.assertEqual(_normalize_cluster_name('  T Cells  '), 't cells')
        self.assertEqual(_normalize_cluster_name('T_cells'), 't cells')  # Underscore -> space
        # Punctuation like + is removed
        self.assertEqual(_normalize_cluster_name('CD4+'), 'cd4')

    def test_parse_sub_celltypes(self):
        """Test sub-celltype parsing."""
        from CASSIA.core.anndata_utils import _parse_sub_celltypes

        result = _parse_sub_celltypes('A, B, C')
        self.assertEqual(result, ['A', 'B', 'C'])

        result = _parse_sub_celltypes('A, B')
        self.assertEqual(result, ['A', 'B', None])

        result = _parse_sub_celltypes('')
        self.assertEqual(result, [None, None, None])

    def test_find_column(self):
        """Test column finding."""
        from CASSIA.core.anndata_utils import _find_column

        df = pd.DataFrame({'col_a': [1], 'col_b': [2]})

        self.assertEqual(_find_column(df, ['col_a', 'col_b']), 'col_a')
        self.assertEqual(_find_column(df, ['col_c', 'col_b']), 'col_b')
        self.assertIsNone(_find_column(df, ['col_c', 'col_d']))


class TestEnhanceScanpyMarkers(unittest.TestCase):
    """Tests for enhance_scanpy_markers function."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        if not ANNDATA_AVAILABLE:
            return

        # Create sample AnnData object with expression data
        n_cells = 100
        n_genes = 20

        np.random.seed(42)

        # Create gene names
        gene_names = [f'Gene{i}' for i in range(n_genes)]

        # Create expression matrix (some genes highly expressed in specific clusters)
        X = np.random.rand(n_cells, n_genes) * 0.1  # Low baseline

        # Assign cells to clusters
        cluster_assignments = np.array(['0'] * 25 + ['1'] * 25 + ['2'] * 25 + ['3'] * 25)

        # Make certain genes highly expressed in specific clusters
        # Cluster 0: Gene0, Gene1 highly expressed
        X[:25, 0] = np.random.rand(25) * 5 + 2
        X[:25, 1] = np.random.rand(25) * 4 + 1

        # Cluster 1: Gene2, Gene3 highly expressed
        X[25:50, 2] = np.random.rand(25) * 5 + 2
        X[25:50, 3] = np.random.rand(25) * 4 + 1

        cls.adata = anndata.AnnData(X=X)
        cls.adata.var_names = gene_names
        cls.adata.obs['leiden'] = pd.Categorical(cluster_assignments)

        # Create mock rank_genes_groups results (structured array format)
        n_ranked = 10
        cls.adata.uns['rank_genes_groups'] = {
            'params': {'groupby': 'leiden', 'method': 'wilcoxon'},
            'names': np.array(
                [(f'Gene{i}', f'Gene{i+2}', f'Gene{i+4}', f'Gene{i+6}') for i in range(n_ranked)],
                dtype=[('0', 'U10'), ('1', 'U10'), ('2', 'U10'), ('3', 'U10')]
            ),
            'scores': np.array(
                [(10-i, 10-i, 10-i, 10-i) for i in range(n_ranked)],
                dtype=[('0', 'f8'), ('1', 'f8'), ('2', 'f8'), ('3', 'f8')]
            ),
            'logfoldchanges': np.array(
                [(2.0-i*0.1, 2.0-i*0.1, 2.0-i*0.1, 2.0-i*0.1) for i in range(n_ranked)],
                dtype=[('0', 'f8'), ('1', 'f8'), ('2', 'f8'), ('3', 'f8')]
            ),
            'pvals_adj': np.array(
                [(0.001*i+0.0001, 0.001*i+0.0001, 0.001*i+0.0001, 0.001*i+0.0001) for i in range(n_ranked)],
                dtype=[('0', 'f8'), ('1', 'f8'), ('2', 'f8'), ('3', 'f8')]
            )
        }

    def setUp(self):
        """Reset adata before each test."""
        if not ANNDATA_AVAILABLE:
            self.skipTest("anndata not installed")

        self.test_adata = self.adata.copy()

    def test_basic_functionality(self):
        """Test basic enhance_scanpy_markers functionality."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=5)

        # Check required columns exist
        self.assertIn('cluster', df.columns)
        self.assertIn('gene', df.columns)
        self.assertIn('pct.1', df.columns)
        self.assertIn('pct.2', df.columns)

        # Check we have results for all clusters
        self.assertEqual(len(df['cluster'].unique()), 4)

    def test_pct_values_range(self):
        """Test that pct values are in valid range [0, 1]."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=5)

        # All pct values should be between 0 and 1
        self.assertTrue((df['pct.1'] >= 0).all())
        self.assertTrue((df['pct.1'] <= 1).all())
        self.assertTrue((df['pct.2'] >= 0).all())
        self.assertTrue((df['pct.2'] <= 1).all())

    def test_n_genes_parameter(self):
        """Test n_genes parameter limits output."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=3)

        # Should have at most 3 genes per cluster (4 clusters * 3 genes = 12)
        self.assertLessEqual(len(df), 12)

    def test_include_stats_true(self):
        """Test include_stats=True includes additional columns."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=5, include_stats=True)

        # Should have additional stat columns
        self.assertIn('avg_log2FC', df.columns)
        self.assertIn('p_val_adj', df.columns)
        self.assertIn('scores', df.columns)

    def test_include_stats_false(self):
        """Test include_stats=False excludes additional columns."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=5, include_stats=False)

        # Should NOT have additional stat columns
        self.assertNotIn('avg_log2FC', df.columns)
        self.assertNotIn('p_val_adj', df.columns)
        self.assertNotIn('scores', df.columns)

    def test_missing_rank_genes_groups_error(self):
        """Test error when rank_genes_groups not found."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        # Remove rank_genes_groups
        del self.test_adata.uns['rank_genes_groups']

        with self.assertRaises(ValueError):
            enhance_scanpy_markers(self.test_adata)

    def test_auto_detect_cluster_col(self):
        """Test auto-detection of cluster column."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        # Should auto-detect 'leiden'
        df = enhance_scanpy_markers(self.test_adata, n_genes=3)
        self.assertIsNotNone(df)

    def test_explicit_cluster_col(self):
        """Test explicit cluster column specification."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        self.test_adata.obs['my_clusters'] = self.test_adata.obs['leiden']

        df = enhance_scanpy_markers(
            self.test_adata,
            cluster_col='my_clusters',
            n_genes=3
        )
        self.assertIsNotNone(df)

    def test_output_format_compatible(self):
        """Test output format is compatible with annotation boost agent."""
        from CASSIA.core.anndata_utils import enhance_scanpy_markers

        df = enhance_scanpy_markers(self.test_adata, n_genes=5)

        # Column order should be: cluster, gene, pct.1, pct.2, [stats...]
        expected_cols = ['cluster', 'gene', 'pct.1', 'pct.2']
        for i, col in enumerate(expected_cols):
            self.assertEqual(df.columns[i], col)


class TestCalculatePctExpressing(unittest.TestCase):
    """Tests for _calculate_pct_expressing helper function."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        if not ANNDATA_AVAILABLE:
            return

    def setUp(self):
        """Set up for each test."""
        if not ANNDATA_AVAILABLE:
            self.skipTest("anndata not installed")

    def test_pct_calculation_basic(self):
        """Test basic percentage calculation."""
        from CASSIA.core.anndata_utils import _calculate_pct_expressing

        # Create simple AnnData
        X = np.array([
            [1, 0],  # Cluster 0, cell 1
            [1, 0],  # Cluster 0, cell 2
            [0, 1],  # Cluster 1, cell 1
            [0, 1],  # Cluster 1, cell 2
        ])

        adata = anndata.AnnData(X=X)
        adata.var_names = ['GeneA', 'GeneB']
        adata.obs['cluster'] = pd.Categorical(['0', '0', '1', '1'])

        # GeneA: 100% in cluster 0, 0% in cluster 1
        pct1, pct2 = _calculate_pct_expressing(adata, 'GeneA', 'cluster', '0')
        self.assertEqual(pct1, 1.0)  # 100% in cluster
        self.assertEqual(pct2, 0.0)  # 0% outside cluster

    def test_pct_with_threshold(self):
        """Test percentage calculation with expression threshold."""
        from CASSIA.core.anndata_utils import _calculate_pct_expressing

        X = np.array([
            [0.5, 0],  # Cluster 0
            [1.5, 0],  # Cluster 0
            [0.5, 1],  # Cluster 1
            [0.5, 1],  # Cluster 1
        ])

        adata = anndata.AnnData(X=X)
        adata.var_names = ['GeneA', 'GeneB']
        adata.obs['cluster'] = pd.Categorical(['0', '0', '1', '1'])

        # With threshold=1.0, only one cell in cluster 0 expresses GeneA
        pct1, pct2 = _calculate_pct_expressing(adata, 'GeneA', 'cluster', '0', min_expression=1.0)
        self.assertEqual(pct1, 0.5)  # 50% in cluster (1 of 2)
        self.assertEqual(pct2, 0.0)  # 0% outside

    def test_gene_not_found(self):
        """Test handling of gene not in adata."""
        from CASSIA.core.anndata_utils import _calculate_pct_expressing

        X = np.array([[1, 0], [0, 1]])
        adata = anndata.AnnData(X=X)
        adata.var_names = ['GeneA', 'GeneB']
        adata.obs['cluster'] = pd.Categorical(['0', '1'])

        pct1, pct2 = _calculate_pct_expressing(adata, 'NonExistentGene', 'cluster', '0')
        self.assertTrue(np.isnan(pct1))
        self.assertTrue(np.isnan(pct2))


if __name__ == '__main__':
    unittest.main()
