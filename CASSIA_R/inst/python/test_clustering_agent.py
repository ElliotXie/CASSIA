#!/usr/bin/env python3
"""
Test script for CASSIA Clustering Agent functions.

This script tests all the clustering agent functions to ensure they work correctly.
"""

import sys
import os
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    
    try:
        import scanpy as sc
        print("✓ scanpy imported successfully")
    except ImportError as e:
        print(f"✗ scanpy import failed: {e}")
        return False
    
    try:
        import matplotlib.pyplot as plt
        print("✓ matplotlib imported successfully")
    except ImportError as e:
        print(f"✗ matplotlib import failed: {e}")
        return False
    
    try:
        from clustering_agent import (
            generate_parameter_grid,
            test_resolution_stability,
            analyze_clustering_with_llm,
            generate_final_visualization,
            preprocess_for_clustering
        )
        print("✓ clustering_agent functions imported successfully")
    except ImportError as e:
        print(f"✗ clustering_agent import failed: {e}")
        return False
    
    try:
        from llm_utils_image import analyze_clustering_images
        print("✓ llm_utils_image imported successfully")
    except ImportError as e:
        print(f"✗ llm_utils_image import failed: {e}")
        print("  (This is expected if testing without LLM functionality)")
    
    return True


def create_test_data():
    """Create synthetic test data for clustering."""
    print("Creating synthetic test data...")
    
    import scanpy as sc
    
    # Create synthetic data
    n_obs = 1000
    n_vars = 2000
    
    # Generate random expression data
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars)).astype(float)
    
    # Add some structure (3 main clusters)
    cluster_centers = [
        (100, 150),  # Cluster 1
        (300, 200),  # Cluster 2  
        (500, 100)   # Cluster 3
    ]
    
    for i, (start, end) in enumerate(cluster_centers):
        cluster_size = n_obs // 3
        if i == 2:  # Last cluster gets remaining cells
            cluster_size = n_obs - (2 * (n_obs // 3))
        
        start_idx = i * (n_obs // 3)
        end_idx = start_idx + cluster_size
        
        # Enhance expression in specific genes for each cluster
        gene_start = start
        gene_end = min(end, n_vars)
        X[start_idx:end_idx, gene_start:gene_end] *= 3
    
    # Create AnnData object
    adata = sc.AnnData(X)
    
    # Add some metadata
    adata.obs['cell_type'] = ['unknown'] * n_obs
    adata.var['gene_names'] = [f'Gene_{i}' for i in range(n_vars)]
    adata.var_names = adata.var['gene_names']
    
    print(f"✓ Created test data: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def test_preprocessing(adata):
    """Test preprocessing function."""
    print("\nTesting preprocessing...")
    
    from clustering_agent import preprocess_for_clustering
    
    try:
        adata_processed = preprocess_for_clustering(adata, n_top_genes=1000, n_pcs=20)
        
        # Check that preprocessing was successful
        assert 'X_pca' in adata_processed.obsm, "PCA not computed"
        assert adata_processed.obsm['X_pca'].shape[1] == 20, "Wrong number of PCs"
        assert hasattr(adata_processed, 'raw'), "Raw data not saved"
        
        print("✓ Preprocessing completed successfully")
        return adata_processed
        
    except Exception as e:
        print(f"✗ Preprocessing failed: {e}")
        return None


def test_parameter_grid(adata):
    """Test parameter grid generation."""
    print("\nTesting parameter grid generation...")
    
    from clustering_agent import generate_parameter_grid
    
    try:
        # Test with small parameter set
        grid_results = generate_parameter_grid(
            adata,
            output_dir="test_grid",
            n_pcs_list=[15, 20],
            min_dist_list=[0.2, 0.3],
            n_neighbors_list=[20, 30],
            resolution=0.5,
            plot_width=6,
            plot_height=5
        )
        
        # Check results
        assert 'parameter_combinations' in grid_results, "Parameter combinations not returned"
        assert 'plot_paths' in grid_results, "Plot paths not returned"
        assert 'combined_plot_path' in grid_results, "Combined plot path not returned"
        assert 'results_summary' in grid_results, "Results summary not returned"
        
        # Check that plots were created
        assert len(grid_results['plot_paths']) == 4, "Wrong number of plots created"
        
        for plot_path in grid_results['plot_paths']:
            assert os.path.exists(plot_path), f"Plot not created: {plot_path}"
        
        assert os.path.exists(grid_results['combined_plot_path']), "Combined plot not created"
        
        print("✓ Parameter grid generation completed successfully")
        print(f"  - Generated {len(grid_results['plot_paths'])} individual plots")
        print(f"  - Created combined plot: {grid_results['combined_plot_path']}")
        
        return grid_results
        
    except Exception as e:
        print(f"✗ Parameter grid generation failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_resolution_stability(adata):
    """Test resolution stability analysis."""
    print("\nTesting resolution stability analysis...")
    
    from clustering_agent import test_resolution_stability
    
    try:
        # Test with small resolution set
        stability_results = test_resolution_stability(
            adata,
            output_dir="test_stability",
            resolution_list=[0.3, 0.5, 0.7],
            n_pcs=20,
            min_dist=0.2,
            n_neighbors=30,
            n_iterations=3,  # Reduced for testing
            plot_width=8,
            plot_height=6
        )
        
        # Check results
        assert 'stability_results' in stability_results, "Stability results not returned"
        assert 'plot_paths' in stability_results, "Plot paths not returned"
        assert 'stability_summary_path' in stability_results, "Stability summary not returned"
        
        # Check that plots were created
        assert len(stability_results['plot_paths']) == 3, "Wrong number of stability plots"
        
        for plot_path in stability_results['plot_paths']:
            assert os.path.exists(plot_path), f"Stability plot not created: {plot_path}"
        
        assert os.path.exists(stability_results['stability_summary_path']), "Stability summary not created"
        
        # Check stability data
        for result in stability_results['stability_results']:
            assert 'resolution' in result, "Resolution not in stability result"
            assert 'mean_stability' in result, "Mean stability not in result"
            assert 'std_stability' in result, "Std stability not in result"
            assert 'mean_n_clusters' in result, "Mean cluster count not in result"
        
        print("✓ Resolution stability analysis completed successfully")
        print(f"  - Tested {len(stability_results['stability_results'])} resolutions")
        print(f"  - Created stability summary: {stability_results['stability_summary_path']}")
        
        return stability_results
        
    except Exception as e:
        print(f"✗ Resolution stability analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_fallback_analysis(grid_results, stability_results):
    """Test fallback analysis (without LLM)."""
    print("\nTesting fallback analysis...")
    
    from clustering_agent import analyze_clustering_with_llm
    
    try:
        # This should use fallback analysis since LLM won't be available
        analysis_results = analyze_clustering_with_llm(
            grid_results,
            stability_results,
            model="google/gemini-2.5-flash",
            provider="openrouter",
            api_key=None  # No API key to trigger fallback
        )
        
        # Check results
        assert 'recommended_parameters' in analysis_results, "Recommended parameters not returned"
        
        recommended = analysis_results['recommended_parameters']
        assert 'n_pcs' in recommended, "n_pcs not in recommended parameters"
        assert 'min_dist' in recommended, "min_dist not in recommended parameters"
        assert 'n_neighbors' in recommended, "n_neighbors not in recommended parameters"
        assert 'resolution' in recommended, "resolution not in recommended parameters"
        
        print("✓ Fallback analysis completed successfully")
        print(f"  - Recommended parameters: {recommended}")
        
        return analysis_results
        
    except Exception as e:
        print(f"✗ Fallback analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_final_visualization(adata, analysis_results):
    """Test final visualization generation."""
    print("\nTesting final visualization generation...")
    
    from clustering_agent import generate_final_visualization
    
    try:
        final_results = generate_final_visualization(
            adata,
            analysis_results['recommended_parameters'],
            output_dir="test_final",
            plot_width=10,
            plot_height=8
        )
        
        # Check results
        assert 'adata_final' in final_results, "Final adata not returned"
        assert 'final_plot_path' in final_results, "Final plot path not returned"
        assert 'cluster_counts' in final_results, "Cluster counts not returned"
        
        # Check that final plot was created
        assert os.path.exists(final_results['final_plot_path']), "Final plot not created"
        
        # Check final clustering
        final_adata = final_results['adata_final']
        assert 'final_clusters' in final_adata.obs, "Final clusters not in adata"
        assert 'X_umap' in final_adata.obsm, "UMAP not computed"
        
        print("✓ Final visualization generation completed successfully")
        print(f"  - Final plot: {final_results['final_plot_path']}")
        print(f"  - Number of final clusters: {len(final_results['cluster_counts'])}")
        
        return final_results
        
    except Exception as e:
        print(f"✗ Final visualization generation failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_llm_utils_image():
    """Test LLM image utilities (without actual API call)."""
    print("\nTesting LLM image utilities...")
    
    try:
        from llm_utils_image import encode_image_to_base64, create_image_url
        
        # Create a small test image
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.scatter([1, 2, 3], [1, 2, 3])
        ax.set_title("Test Image")
        test_image_path = "test_image.png"
        plt.savefig(test_image_path)
        plt.close()
        
        # Test image encoding
        base64_str = encode_image_to_base64(test_image_path)
        assert len(base64_str) > 0, "Base64 encoding failed"
        
        # Test image URL creation
        image_url = create_image_url(test_image_path)
        assert image_url.startswith("data:image/png;base64,"), "Image URL format incorrect"
        
        # Clean up
        os.remove(test_image_path)
        
        print("✓ LLM image utilities work correctly")
        return True
        
    except Exception as e:
        print(f"✗ LLM image utilities failed: {e}")
        return False


def main():
    """Run all tests."""
    print("="*60)
    print("CASSIA Clustering Agent - Test Suite")
    print("="*60)
    
    # Test imports
    if not test_imports():
        print("\n✗ Import tests failed. Cannot continue.")
        return False
    
    # Create test data
    adata = create_test_data()
    if adata is None:
        print("\n✗ Test data creation failed. Cannot continue.")
        return False
    
    # Test preprocessing
    adata_processed = test_preprocessing(adata)
    if adata_processed is None:
        print("\n✗ Preprocessing failed. Cannot continue.")
        return False
    
    # Test parameter grid
    grid_results = test_parameter_grid(adata_processed)
    if grid_results is None:
        print("\n✗ Parameter grid test failed.")
        return False
    
    # Test resolution stability
    stability_results = test_resolution_stability(adata_processed)
    if stability_results is None:
        print("\n✗ Resolution stability test failed.")
        return False
    
    # Test fallback analysis
    analysis_results = test_fallback_analysis(grid_results, stability_results)
    if analysis_results is None:
        print("\n✗ Fallback analysis test failed.")
        return False
    
    # Test final visualization
    final_results = test_final_visualization(adata_processed, analysis_results)
    if final_results is None:
        print("\n✗ Final visualization test failed.")
        return False
    
    # Test LLM utilities
    if not test_llm_utils_image():
        print("\n✗ LLM utilities test failed.")
        return False
    
    print("\n" + "="*60)
    print("🎉 ALL TESTS PASSED! 🎉")
    print("="*60)
    print("\nThe clustering agent functions are working correctly!")
    print("You can now use them in your analysis pipeline.")
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)