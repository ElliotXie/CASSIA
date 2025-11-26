#!/usr/bin/env python3
"""
CASSIA Clustering Agent - Automated Parameter Optimization for Single-Cell Clustering

This module provides intelligent clustering functions that:
1. Generate parameter grids with different clustering parameters
2. Test resolution stability 
3. Use LLM to analyze plots and select optimal parameters
4. Generate final visualizations

Author: CASSIA Team
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
import json
import warnings
from typing import Dict, List, Tuple, Optional, Union
from itertools import product
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor='white')


def generate_parameter_grid(
    adata,
    output_dir: str = "clustering_grid",
    n_pcs_list: List[int] = [20, 30, 40],
    min_dist_list: List[float] = [0.2, 0.3],
    n_neighbors_list: List[int] = [30, 40],
    resolution: float = 0.7,
    plot_width: int = 8,
    plot_height: int = 6
) -> Dict:
    """
    Generate parameter grid plots for different clustering parameter combinations.
    
    Parameters:
    -----------
    adata : AnnData
        Preprocessed AnnData object (normalized, scaled, with PCA)
    output_dir : str
        Directory to save plots
    n_pcs_list : List[int]
        List of PC dimensions to test
    min_dist_list : List[float]
        List of min_dist values for UMAP
    n_neighbors_list : List[int]
        List of n_neighbors values for UMAP
    resolution : float
        Fixed clustering resolution for parameter grid
    plot_width : int
        Plot width in inches
    plot_height : int
        Plot height in inches
        
    Returns:
    --------
    Dict containing:
        - parameter_combinations: List of all parameter combinations tested
        - plot_paths: List of generated plot paths
        - combined_plot_path: Path to combined grid plot
        - results_summary: Summary of results for each combination
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Check preprocessing
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA must be computed before generating parameter grid. Run sc.pp.pca() first.")
    
    # Generate all parameter combinations
    param_combinations = list(product(n_pcs_list, min_dist_list, n_neighbors_list))
    
    print(f"Generating parameter grid with {len(param_combinations)} combinations...")
    
    # Store results
    results_summary = []
    plot_paths = []
    
    # Process each combination
    for i, (n_pcs, min_dist, n_neighbors) in enumerate(param_combinations):
        print(f"Processing {i+1}/{len(param_combinations)}: "
              f"n_pcs={n_pcs}, min_dist={min_dist}, n_neighbors={n_neighbors}")
        
        # Create copy for this combination
        adata_temp = adata.copy()
        
        # Clustering pipeline
        sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.leiden(adata_temp, resolution=resolution, key_added='clusters')
        sc.tl.umap(adata_temp, min_dist=min_dist)
        
        # Generate plot
        reduction_name = f"n_pcs_{n_pcs}_min_dist_{min_dist}_n_neighbors_{n_neighbors}"
        plot_path = os.path.join(output_dir, f"{reduction_name}.png")
        
        fig, ax = plt.subplots(figsize=(plot_width, plot_height))
        sc.pl.umap(adata_temp, color='clusters', legend_loc='on data', 
                   legend_fontsize=12, legend_fontoutline=2, ax=ax, show=False)
        
        ax.set_title(f"n_pcs: {n_pcs} | min_dist: {min_dist} | n_neighbors: {n_neighbors}\n"
                    f"Clusters: {len(adata_temp.obs['clusters'].cat.categories)}", 
                    fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        # Store results
        n_clusters = len(adata_temp.obs['clusters'].cat.categories)
        results_summary.append({
            'n_pcs': n_pcs,
            'min_dist': min_dist,
            'n_neighbors': n_neighbors,
            'n_clusters': n_clusters,
            'plot_path': plot_path,
            'reduction_name': reduction_name
        })
        plot_paths.append(plot_path)
    
    # Create combined plot
    print("Creating combined parameter grid plot...")
    n_combinations = len(param_combinations)
    n_cols = min(4, n_combinations)
    n_rows = (n_combinations + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 5))
    if n_rows == 1:
        axes = axes.reshape(1, -1) if n_cols > 1 else [axes]
    if n_cols == 1:
        axes = axes.reshape(-1, 1) if n_rows > 1 else [axes]
    
    for i, result in enumerate(results_summary):
        row = i // n_cols
        col = i % n_cols
        
        # Load and display individual plot
        img = plt.imread(result['plot_path'])
        if n_rows == 1 and n_cols == 1:
            axes.imshow(img)
            axes.axis('off')
            axes.set_title(f"{result['reduction_name']}\n({result['n_clusters']} clusters)", 
                          fontsize=10)
        else:
            axes[row, col].imshow(img)
            axes[row, col].axis('off')
            axes[row, col].set_title(f"{result['reduction_name']}\n({result['n_clusters']} clusters)", 
                                    fontsize=10)
    
    # Hide unused subplots
    for i in range(n_combinations, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        if n_rows > 1 or n_cols > 1:
            axes[row, col].axis('off')
    
    plt.tight_layout()
    combined_plot_path = os.path.join(output_dir, "parameter_grid_combined.png")
    plt.savefig(combined_plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Parameter grid completed. Combined plot: {combined_plot_path}")
    
    return {
        'parameter_combinations': param_combinations,
        'plot_paths': plot_paths,
        'combined_plot_path': combined_plot_path,
        'results_summary': results_summary,
        'output_dir': output_dir
    }


def test_resolution_stability(
    adata,
    output_dir: str = "resolution_stability",
    resolution_list: List[float] = [0.3, 0.5, 0.7, 0.9, 1.1, 1.3],
    n_pcs: int = 30,
    min_dist: float = 0.2,
    n_neighbors: int = 30,
    n_iterations: int = 5,
    plot_width: int = 10,
    plot_height: int = 6
) -> Dict:
    """
    Test clustering stability across different resolutions.
    
    Parameters:
    -----------
    adata : AnnData
        Preprocessed AnnData object
    output_dir : str
        Directory to save stability analysis
    resolution_list : List[float]
        List of resolutions to test
    n_pcs : int
        Number of PCs to use
    min_dist : float
        UMAP min_dist parameter
    n_neighbors : int
        Number of neighbors
    n_iterations : int
        Number of iterations for stability testing
    plot_width : int
        Plot width
    plot_height : int
        Plot height
        
    Returns:
    --------
    Dict containing stability analysis results
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Testing resolution stability with {len(resolution_list)} resolutions...")
    
    stability_results = []
    plot_paths = []
    
    for resolution in resolution_list:
        print(f"Testing resolution {resolution}...")
        
        # Test stability with multiple iterations
        clustering_results = []
        n_clusters_list = []
        
        for iteration in range(n_iterations):
            adata_temp = adata.copy()
            
            # Add small amount of noise for stability testing
            if iteration > 0:
                noise_factor = 0.01
                adata_temp.X = adata_temp.X + np.random.normal(0, noise_factor, adata_temp.X.shape)
            
            # Clustering pipeline
            sc.pp.neighbors(adata_temp, n_neighbors=n_neighbors, n_pcs=n_pcs)
            sc.tl.leiden(adata_temp, resolution=resolution, key_added='clusters')
            
            clustering_results.append(adata_temp.obs['clusters'].values)
            n_clusters_list.append(len(adata_temp.obs['clusters'].cat.categories))
        
        # Calculate stability metrics
        stability_scores = []
        for i in range(len(clustering_results)):
            for j in range(i+1, len(clustering_results)):
                ari = adjusted_rand_score(clustering_results[i], clustering_results[j])
                stability_scores.append(ari)
        
        mean_stability = np.mean(stability_scores)
        std_stability = np.std(stability_scores)
        mean_n_clusters = np.mean(n_clusters_list)
        std_n_clusters = np.std(n_clusters_list)
        
        # Generate representative clustering for visualization
        adata_repr = adata.copy()
        sc.pp.neighbors(adata_repr, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.leiden(adata_repr, resolution=resolution, key_added='clusters')
        sc.tl.umap(adata_repr, min_dist=min_dist)
        
        # Plot
        plot_path = os.path.join(output_dir, f"resolution_{resolution}.png")
        fig, ax = plt.subplots(figsize=(plot_width, plot_height))
        sc.pl.umap(adata_repr, color='clusters', legend_loc='on data', 
                   legend_fontsize=12, legend_fontoutline=2, ax=ax, show=False)
        
        ax.set_title(f"Resolution: {resolution}\n"
                    f"Clusters: {len(adata_repr.obs['clusters'].cat.categories)} "
                    f"(±{std_n_clusters:.1f}), Stability: {mean_stability:.3f} (±{std_stability:.3f})", 
                    fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        stability_results.append({
            'resolution': resolution,
            'mean_stability': mean_stability,
            'std_stability': std_stability,
            'mean_n_clusters': mean_n_clusters,
            'std_n_clusters': std_n_clusters,
            'plot_path': plot_path
        })
        plot_paths.append(plot_path)
    
    # Create stability summary plot
    resolutions = [r['resolution'] for r in stability_results]
    stabilities = [r['mean_stability'] for r in stability_results]
    stability_stds = [r['std_stability'] for r in stability_results]
    cluster_counts = [r['mean_n_clusters'] for r in stability_results]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Stability plot
    ax1.errorbar(resolutions, stabilities, yerr=stability_stds, 
                marker='o', capsize=5, capthick=2, linewidth=2)
    ax1.set_xlabel('Resolution')
    ax1.set_ylabel('Stability (ARI)')
    ax1.set_title('Clustering Stability vs Resolution')
    ax1.grid(True, alpha=0.3)
    
    # Cluster count plot
    ax2.plot(resolutions, cluster_counts, marker='s', linewidth=2, markersize=8)
    ax2.set_xlabel('Resolution')
    ax2.set_ylabel('Number of Clusters')
    ax2.set_title('Number of Clusters vs Resolution')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    stability_summary_path = os.path.join(output_dir, "stability_summary.png")
    plt.savefig(stability_summary_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Create combined resolution plots
    n_resolutions = len(resolution_list)
    n_cols = min(3, n_resolutions)
    n_rows = (n_resolutions + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 6, n_rows * 5))
    if n_rows == 1:
        axes = axes.reshape(1, -1) if n_cols > 1 else [axes]
    if n_cols == 1:
        axes = axes.reshape(-1, 1) if n_rows > 1 else [axes]
    
    for i, result in enumerate(stability_results):
        row = i // n_cols
        col = i % n_cols
        
        img = plt.imread(result['plot_path'])
        if n_rows == 1 and n_cols == 1:
            axes.imshow(img)
            axes.axis('off')
        else:
            axes[row, col].imshow(img)
            axes[row, col].axis('off')
    
    # Hide unused subplots
    for i in range(n_resolutions, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        if n_rows > 1 or n_cols > 1:
            axes[row, col].axis('off')
    
    plt.tight_layout()
    combined_resolution_path = os.path.join(output_dir, "resolution_grid_combined.png")
    plt.savefig(combined_resolution_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Resolution stability analysis completed. Summary: {stability_summary_path}")
    
    return {
        'stability_results': stability_results,
        'plot_paths': plot_paths,
        'stability_summary_path': stability_summary_path,
        'combined_resolution_path': combined_resolution_path,
        'output_dir': output_dir
    }


def analyze_clustering_with_llm(
    parameter_grid_results: Dict,
    resolution_stability_results: Dict,
    model: str = "google/gemini-2.5-flash",
    provider: str = "openrouter",
    api_key: Optional[str] = None,
    max_images_per_call: int = 6
) -> Dict:
    """
    Use LLM to analyze clustering results and select optimal parameters.
    
    Parameters:
    -----------
    parameter_grid_results : Dict
        Results from generate_parameter_grid()
    resolution_stability_results : Dict
        Results from test_resolution_stability()
    model : str
        LLM model to use
    provider : str
        API provider
    api_key : Optional[str]
        API key
    max_images_per_call : int
        Maximum images per LLM call
        
    Returns:
    --------
    Dict containing LLM analysis and recommended parameters
    """
    
    print("Analyzing clustering results with LLM...")
    
    try:
        # Import image analysis utilities
        from llm_utils_image import analyze_clustering_images
        
        # Analyze parameter grid
        parameter_analysis = analyze_clustering_images(
            image_paths=parameter_grid_results['plot_paths'],
            combined_image_path=parameter_grid_results['combined_plot_path'],
            results_data=parameter_grid_results['results_summary'],
            analysis_type="parameter_grid",
            model=model,
            provider=provider,
            api_key=api_key,
            max_images_per_call=max_images_per_call
        )
        
        # Analyze resolution stability
        resolution_analysis = analyze_clustering_images(
            image_paths=resolution_stability_results['plot_paths'],
            combined_image_path=resolution_stability_results['combined_resolution_path'],
            results_data=resolution_stability_results['stability_results'],
            analysis_type="resolution_stability",
            model=model,
            provider=provider,
            api_key=api_key,
            max_images_per_call=max_images_per_call
        )
        
        # Combine analyses to get final recommendations
        final_analysis = analyze_clustering_images(
            image_paths=[parameter_grid_results['combined_plot_path'], 
                        resolution_stability_results['stability_summary_path']],
            combined_image_path=None,
            results_data={
                'parameter_results': parameter_grid_results['results_summary'],
                'stability_results': resolution_stability_results['stability_results']
            },
            analysis_type="final_recommendation",
            model=model,
            provider=provider,
            api_key=api_key,
            max_images_per_call=max_images_per_call
        )
        
        return {
            'parameter_analysis': parameter_analysis,
            'resolution_analysis': resolution_analysis,
            'final_recommendation': final_analysis,
            'recommended_parameters': final_analysis.get('recommended_parameters', {}),
            'analysis_method': 'llm'
        }
        
    except ImportError:
        print("Warning: llm_utils_image not found. Using fallback analysis.")
        return _fallback_analysis(parameter_grid_results, resolution_stability_results)


def generate_final_visualization(
    adata,
    recommended_parameters: Dict,
    output_dir: str = "final_clustering",
    plot_width: int = 12,
    plot_height: int = 8
) -> Dict:
    """
    Generate final clustering visualization with recommended parameters.
    
    Parameters:
    -----------
    adata : AnnData
        Preprocessed AnnData object
    recommended_parameters : Dict
        Recommended parameters from LLM analysis
    output_dir : str
        Output directory
    plot_width : int
        Plot width
    plot_height : int
        Plot height
        
    Returns:
    --------
    Dict containing final clustering results
    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Generating final clustering with parameters: {recommended_parameters}")
    
    # Apply recommended parameters
    adata_final = adata.copy()
    
    sc.pp.neighbors(adata_final, 
                    n_neighbors=recommended_parameters['n_neighbors'],
                    n_pcs=recommended_parameters['n_pcs'])
    
    sc.tl.leiden(adata_final, 
                 resolution=recommended_parameters['resolution'],
                 key_added='final_clusters')
    
    sc.tl.umap(adata_final, min_dist=recommended_parameters['min_dist'])
    
    # Generate comprehensive visualization
    fig, axes = plt.subplots(2, 2, figsize=(plot_width*2, plot_height*2))
    
    # Main clustering plot
    sc.pl.umap(adata_final, color='final_clusters', legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=2, ax=axes[0,0], show=False)
    axes[0,0].set_title('Final Clustering Results', fontsize=14, fontweight='bold')
    
    # Cluster proportions
    cluster_counts = adata_final.obs['final_clusters'].value_counts().sort_index()
    axes[0,1].bar(range(len(cluster_counts)), cluster_counts.values)
    axes[0,1].set_xlabel('Cluster')
    axes[0,1].set_ylabel('Number of Cells')
    axes[0,1].set_title('Cluster Sizes', fontsize=14, fontweight='bold')
    
    # Parameter summary
    param_text = f"""Final Parameters:
    
n_pcs: {recommended_parameters['n_pcs']}
min_dist: {recommended_parameters['min_dist']}
n_neighbors: {recommended_parameters['n_neighbors']}
resolution: {recommended_parameters['resolution']}

Total clusters: {len(cluster_counts)}
Total cells: {adata_final.n_obs}
"""
    axes[1,0].text(0.1, 0.5, param_text, fontsize=12, verticalalignment='center',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
    axes[1,0].set_xlim(0, 1)
    axes[1,0].set_ylim(0, 1)
    axes[1,0].axis('off')
    axes[1,0].set_title('Parameter Summary', fontsize=14, fontweight='bold')
    
    # Quality metrics (if available)
    try:
        # Calculate some basic quality metrics
        # Silhouette score, connectivity, etc.
        from sklearn.metrics import silhouette_score
        if hasattr(adata_final.obsm, 'X_pca'):
            silhouette = silhouette_score(adata_final.obsm['X_pca'][:, :recommended_parameters['n_pcs']], 
                                        adata_final.obs['final_clusters'])
            
            quality_text = f"""Quality Metrics:
            
Silhouette Score: {silhouette:.3f}
Number of Clusters: {len(cluster_counts)}
Median Cluster Size: {cluster_counts.median():.0f}
"""
            axes[1,1].text(0.1, 0.5, quality_text, fontsize=12, verticalalignment='center',
                          bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        else:
            axes[1,1].text(0.1, 0.5, "Quality metrics\nnot available", fontsize=12, 
                          verticalalignment='center')
    except Exception:
        axes[1,1].text(0.1, 0.5, "Quality metrics\nnot available", fontsize=12,
                      verticalalignment='center')
    
    axes[1,1].set_xlim(0, 1)
    axes[1,1].set_ylim(0, 1)
    axes[1,1].axis('off')
    axes[1,1].set_title('Quality Assessment', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    final_plot_path = os.path.join(output_dir, "final_clustering_summary.png")
    plt.savefig(final_plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save the final clustering results
    results_path = os.path.join(output_dir, "final_clustering_results.h5ad")
    adata_final.write(results_path)
    
    print(f"Final clustering completed. Results saved to: {results_path}")
    
    return {
        'adata_final': adata_final,
        'final_plot_path': final_plot_path,
        'results_path': results_path,
        'cluster_counts': cluster_counts.to_dict(),
        'recommended_parameters': recommended_parameters
    }


def _fallback_analysis(parameter_grid_results: Dict, resolution_stability_results: Dict) -> Dict:
    """
    Fallback analysis when LLM is not available.
    """
    
    # Select parameters based on heuristics
    param_results = parameter_grid_results['results_summary']
    stability_results = resolution_stability_results['stability_results']
    
    # Find most stable resolution
    best_stability_idx = np.argmax([r['mean_stability'] for r in stability_results])
    best_resolution = stability_results[best_stability_idx]['resolution']
    
    # Find moderate cluster count from parameter grid
    cluster_counts = [r['n_clusters'] for r in param_results]
    median_clusters = np.median(cluster_counts)
    best_param_idx = np.argmin([abs(count - median_clusters) for count in cluster_counts])
    best_params = param_results[best_param_idx]
    
    return {
        'recommended_parameters': {
            'n_pcs': best_params['n_pcs'],
            'min_dist': best_params['min_dist'],
            'n_neighbors': best_params['n_neighbors'],
            'resolution': best_resolution
        },
        'analysis_method': 'heuristic_fallback',
        'reasoning': f"Selected most stable resolution ({best_resolution}) and "
                    f"parameters producing moderate cluster count ({best_params['n_clusters']})"
    }


def preprocess_for_clustering(adata, n_top_genes: int = 2000, n_pcs: int = 50):
    """
    Preprocessing pipeline for clustering analysis.
    """
    
    print("Preprocessing data for clustering...")
    
    adata = adata.copy()
    
    # Basic preprocessing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    
    print(f"Preprocessing complete. Shape: {adata.shape}")
    
    return adata


if __name__ == "__main__":
    print("CASSIA Clustering Agent - Modular Functions")
    print("Use individual functions in your analysis pipeline:")
    print("1. generate_parameter_grid()")
    print("2. test_resolution_stability()")
    print("3. analyze_clustering_with_llm()")
    print("4. generate_final_visualization()")