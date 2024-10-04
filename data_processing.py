import scanpy as sc
import numpy as np
import pandas as pd

def process_tissue_data(adata, tissue_type):
    """
    Process the AnnData object for a specific tissue type.
    
    Parameters:
    adata (AnnData): The original AnnData object containing all tissues.
    tissue_type (str): The tissue type to process (e.g., 'prostate', 'lung').
    
    Returns:
    AnnData: Processed AnnData object for the specified tissue.
    """
    # Subset the AnnData object to include only the specified tissue
    tissue_adata = adata[adata.obs['tissue'] == tissue_type].copy()
    
    # Print unique values for different cell type classifications
    print(f"{tissue_type.capitalize()} Broad cell types:")
    print(tissue_adata.obs['Broad cell type'].unique())
    print(f"\n{tissue_type.capitalize()} Granular cell types:")
    print(tissue_adata.obs['Granular cell type'].unique())
    print(f"\n{tissue_type.capitalize()} Cell types level 2:")
    print(tissue_adata.obs['Cell types level 2'].unique())
    print(f"\n{tissue_type.capitalize()} Cell types level 3:")
    print(tissue_adata.obs['Cell types level 3'].unique())
    
    # Normalize the data by total counts per cell and scale to 10,000 reads per cell
    sc.pp.normalize_total(tissue_adata, target_sum=1e4)
    
    # Log-transform the data after adding a pseudocount of 1
    sc.pp.log1p(tissue_adata)
    
    # Perform batch correction
    sc.pp.combat(tissue_adata, key='batch')
    
    return tissue_adata

def process_multiple_tissues(adata, tissues, n_genes=50):
    """
    Process multiple tissues and generate marker genes for each.
    
    Parameters:
    adata (AnnData): The original AnnData object containing all tissues.
    tissues (list): List of tissue names to process.
    n_genes (int): Number of top genes to include in the analysis (default: 50).
    
    Returns:
    dict: A dictionary containing results for each tissue.
    """
    annotation_levels = ['Broad cell type', 'Granular cell type', 'Cell types level 2', 'Cell types level 3']
    results = {}

    for tissue in tissues:
        print(f"\nProcessing {tissue}...")
        processed_adata = process_tissue_data(adata, tissue)
        tissue_results = analyze_and_export_markers(processed_adata, annotation_levels, n_genes=n_genes)
        results[tissue] = tissue_results
        print(f"Processing and analysis completed for {tissue}.")

    return results