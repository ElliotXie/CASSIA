"""
Generate Scanpy rank_genes_groups Output for Testing
=====================================================
This script generates and saves Scanpy's rank_genes_groups output from the
PBMC3k dataset for use in CASSIA testing.

Run this script once to generate the test data:
    python generate_scanpy_markers.py

Output files saved to ./data/:
    - scanpy_rank_genes_groups.pkl  (the rank_genes_groups dictionary - structured array format)
    - scanpy_cluster_info.pkl       (cluster column name and cell counts)
    - scanpy_rank_genes_df.csv      (flat DataFrame from sc.get.rank_genes_groups_df)
"""

import pickle
import numpy as np
import pandas as pd
from pathlib import Path


def generate_and_save_scanpy_markers():
    """Generate Scanpy rank_genes_groups from PBMC3k and save to files."""
    import scanpy as sc

    print("=" * 60)
    print("Generating Scanpy rank_genes_groups test data")
    print("=" * 60)

    # Create data directory
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)

    # Load PBMC3k dataset
    print("\nLoading PBMC3k dataset from Scanpy...")
    adata = sc.datasets.pbmc3k_processed()

    print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"Cluster column: 'louvain' with {adata.obs['louvain'].nunique()} clusters")
    print(f"Clusters: {adata.obs['louvain'].unique().tolist()}")

    # Run rank_genes_groups
    print("\nRunning sc.tl.rank_genes_groups (Wilcoxon test)...")
    sc.tl.rank_genes_groups(
        adata,
        groupby='louvain',
        method='wilcoxon',
        use_raw=False,
        key_added='rank_genes_groups'
    )

    # Extract results
    rank_genes_dict = dict(adata.uns['rank_genes_groups'])
    print(f"rank_genes_groups completed!")
    print(f"Keys in result: {list(rank_genes_dict.keys())}")

    # Calculate pcts (percentage expressed) for each cluster
    print("\nCalculating percentage expressed per cluster...")
    cluster_col = 'louvain'
    clusters = rank_genes_dict['names'].dtype.names
    n_genes = len(rank_genes_dict['names'])

    # Create dtype for pcts
    dt_float = np.dtype([(name, 'f8') for name in clusters])
    pcts_data = np.zeros(n_genes, dtype=dt_float)

    # Get the expression matrix
    if adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.X
        var_names = adata.var_names

    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()

    var_names_list = list(var_names)

    for cluster in clusters:
        # Get cells in this cluster
        cluster_mask = adata.obs[cluster_col].astype(str) == str(cluster)
        X_cluster = X[cluster_mask, :]

        # For each gene in the ranked list
        for i, gene in enumerate(rank_genes_dict['names'][cluster]):
            if gene in var_names_list:
                gene_idx = var_names_list.index(gene)
                # Calculate percentage of cells expressing this gene (> 0)
                pct = np.mean(X_cluster[:, gene_idx] > 0)
                pcts_data[cluster][i] = pct
            else:
                pcts_data[cluster][i] = 0.0

    # Add pcts to the dictionary
    rank_genes_dict['pcts'] = pcts_data

    # Save rank_genes_groups dictionary
    rank_genes_path = data_dir / "scanpy_rank_genes_groups.pkl"
    with open(rank_genes_path, 'wb') as f:
        pickle.dump(rank_genes_dict, f)
    print(f"\nSaved rank_genes_groups to: {rank_genes_path}")

    # Save cluster info
    cluster_info = {
        'cluster_col': cluster_col,
        'clusters': list(clusters),
        'n_cells': adata.shape[0],
        'n_genes': adata.shape[1],
        'cell_counts': adata.obs[cluster_col].value_counts().to_dict()
    }
    cluster_info_path = data_dir / "scanpy_cluster_info.pkl"
    with open(cluster_info_path, 'wb') as f:
        pickle.dump(cluster_info, f)
    print(f"Saved cluster info to: {cluster_info_path}")

    # Also save the flat DataFrame format from sc.get.rank_genes_groups_df()
    # This is the modern scanpy API that returns a flat DataFrame
    print("\nGenerating flat DataFrame format (sc.get.rank_genes_groups_df)...")
    scanpy_df = sc.get.rank_genes_groups_df(adata, group=None)
    scanpy_df_path = data_dir / "scanpy_rank_genes_df.csv"
    scanpy_df.to_csv(scanpy_df_path, index=False)
    print(f"Saved flat DataFrame to: {scanpy_df_path}")
    print(f"  Shape: {scanpy_df.shape}")
    print(f"  Columns: {scanpy_df.columns.tolist()}")
    print(f"  Sample rows:")
    print(scanpy_df.head(3).to_string(index=False))

    # Print summary
    print("\n" + "=" * 60)
    print("Summary of saved data:")
    print("=" * 60)
    print(f"Clusters: {clusters}")
    print(f"Genes per cluster: {n_genes}")
    print(f"\nCell counts per cluster:")
    for cluster, count in cluster_info['cell_counts'].items():
        print(f"  {cluster}: {count} cells")

    print(f"\nKeys in rank_genes_groups: {list(rank_genes_dict.keys())}")

    # Verify the data
    print("\n" + "=" * 60)
    print("Verification - Top 5 markers per cluster:")
    print("=" * 60)
    for cluster in clusters:
        genes = rank_genes_dict['names'][cluster][:5]
        scores = rank_genes_dict['scores'][cluster][:5]
        logfc = rank_genes_dict['logfoldchanges'][cluster][:5]
        pcts = rank_genes_dict['pcts'][cluster][:5]
        print(f"\n{cluster}:")
        for i in range(5):
            print(f"  {genes[i]}: score={scores[i]:.2f}, logFC={logfc[i]:.2f}, pct={pcts[i]:.2f}")

    print("\n" + "=" * 60)
    print("Data generation complete!")
    print("=" * 60)
    print("\nGenerated files:")
    print(f"  1. {rank_genes_path.name} - Structured array format (legacy)")
    print(f"  2. {cluster_info_path.name} - Cluster metadata")
    print(f"  3. {scanpy_df_path.name} - Flat DataFrame format (modern)")

    return rank_genes_dict, cluster_info, scanpy_df


if __name__ == "__main__":
    generate_and_save_scanpy_markers()
