"""
CASSIA Test Suite - Data Fixtures
=================================
Functions for loading and preparing test data.
"""

import pandas as pd
from pathlib import Path


def get_test_root() -> Path:
    """Get the root directory of the test suite."""
    return Path(__file__).parent.parent.parent


def get_marker_file_path() -> Path:
    """Get the path to the marker data file."""
    return get_test_root() / "data" / "markers" / "processed.csv"


def load_markers() -> pd.DataFrame:
    """
    Load the processed marker data.

    Returns:
        pd.DataFrame: DataFrame with columns ['Broad.cell.type', 'Top.Markers']
    """
    marker_path = get_marker_file_path()
    if not marker_path.exists():
        raise FileNotFoundError(f"Marker file not found: {marker_path}")

    df = pd.read_csv(marker_path, index_col=0)
    return df


def get_all_clusters() -> list:
    """
    Get list of all available cell type clusters.

    Returns:
        list: List of cluster names
    """
    df = load_markers()
    return df['Broad.cell.type'].tolist()


def get_cluster_markers(cluster_name: str, n_genes: int = None) -> list:
    """
    Get marker genes for a specific cluster.

    Args:
        cluster_name: Name of the cell type cluster
        n_genes: Number of top genes to return (None = all)

    Returns:
        list: List of marker gene names
    """
    df = load_markers()

    # Find the cluster
    cluster_row = df[df['Broad.cell.type'] == cluster_name]
    if cluster_row.empty:
        available = get_all_clusters()
        raise ValueError(f"Cluster '{cluster_name}' not found. Available: {available}")

    # Get markers string and split
    markers_str = cluster_row['Top.Markers'].values[0]
    markers = [m.strip() for m in markers_str.split(',')]

    if n_genes is not None:
        markers = markers[:n_genes]

    return markers


def get_marker_dataframe_for_cluster(cluster_name: str, n_genes: int = 30) -> pd.DataFrame:
    """
    Get a DataFrame formatted for runCASSIA for a single cluster.

    Args:
        cluster_name: Name of the cell type cluster
        n_genes: Number of top genes to include

    Returns:
        pd.DataFrame: DataFrame with columns ['gene', 'cell_type']
    """
    markers = get_cluster_markers(cluster_name, n_genes)
    return pd.DataFrame({
        'gene': markers,
        'cell_type': [cluster_name] * len(markers)
    })


def get_full_marker_dataframe() -> pd.DataFrame:
    """
    Get the full marker DataFrame for batch processing.

    Returns:
        pd.DataFrame: DataFrame with all clusters and their markers
    """
    return load_markers()
