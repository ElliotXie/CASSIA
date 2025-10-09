"""
Sample data loading utilities for CASSIA tests.

This module provides utilities for loading and managing sample datasets.
"""

import pandas as pd
from pathlib import Path
from typing import Optional, List, Union


class SampleDataLoader:
    """
    Loader for sample CASSIA data.

    This class provides convenient methods to load various sample datasets
    included with CASSIA for testing purposes.

    Example:
        >>> loader = SampleDataLoader()
        >>> data = loader.load_processed()
        >>> print(f"Loaded {len(data)} rows")
    """

    def __init__(self):
        """Initialize loader with path to data directory."""
        # Navigate from test/shared/ to CASSIA/data/
        self.data_dir = Path(__file__).parent.parent.parent / "data"

        if not self.data_dir.exists():
            raise FileNotFoundError(
                f"Data directory not found: {self.data_dir}\n"
                "Expected structure: CASSIA_python/CASSIA/data/"
            )

    def load_processed(self) -> pd.DataFrame:
        """
        Load processed.csv (6 clean clusters for standard testing).

        Returns:
            DataFrame with processed marker data

        Example:
            >>> data = loader.load_processed()
            >>> print(data.columns.tolist())
        """
        file_path = self.data_dir / "processed.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        return pd.read_csv(file_path)

    def load_unprocessed(self) -> pd.DataFrame:
        """
        Load unprocessed.csv (large dataset for stress testing).

        Returns:
            DataFrame with unprocessed marker data

        Example:
            >>> data = loader.load_unprocessed()
            >>> print(f"Shape: {data.shape}")
        """
        file_path = self.data_dir / "unprocessed.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        return pd.read_csv(file_path)

    def load_subcluster_results(self) -> pd.DataFrame:
        """
        Load subcluster_results.csv (subclustering test data).

        Returns:
            DataFrame with subcluster results

        Example:
            >>> data = loader.load_subcluster_results()
        """
        file_path = self.data_dir / "subcluster_results.csv"
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        return pd.read_csv(file_path)

    def load_custom(self, filename: str) -> pd.DataFrame:
        """
        Load custom data file from data directory.

        Args:
            filename: Name of file to load

        Returns:
            DataFrame with loaded data

        Example:
            >>> data = loader.load_custom("my_markers.csv")
        """
        file_path = self.data_dir / filename
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")

        return pd.read_csv(file_path)

    def load_subset(
        self,
        dataset: str = "processed",
        clusters: Optional[List[str]] = None,
        n_clusters: Optional[int] = None,
        random_state: int = 42
    ) -> pd.DataFrame:
        """
        Load a subset of data.

        Args:
            dataset: 'processed', 'unprocessed', or 'subcluster'
            clusters: List of specific cluster names to include
            n_clusters: Number of random clusters to include
            random_state: Random seed for reproducibility

        Returns:
            DataFrame with subset of data

        Example:
            >>> # Get specific clusters
            >>> data = loader.load_subset(
            ...     dataset="processed",
            ...     clusters=["monocyte", "plasma cell"]
            ... )
            >>>
            >>> # Get random N clusters
            >>> data = loader.load_subset(
            ...     dataset="processed",
            ...     n_clusters=3
            ... )
        """
        # Load base dataset
        if dataset == "processed":
            df = self.load_processed()
        elif dataset == "unprocessed":
            df = self.load_unprocessed()
        elif dataset == "subcluster":
            df = self.load_subcluster_results()
        else:
            raise ValueError(
                f"Unknown dataset: {dataset}. "
                "Use 'processed', 'unprocessed', or 'subcluster'"
            )

        # If no filtering requested, return full dataset
        if clusters is None and n_clusters is None:
            return df

        # Determine cluster column name
        # Usually second column has cell types
        if len(df.columns) > 1:
            cluster_col = df.columns[1]
        else:
            raise ValueError("DataFrame must have at least 2 columns")

        # Filter by specific clusters
        if clusters is not None:
            df = df[df[cluster_col].isin(clusters)].copy()

        # OR sample N clusters
        elif n_clusters is not None:
            unique_clusters = df[cluster_col].unique()
            if n_clusters > len(unique_clusters):
                raise ValueError(
                    f"Requested {n_clusters} clusters, "
                    f"but only {len(unique_clusters)} available"
                )

            # Sample clusters
            import numpy as np
            np.random.seed(random_state)
            selected_clusters = np.random.choice(
                unique_clusters,
                size=n_clusters,
                replace=False
            )
            df = df[df[cluster_col].isin(selected_clusters)].copy()

        return df

    def get_single_cluster(
        self,
        cluster_name: str,
        dataset: str = "processed"
    ) -> pd.DataFrame:
        """
        Get data for a single cluster.

        Args:
            cluster_name: Name of cluster to extract
            dataset: Which dataset to use

        Returns:
            DataFrame with data for one cluster

        Example:
            >>> data = loader.get_single_cluster("monocyte")
        """
        return self.load_subset(dataset=dataset, clusters=[cluster_name])

    def get_marker_list(
        self,
        cluster_name: str,
        dataset: str = "processed",
        max_genes: Optional[int] = None
    ) -> List[str]:
        """
        Get marker list for a specific cluster.

        Args:
            cluster_name: Name of cluster
            dataset: Which dataset to use
            max_genes: Maximum number of genes to return

        Returns:
            List of marker gene names

        Example:
            >>> markers = loader.get_marker_list("monocyte", max_genes=50)
            >>> print(f"Top markers: {', '.join(markers[:5])}")
        """
        # Get cluster data
        cluster_data = self.get_single_cluster(cluster_name, dataset)

        if cluster_data.empty:
            return []

        # Extract marker column (usually third column)
        if len(cluster_data.columns) > 2:
            marker_col = cluster_data.columns[2]

            # Get first row (should be comma-separated markers)
            marker_string = cluster_data[marker_col].iloc[0]

            # Split by comma
            markers = [m.strip() for m in marker_string.split(',')]

            # Limit if requested
            if max_genes is not None:
                markers = markers[:max_genes]

            return markers
        else:
            return []

    def list_available_clusters(self, dataset: str = "processed") -> List[str]:
        """
        List all available cluster names in dataset.

        Args:
            dataset: Which dataset to check

        Returns:
            List of cluster names

        Example:
            >>> clusters = loader.list_available_clusters()
            >>> print(f"Available clusters: {clusters}")
        """
        df = self.load_subset(dataset=dataset)

        if len(df.columns) > 1:
            cluster_col = df.columns[1]
            return df[cluster_col].unique().tolist()
        else:
            return []

    def get_dataset_info(self, dataset: str = "processed") -> dict:
        """
        Get information about a dataset.

        Args:
            dataset: Which dataset to analyze

        Returns:
            Dictionary with dataset statistics

        Example:
            >>> info = loader.get_dataset_info("processed")
            >>> print(f"Clusters: {info['n_clusters']}")
        """
        df = self.load_subset(dataset=dataset)

        info = {
            'n_rows': len(df),
            'n_columns': len(df.columns),
            'columns': df.columns.tolist(),
        }

        if len(df.columns) > 1:
            cluster_col = df.columns[1]
            clusters = df[cluster_col].unique()
            info['n_clusters'] = len(clusters)
            info['clusters'] = clusters.tolist()

        return info


def load_sample_data(dataset: str = "processed") -> pd.DataFrame:
    """
    Quick load sample data.

    Convenience function for quickly loading a dataset without
    creating a SampleDataLoader instance.

    Args:
        dataset: 'processed', 'unprocessed', or 'subcluster'

    Returns:
        DataFrame with sample data

    Example:
        >>> data = load_sample_data("processed")
        >>> print(data.shape)
    """
    loader = SampleDataLoader()

    if dataset == "processed":
        return loader.load_processed()
    elif dataset == "unprocessed":
        return loader.load_unprocessed()
    elif dataset == "subcluster":
        return loader.load_subcluster_results()
    else:
        raise ValueError(
            f"Unknown dataset: {dataset}. "
            "Use 'processed', 'unprocessed', or 'subcluster'"
        )


def load_marker_list(
    cluster_name: str,
    dataset: str = "processed",
    max_genes: Optional[int] = 50
) -> List[str]:
    """
    Quick load marker list for a cluster.

    Args:
        cluster_name: Name of cluster
        dataset: Which dataset to use
        max_genes: Maximum number of genes

    Returns:
        List of marker genes

    Example:
        >>> markers = load_marker_list("monocyte", max_genes=30)
    """
    loader = SampleDataLoader()
    return loader.get_marker_list(cluster_name, dataset, max_genes)
