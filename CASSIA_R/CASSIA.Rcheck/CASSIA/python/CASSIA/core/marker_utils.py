"""
Marker processing utilities for CASSIA.

This module provides functions for loading, parsing, and processing
marker gene data from various formats (Seurat, Scanpy).
"""

import re
import os
import pandas as pd
import numpy as np
from importlib import resources


def split_markers(marker_string):
    """
    Split a marker string into individual marker names.

    Handles various separator formats: comma+space, comma only, or space.

    Args:
        marker_string (str): String containing marker names

    Returns:
        list: List of individual marker names
    """
    # First, try splitting by comma and space
    markers = re.split(r',\s*', marker_string)

    # If that results in only one marker, try splitting by comma only
    if len(markers) == 1:
        markers = marker_string.split(',')

    # If still only one marker, try splitting by space
    if len(markers) == 1:
        markers = marker_string.split()

    # Remove any empty strings
    markers = [m.strip() for m in markers if m.strip()]

    return markers


def _validate_ranking_parameters(df, ranking_method, ascending):
    """
    Validate ranking method and column existence.

    Args:
        df: DataFrame to validate
        ranking_method (str): Ranking method to use
        ascending: Sort direction

    Raises:
        ValueError: If ranking method is invalid or required column is missing
    """
    valid_methods = ["avg_log2FC", "p_val_adj", "pct_diff", "Score"]
    if ranking_method not in valid_methods:
        raise ValueError(f"ranking_method must be one of {valid_methods}")

    if ranking_method == "Score" and "Score" not in df.columns:
        available_cols = [col for col in df.columns if col.lower() in ['score', 'scores']]
        if available_cols:
            suggestion = f". Did you mean '{available_cols[0]}'?"
        else:
            suggestion = ". Available numeric columns: " + ", ".join(df.select_dtypes(include=[np.number]).columns.tolist())
        raise ValueError(f"Column 'Score' not found in DataFrame{suggestion}")

    if ranking_method == "p_val_adj" and "p_val_adj" not in df.columns:
        raise ValueError("Column 'p_val_adj' not found in DataFrame")


def _prepare_ranking_column(df, ranking_method):
    """
    Prepare the ranking column, calculating if necessary.

    Args:
        df: DataFrame to prepare
        ranking_method (str): Ranking method to use

    Returns:
        tuple: (prepared DataFrame, column name to sort by)
    """
    df_copy = df.copy()

    if ranking_method == "pct_diff":
        if "pct.1" not in df_copy.columns or "pct.2" not in df_copy.columns:
            raise ValueError("Columns 'pct.1' and 'pct.2' required for pct_diff ranking")
        df_copy["pct_diff"] = df_copy["pct.1"] - df_copy["pct.2"]
        return df_copy, "pct_diff"

    return df_copy, ranking_method


def _get_sort_direction(ranking_method, ascending):
    """
    Get the sort direction for the ranking method.

    Args:
        ranking_method (str): Ranking method being used
        ascending: User-specified sort direction, or None for default

    Returns:
        bool: True for ascending, False for descending
    """
    DEFAULT_SORT_DIRECTIONS = {
        "avg_log2FC": False,    # Higher is better
        "p_val_adj": True,      # Lower p-value is better
        "pct_diff": False,      # Higher difference is better
        "Score": False          # Higher score is better (default)
    }

    if ascending is not None:
        return ascending
    return DEFAULT_SORT_DIRECTIONS.get(ranking_method, False)


def get_top_markers(df, n_genes=10, format_type=None, ranking_method="avg_log2FC", ascending=None):
    """
    Get top markers from either Seurat or Scanpy differential expression results.

    Args:
        df: Either a pandas DataFrame (Seurat format, or Scanpy flat DataFrame from
            sc.get.rank_genes_groups_df), or dictionary (Scanpy structured array format)
        n_genes: Number of top genes to return per cluster
        format_type: Either 'seurat', 'scanpy', 'scanpy_df', or None (auto-detect)
        ranking_method: Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending: Sort direction (None uses default for each method)

    Returns:
        pandas DataFrame with cluster and marker columns
    """
    # Auto-detect format if not specified
    if format_type is None:
        if isinstance(df, pd.DataFrame):
            # Check for scanpy flat DataFrame format from sc.get.rank_genes_groups_df()
            # Has columns: group, names, scores, logfoldchanges, pvals, pvals_adj
            if 'group' in df.columns and 'names' in df.columns and 'logfoldchanges' in df.columns:
                format_type = 'scanpy_df'
            # Check for Seurat format (has cluster and gene columns)
            elif 'cluster' in df.columns and 'gene' in df.columns:
                format_type = 'seurat'
            # Fallback: try scanpy structured array format
            elif 'names' in df and hasattr(df['names'], 'dtype') and hasattr(df['names'].dtype, 'names') and df['names'].dtype.names is not None:
                format_type = 'scanpy'
            else:
                # Default to seurat if it looks like a DataFrame
                format_type = 'seurat'
        elif hasattr(df, '__getitem__') and 'names' in df:
            # Dictionary-like structure (scanpy structured array)
            format_type = 'scanpy'
        else:
            format_type = 'seurat'

    if format_type == 'scanpy_df':
        # Process Scanpy flat DataFrame format from sc.get.rank_genes_groups_df()
        # Columns: group, names, scores, logfoldchanges, pvals, pvals_adj
        print("Detected scanpy DataFrame format (from sc.get.rank_genes_groups_df)")

        # Create a copy to avoid modifying the original
        df_work = df.copy()

        # Rename columns to match Seurat format
        column_mapping = {
            'group': 'cluster',
            'names': 'gene',
            'logfoldchanges': 'avg_log2FC',
            'pvals_adj': 'p_val_adj'
        }
        df_work = df_work.rename(columns=column_mapping)

        # Handle inf values in avg_log2FC
        df_work['avg_log2FC'] = pd.to_numeric(df_work['avg_log2FC'], errors='coerce')
        max_non_inf = df_work['avg_log2FC'].replace([np.inf, -np.inf], np.nan).max()
        min_non_inf = df_work['avg_log2FC'].replace([np.inf, -np.inf], np.nan).min()
        df_work['avg_log2FC'] = df_work['avg_log2FC'].replace([np.inf, -np.inf], [max_non_inf, min_non_inf])

        # Check if PCT columns exist (from enhance_scanpy_markers which produces pct.1/pct.2)
        has_pct = 'pct.1' in df.columns
        if has_pct:
            # Copy pct columns to working dataframe
            df_work['pct.1'] = df['pct.1']
            if 'pct.2' in df.columns:
                df_work['pct.2'] = df['pct.2']
            else:
                df_work['pct.2'] = 0
            # Filter with PCT
            df_filtered = df_work[
                (df_work['p_val_adj'] < 0.05) &
                (df_work['avg_log2FC'] > 0.25) &
                (df_work['pct.1'] >= 0.1)
            ].copy()
        else:
            # No PCT columns - filter without PCT threshold
            df_filtered = df_work[
                (df_work['p_val_adj'] < 0.05) &
                (df_work['avg_log2FC'] > 0.25)
            ].copy()

        # Use scores column for ranking if available and ranking_method is avg_log2FC
        if ranking_method == "avg_log2FC" and 'scores' in df.columns:
            df_filtered['Score'] = df.loc[df_filtered.index, 'scores']

        # Handle pct_diff ranking method when PCT columns are missing
        if ranking_method == "pct_diff" and not has_pct:
            print("Warning: pct_diff ranking requires PCT columns. Falling back to avg_log2FC.")
            ranking_method = "avg_log2FC"

        # Validate parameters and prepare ranking
        _validate_ranking_parameters(df_filtered, ranking_method, ascending)
        df_prepared, sort_column = _prepare_ranking_column(df_filtered, ranking_method)
        sort_ascending = _get_sort_direction(ranking_method, ascending)

        # Sort within each cluster by specified method and get top n genes
        top_markers = []

        for cluster in df_prepared['cluster'].unique():
            cluster_data = df_prepared[df_prepared['cluster'] == cluster]
            # Sort by specified column and direction, then take top n
            top_n = (cluster_data
                    .sort_values(sort_column, ascending=sort_ascending, na_position='last')
                    .head(n_genes))

            top_markers.append(top_n)

        # Combine all results
        if top_markers:
            top_markers = pd.concat(top_markers, ignore_index=True)

            # Create markers column by concatenating genes in order
            result = (top_markers
                     .groupby('cluster', observed=True)
                     .agg({'gene': lambda x: ','.join(x)})
                     .rename(columns={'gene': 'markers'})
                     .reset_index())

            return result
        else:
            return pd.DataFrame(columns=['cluster', 'markers'])

    elif format_type == 'scanpy':
        # Process Scanpy format
        clusters = df['names'].dtype.names
        top_markers = []

        for cluster in clusters:
            # Get data for this cluster
            genes = df['names'][cluster]
            logfc = df['logfoldchanges'][cluster].astype(float).copy()
            pvals_adj = df['pvals_adj'][cluster].astype(float).copy()
            pcts = df['pcts'][cluster].astype(float).copy()

            # Handle NaN and inf values in logfc (like Seurat format does)
            # Replace inf/-inf and NaN with max/min finite values
            finite_mask = np.isfinite(logfc)
            if finite_mask.any():
                max_finite = logfc[finite_mask].max()
                min_finite = logfc[finite_mask].min()
                # Replace inf with max, -inf with min, NaN with max (upregulated markers)
                logfc = np.where(np.isnan(logfc), max_finite, logfc)
                logfc = np.where(np.isposinf(logfc), max_finite, logfc)
                logfc = np.where(np.isneginf(logfc), min_finite, logfc)
            else:
                # If all values are non-finite, set to default
                logfc = np.ones_like(logfc) * 1.0

            # Create temporary DataFrame for sorting
            cluster_df = pd.DataFrame({
                'gene': genes,
                'avg_log2FC': logfc,
                'p_val_adj': pvals_adj,
                'pct.1': pcts,  # Assuming this represents pct.1
                'pct.2': np.zeros_like(pcts)  # May need to adjust based on data structure
            })

            # Filter for significant upregulated genes with PCT threshold
            mask = (cluster_df['p_val_adj'] < 0.05) & (cluster_df['avg_log2FC'] > 0.25) & (cluster_df['pct.1'] >= 0.1)
            filtered_df = cluster_df[mask]

            if not filtered_df.empty:
                # Validate parameters and prepare ranking
                _validate_ranking_parameters(filtered_df, ranking_method, ascending)
                df_prepared, sort_column = _prepare_ranking_column(filtered_df, ranking_method)
                sort_ascending = _get_sort_direction(ranking_method, ascending)

                # Sort and get top n genes
                top_genes_df = (df_prepared
                               .sort_values(sort_column, ascending=sort_ascending, na_position='last')
                               .head(n_genes))
                valid_genes = top_genes_df['gene'].values

                # Join genes with commas
                markers = ','.join(valid_genes)
                top_markers.append({
                    'cluster': cluster,
                    'markers': markers
                })

        return pd.DataFrame(top_markers)

    else:  # Seurat format
        # Convert string 'inf' and '-inf' to numeric values first
        df['avg_log2FC'] = pd.to_numeric(df['avg_log2FC'].replace({'inf': np.inf, '-inf': -np.inf}), errors='coerce')

        # Replace inf and -inf values with max and min non-inf values respectively
        max_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).max()
        min_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).min()
        df['avg_log2FC'] = df['avg_log2FC'].replace([np.inf, -np.inf], [max_non_inf, min_non_inf])

        # Filter by adjusted p-value, positive log2FC, and PCT
        df_filtered = df[
            (df['p_val_adj'] < 0.05) &
            (df['avg_log2FC'] > 0.25) &
            ((df['pct.1'] >= 0.1) | (df['pct.2'] >= 0.1))  # Add PCT filter
        ].copy()

        # Validate parameters and prepare ranking
        _validate_ranking_parameters(df_filtered, ranking_method, ascending)
        df_prepared, sort_column = _prepare_ranking_column(df_filtered, ranking_method)
        sort_ascending = _get_sort_direction(ranking_method, ascending)

        # Sort within each cluster by specified method and get top n genes
        top_markers = []

        for cluster in df_prepared['cluster'].unique():
            cluster_data = df_prepared[df_prepared['cluster'] == cluster]
            # Sort by specified column and direction, then take top n
            top_n = (cluster_data
                    .sort_values(sort_column, ascending=sort_ascending, na_position='last')
                    .head(n_genes))

            # Handle NaN values warning
            nan_count = cluster_data[sort_column].isna().sum()
            if nan_count > 0:
                print(f"Warning: {nan_count} NaN values found in {sort_column} column for cluster {cluster}, they will be placed at the end")

            top_markers.append(top_n)

        # Combine all results
        if top_markers:
            top_markers = pd.concat(top_markers, ignore_index=True)

            # Create markers column by concatenating genes in order
            result = (top_markers
                     .groupby('cluster', observed=True)
                     .agg({'gene': lambda x: ','.join(x)})
                     .rename(columns={'gene': 'markers'})
                     .reset_index())

            return result
        else:
            return pd.DataFrame(columns=['cluster', 'markers'])


def loadmarker(marker_type="processed"):
    """
    Load built-in marker files.

    Args:
        marker_type (str): Type of markers to load. Options:
            - "processed": For processed marker data
            - "unprocessed": For raw unprocessed marker data
            - "subcluster_results": For subcluster analysis results

    Returns:
        pandas.DataFrame: Marker data

    Raises:
        ValueError: If marker_type is not recognized
    """
    marker_files = {
        "processed": "processed.csv",
        "unprocessed": "unprocessed.csv",
        "subcluster_results": "subcluster_results.csv"
    }

    if marker_type not in marker_files:
        raise ValueError(f"Unknown marker type: {marker_type}. Available types: {list(marker_files.keys())}")

    filename = marker_files[marker_type]

    try:
        # Using importlib.resources for Python 3.7+
        with resources.path('CASSIA.data', filename) as file_path:
            return pd.read_csv(file_path)
    except Exception as e:
        raise Exception(f"Error loading marker file: {str(e)}")


def list_available_markers():
    """
    List all available built-in marker sets.

    Returns:
        list: Names of available marker sets (without .csv extension)
    """
    try:
        with resources.path('CASSIA.data', '') as data_path:
            marker_files = [f for f in os.listdir(data_path) if f.endswith('.csv')]
        return [f.replace('.csv', '') for f in marker_files]
    except Exception as e:
        raise Exception(f"Error listing marker files: {str(e)}")
