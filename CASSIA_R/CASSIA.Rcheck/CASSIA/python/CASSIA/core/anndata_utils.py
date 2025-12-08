# CASSIA - AnnData Integration Utilities
# Functions for integrating CASSIA annotation results with Scanpy AnnData objects

"""
This module provides functionality to add CASSIA cell type annotations
to Scanpy AnnData objects, similar to the R function add_cassia_to_seurat().
"""

import re
import warnings
from typing import Union, Optional, List, Dict, Any

import pandas as pd
import numpy as np

# Column name mappings (mirrors R implementation in utils.R)
COLUMN_MAPPING = {
    'general_celltype': ['Predicted General Cell Type', 'Predicted Main Cell Type',
                         'General Cell Type', 'General_Cell_Type'],
    'sub_celltype': ['Predicted Detailed Cell Type', 'Predicted Sub Cell Types',
                     'Sub Cell Type', 'Sub_Cell_Type'],
    'mixed_celltype': ['Possible Mixed Cell Types', 'Mixed Cell Types',
                       'Possible_Mixed_Cell_Types'],
    'score': ['Score', 'Consensus Score', 'Consensus_Score']
}

MERGED_GROUPING_MAPPING = {
    'merged_grouping_1': ['Merged Grouping 1', 'Merged_Grouping_1'],
    'merged_grouping_2': ['Merged Grouping 2', 'Merged_Grouping_2'],
    'merged_grouping_3': ['Merged Grouping 3', 'Merged_Grouping_3']
}

# Alternative cluster column names to try
ALT_CLUSTER_COLS = ['Cluster ID', 'True Cell Type', 'Cluster', 'cluster_id',
                    'ClusterID', 'Cluster_ID', 'cluster']


def _normalize_cluster_name(name: str) -> str:
    """
    Normalize cluster name for comparison.

    - Convert to lowercase
    - Replace underscores with spaces (for equivalence)
    - Remove punctuation
    - Collapse whitespace
    - Strip leading/trailing whitespace
    """
    name = str(name).lower()
    name = name.replace('_', ' ')        # Treat underscores as spaces
    name = re.sub(r'[^\w\s]', '', name)  # Remove punctuation
    name = re.sub(r'\s+', ' ', name)     # Collapse whitespace
    return name.strip()


def _find_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """
    Find the first matching column from a list of candidates.

    Parameters
    ----------
    df : DataFrame
        DataFrame to search in
    candidates : list
        List of possible column names

    Returns
    -------
    str or None
        First matching column name, or None if not found
    """
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _detect_cluster_col(adata) -> str:
    """
    Auto-detect cluster column in AnnData object.

    Tries 'leiden' first, then 'louvain'.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix

    Returns
    -------
    str
        Name of detected cluster column

    Raises
    ------
    ValueError
        If no cluster column found
    """
    for col in ['leiden', 'louvain']:
        if col in adata.obs.columns:
            return col

    raise ValueError(
        "Could not auto-detect cluster column. "
        "Neither 'leiden' nor 'louvain' found in adata.obs. "
        "Please specify cluster_col explicitly."
    )


def _load_cassia_results(cassia_results: Union[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Load CASSIA results from file path or DataFrame.

    Parameters
    ----------
    cassia_results : str or DataFrame
        Path to CSV file or DataFrame

    Returns
    -------
    DataFrame
        CASSIA results with normalized column names
    """
    if isinstance(cassia_results, pd.DataFrame):
        df = cassia_results.copy()
    elif isinstance(cassia_results, str):
        df = pd.read_csv(cassia_results)
    else:
        raise TypeError(
            f"cassia_results must be a file path (str) or DataFrame, "
            f"got {type(cassia_results)}"
        )

    # Fix column names (replace dots with spaces, like R implementation)
    df.columns = [col.replace('.', ' ') for col in df.columns]

    return df


def _fuzzy_match_clusters(
    adata_clusters: List[str],
    cassia_clusters: List[str]
) -> Dict[str, str]:
    """
    Create fuzzy mapping between AnnData clusters and CASSIA clusters.

    Handles variations like "0" vs "cluster_0" vs "Cluster 0".

    Parameters
    ----------
    adata_clusters : list
        Unique cluster values from AnnData
    cassia_clusters : list
        Unique cluster values from CASSIA results

    Returns
    -------
    dict
        Mapping from adata cluster names to CASSIA cluster names
    """
    # Normalize all cluster names
    adata_norm = {_normalize_cluster_name(c): c for c in adata_clusters}
    cassia_norm = {_normalize_cluster_name(c): c for c in cassia_clusters}

    mapping = {}

    for adata_cluster in adata_clusters:
        norm_adata = _normalize_cluster_name(adata_cluster)

        # Try exact match first
        if adata_cluster in cassia_clusters:
            mapping[adata_cluster] = adata_cluster
            continue

        # Try normalized exact match
        if norm_adata in cassia_norm:
            mapping[adata_cluster] = cassia_norm[norm_adata]
            continue

        # Try substring matching (e.g., "0" in "cluster_0")
        found = False
        for norm_cassia, cassia_cluster in cassia_norm.items():
            if norm_adata in norm_cassia or norm_cassia in norm_adata:
                mapping[adata_cluster] = cassia_cluster
                found = True
                break

        if not found:
            # No match found - will result in NaN for this cluster
            mapping[adata_cluster] = None

    return mapping


def _parse_sub_celltypes(value: str) -> List[Optional[str]]:
    """
    Parse comma-separated sub-celltype string into ranked list.

    Parameters
    ----------
    value : str
        Comma-separated sub-celltype string

    Returns
    -------
    list
        List of [subtype_1, subtype_2, subtype_3] (padded with None)
    """
    if pd.isna(value) or value == '':
        return [None, None, None]

    # Remove brackets and clean up
    clean = re.sub(r'^\s*\[?\s*|\s*\]?\s*$', '', str(value))

    # Split by comma and trim
    parts = [p.strip() for p in clean.split(',')]

    # Pad to 3 elements
    while len(parts) < 3:
        parts.append(None)

    return parts[:3]


def add_cassia_to_anndata(
    adata,
    cassia_results: Union[str, pd.DataFrame],
    cluster_col: Optional[str] = None,
    cassia_cluster_col: str = "Cluster ID",
    prefix: str = "CASSIA_",
    replace_existing: bool = False,
    fuzzy_match: bool = True,
    columns_to_include: int = 2,
    inplace: bool = True
):
    """
    Integrate CASSIA annotation results into AnnData object.

    This function adds CASSIA cell type annotations as columns in adata.obs,
    similar to the R function add_cassia_to_seurat().

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    cassia_results : str or DataFrame
        Path to CASSIA results CSV file OR DataFrame from runCASSIA_batch()
    cluster_col : str or None
        Column name in adata.obs containing cluster assignments.
        If None, auto-detects: tries 'leiden', then 'louvain'
    cassia_cluster_col : str
        Column name in CASSIA results for cluster IDs (default: "Cluster ID")
    prefix : str
        Prefix for new column names (default: "CASSIA_")
    replace_existing : bool
        Whether to overwrite existing CASSIA columns (default: False)
    fuzzy_match : bool
        Enable fuzzy matching for cluster ID alignment (default: True)
    columns_to_include : int
        1 = merged groupings only, 2 = all metrics (default: 2)
    inplace : bool
        If True, modify adata in place. If False, return copy (default: True)

    Returns
    -------
    AnnData or None
        Modified AnnData if inplace=False, else None

    Notes
    -----
    Adds the following columns to adata.obs (with prefix):
    - general_celltype: Main cell type prediction
    - sub_celltype: Primary sub-celltype (first from ranked list)
    - sub_celltype_all: Full comma-separated string
    - sub_celltype_1/2/3: Split ranked alternatives
    - mixed_celltype: Possible mixed cell types
    - score: Quality/consensus score (0-100)
    - merged_grouping_1/2/3: Hierarchical groupings (if available)
    - combined_celltype: Format "General :: Subtype"

    Also stores cluster-level summary in adata.uns['CASSIA']

    Examples
    --------
    >>> import scanpy as sc
    >>> from CASSIA import add_cassia_to_anndata, runCASSIA_batch
    >>>
    >>> adata = sc.read_h5ad("my_data.h5ad")
    >>>
    >>> # Option 1: From file path
    >>> add_cassia_to_anndata(adata, "cassia_results.csv")
    >>>
    >>> # Option 2: From DataFrame
    >>> results_df = runCASSIA_batch(markers_df, tissue="lung", species="human")
    >>> add_cassia_to_anndata(adata, results_df)
    >>>
    >>> # Access annotations
    >>> print(adata.obs['CASSIA_general_celltype'])
    >>> print(adata.uns['CASSIA'])  # Cluster-level summary
    """
    # Import anndata here to make it an optional dependency
    try:
        import anndata
    except ImportError:
        raise ImportError(
            "The anndata package is required for this function. "
            "Install it with: pip install anndata"
        )

    # Validate columns_to_include
    if columns_to_include not in [1, 2]:
        warnings.warn(
            f"Invalid columns_to_include={columns_to_include}. Using default 2 (all metrics)."
        )
        columns_to_include = 2

    # Handle inplace
    if not inplace:
        adata = adata.copy()

    # Auto-detect cluster column if not specified
    if cluster_col is None:
        cluster_col = _detect_cluster_col(adata)
        print(f"Auto-detected cluster column: '{cluster_col}'")

    # Validate cluster column exists
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"Cluster column '{cluster_col}' not found in adata.obs")

    # Load CASSIA results
    cassia_df = _load_cassia_results(cassia_results)

    # Fix cassia_cluster_col if it has dots
    cassia_cluster_col = cassia_cluster_col.replace('.', ' ')

    # Find cluster column in CASSIA results
    if cassia_cluster_col not in cassia_df.columns:
        found_col = _find_column(cassia_df, ALT_CLUSTER_COLS)
        if found_col is None:
            raise ValueError(
                f"Could not find cluster column in CASSIA results. "
                f"Available columns: {list(cassia_df.columns)}"
            )
        cassia_cluster_col = found_col
        print(f"Using '{cassia_cluster_col}' as cluster column from CASSIA results")

    # Get unique clusters
    adata_clusters = adata.obs[cluster_col].astype(str).unique().tolist()
    cassia_clusters = cassia_df[cassia_cluster_col].astype(str).unique().tolist()

    # Build cluster mapping
    if fuzzy_match:
        cluster_map = _fuzzy_match_clusters(adata_clusters, cassia_clusters)

        # Report any mappings created
        mapped = {k: v for k, v in cluster_map.items() if v is not None and k != v}
        if mapped:
            print("Created cluster mapping:")
            for k, v in mapped.items():
                print(f"  '{k}' -> '{v}'")

        # Warn about unmatched
        unmatched = [k for k, v in cluster_map.items() if v is None]
        if unmatched:
            warnings.warn(f"Could not find matches for clusters: {unmatched}")
    else:
        cluster_map = {c: c if c in cassia_clusters else None for c in adata_clusters}

    # Create lookup dictionaries from CASSIA results
    cassia_df_indexed = cassia_df.set_index(cassia_df[cassia_cluster_col].astype(str))

    # Determine which columns to process
    all_mappings = {**COLUMN_MAPPING, **MERGED_GROUPING_MAPPING}

    if columns_to_include == 1:
        # Only merged groupings
        available_merged = {}
        for key, candidates in MERGED_GROUPING_MAPPING.items():
            col = _find_column(cassia_df, candidates)
            if col:
                available_merged[key] = col

        if not available_merged:
            warnings.warn(
                "No merged grouping columns found. Using all metrics instead."
            )
            columns_to_include = 2
        else:
            columns_to_process = available_merged

    if columns_to_include == 2:
        # All metrics
        columns_to_process = {}
        for key, candidates in all_mappings.items():
            col = _find_column(cassia_df, candidates)
            if col:
                columns_to_process[key] = col

    # Map clusters for each cell
    cell_clusters = adata.obs[cluster_col].astype(str)
    mapped_clusters = cell_clusters.map(lambda x: cluster_map.get(x))

    # Check for existing columns
    existing_cols = [f"{prefix}{k}" for k in columns_to_process.keys()
                     if f"{prefix}{k}" in adata.obs.columns]
    if existing_cols and not replace_existing:
        raise ValueError(
            f"CASSIA columns already exist: {existing_cols}. "
            f"Set replace_existing=True to overwrite."
        )

    # Add annotation columns to adata.obs
    for key, source_col in columns_to_process.items():
        if key == 'sub_celltype':
            # Handle sub_celltype separately (needs splitting)
            continue

        col_name = f"{prefix}{key}"

        # Create lookup from CASSIA cluster to value
        lookup = cassia_df_indexed[source_col].to_dict()

        # Map values for each cell
        adata.obs[col_name] = mapped_clusters.map(lookup)

    # Special handling for sub_celltype (split into multiple columns)
    if columns_to_include == 2:
        sub_col = _find_column(cassia_df, COLUMN_MAPPING['sub_celltype'])
        if sub_col:
            # Store full string
            lookup_all = cassia_df_indexed[sub_col].to_dict()
            adata.obs[f"{prefix}sub_celltype_all"] = mapped_clusters.map(lookup_all)

            # Parse and split
            parsed = cassia_df_indexed[sub_col].apply(_parse_sub_celltypes)

            for i in range(3):
                lookup_i = {k: v[i] for k, v in parsed.items()}
                adata.obs[f"{prefix}sub_celltype_{i+1}"] = mapped_clusters.map(lookup_i)

            # Main sub_celltype is the first one
            adata.obs[f"{prefix}sub_celltype"] = adata.obs[f"{prefix}sub_celltype_1"]

    # Add combined celltype column (General :: Subtype)
    if columns_to_include == 2:
        general_col = f"{prefix}general_celltype"
        sub_col = f"{prefix}sub_celltype_1"
        combined_col = f"{prefix}combined_celltype"

        if general_col in adata.obs.columns:
            if sub_col in adata.obs.columns:
                adata.obs[combined_col] = adata.obs.apply(
                    lambda row: (
                        f"{row[general_col]} :: {row[sub_col]}"
                        if pd.notna(row[sub_col]) and row[sub_col] != ''
                        else row[general_col]
                    ),
                    axis=1
                )
            else:
                adata.obs[combined_col] = adata.obs[general_col]

            print(f"Added combined cell type column: {combined_col}")

    # Store cluster-level summary in adata.uns
    # Select relevant columns for summary
    summary_cols = [cassia_cluster_col]
    for key, candidates in all_mappings.items():
        col = _find_column(cassia_df, candidates)
        if col and col not in summary_cols:
            summary_cols.append(col)

    adata.uns['CASSIA'] = cassia_df[summary_cols].copy()
    print(f"Stored cluster-level summary in adata.uns['CASSIA']")

    # Report summary
    added_cols = [c for c in adata.obs.columns if c.startswith(prefix)]
    print(f"Added {len(added_cols)} columns to adata.obs: {added_cols}")

    if not inplace:
        return adata


def _calculate_pct_expressing(
    adata,
    gene: str,
    cluster_col: str,
    cluster_id,
    min_expression: float = 0.0
) -> tuple:
    """
    Calculate percentage of cells expressing a gene in and outside a cluster.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    gene : str
        Gene name to check expression for
    cluster_col : str
        Column name in adata.obs containing cluster assignments
    cluster_id
        Cluster ID to calculate percentages for
    min_expression : float
        Threshold above which a cell is considered "expressing" the gene

    Returns
    -------
    tuple
        (pct1, pct2) where:
        - pct1: fraction of cells in cluster expressing gene
        - pct2: fraction of cells outside cluster expressing gene
    """
    try:
        gene_idx = list(adata.var_names).index(gene)
    except ValueError:
        # Gene not found
        return np.nan, np.nan

    # Get expression values for this gene
    expr = adata.X[:, gene_idx]

    # Handle sparse matrices
    if hasattr(expr, 'toarray'):
        expr = expr.toarray().flatten()
    elif hasattr(expr, 'A'):
        expr = expr.A.flatten()
    else:
        expr = np.asarray(expr).flatten()

    # Get cluster mask
    is_cluster = (adata.obs[cluster_col].astype(str) == str(cluster_id)).values

    # Calculate percentages
    n_in_cluster = is_cluster.sum()
    n_outside_cluster = (~is_cluster).sum()

    if n_in_cluster > 0:
        pct1 = (expr[is_cluster] > min_expression).sum() / n_in_cluster
    else:
        pct1 = np.nan

    if n_outside_cluster > 0:
        pct2 = (expr[~is_cluster] > min_expression).sum() / n_outside_cluster
    else:
        pct2 = np.nan

    return float(pct1), float(pct2)


def enhance_scanpy_markers(
    adata,
    cluster_col: Optional[str] = None,
    n_genes: Optional[int] = None,
    min_expression: float = 0.0,
    include_stats: bool = True,
    key: str = "rank_genes_groups"
) -> pd.DataFrame:
    """
    Extract and enhance Scanpy marker genes with pct.1 and pct.2 values.

    This function extracts marker genes from Scanpy's rank_genes_groups results
    and calculates the percentage of cells expressing each gene within and
    outside each cluster. This format is compatible with CASSIA's annotation
    boost agent.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with rank_genes_groups results in .uns
    cluster_col : str or None
        Column name in adata.obs containing cluster assignments.
        If None, auto-detects: tries 'leiden', then 'louvain'
    n_genes : int or None
        Number of top genes per cluster to include. If None, includes all genes.
    min_expression : float
        Threshold above which a cell is considered "expressing" a gene.
        Default is 0.0 (any expression counts).
    include_stats : bool
        Whether to include additional statistics (logfoldchanges, pvals, scores).
        Default is True.
    key : str
        Key in adata.uns where rank_genes_groups results are stored.
        Default is "rank_genes_groups".

    Returns
    -------
    DataFrame
        DataFrame with columns:
        - cluster: Cluster ID
        - gene: Gene name
        - pct.1: Fraction of cells in cluster expressing gene (0.0-1.0)
        - pct.2: Fraction of cells outside cluster expressing gene (0.0-1.0)
        - avg_log2FC: Log fold change (if include_stats=True and available)
        - p_val_adj: Adjusted p-value (if include_stats=True and available)
        - scores: Scanpy score (if include_stats=True and available)

    Raises
    ------
    ValueError
        If rank_genes_groups results not found in adata.uns
    ImportError
        If anndata is not installed

    Examples
    --------
    >>> import scanpy as sc
    >>> from CASSIA import enhance_scanpy_markers, runCASSIA_batch
    >>>
    >>> # Standard Scanpy workflow
    >>> adata = sc.read_h5ad("pbmc3k.h5ad")
    >>> sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    >>>
    >>> # Enhance markers with pct values
    >>> markers_df = enhance_scanpy_markers(adata, n_genes=50)
    >>>
    >>> # Now compatible with annotation boost
    >>> results = runCASSIA_batch(markers_df, tissue="blood", species="human")
    """
    # Import anndata here to make it an optional dependency
    try:
        import anndata
    except ImportError:
        raise ImportError(
            "The anndata package is required for this function. "
            "Install it with: pip install anndata"
        )

    # Check for rank_genes_groups results
    if key not in adata.uns:
        raise ValueError(
            f"rank_genes_groups results not found in adata.uns['{key}']. "
            f"Please run sc.tl.rank_genes_groups() first."
        )

    rgg = adata.uns[key]

    # Auto-detect cluster column if not specified
    if cluster_col is None:
        cluster_col = _detect_cluster_col(adata)
        print(f"Auto-detected cluster column: '{cluster_col}'")

    # Get cluster names from rank_genes_groups results
    if 'names' not in rgg:
        raise ValueError("rank_genes_groups results missing 'names' field")

    # Get clusters - handle both structured array and dict formats
    if hasattr(rgg['names'], 'dtype') and rgg['names'].dtype.names:
        # Structured array format
        clusters = list(rgg['names'].dtype.names)
    else:
        # Try to get from params or infer from data
        if 'params' in rgg and 'groupby' in rgg['params']:
            clusters = adata.obs[rgg['params']['groupby']].unique().tolist()
        else:
            clusters = adata.obs[cluster_col].unique().tolist()

    # Determine number of genes per cluster
    if hasattr(rgg['names'], 'dtype') and rgg['names'].dtype.names:
        max_genes = len(rgg['names'])
    else:
        # Assume it's a dict-like structure
        first_cluster = clusters[0]
        max_genes = len(rgg['names'][first_cluster]) if isinstance(rgg['names'], dict) else len(rgg['names'])

    if n_genes is None:
        n_genes = max_genes
    else:
        n_genes = min(n_genes, max_genes)

    # Build results
    results = []

    for cluster in clusters:
        # Extract gene names for this cluster
        if hasattr(rgg['names'], 'dtype') and rgg['names'].dtype.names:
            # Structured array format
            genes = rgg['names'][cluster][:n_genes]
        elif isinstance(rgg['names'], dict):
            genes = rgg['names'][cluster][:n_genes]
        else:
            # Assume it's a 2D array with cluster as column index
            cluster_idx = clusters.index(cluster)
            genes = rgg['names'][:n_genes, cluster_idx] if rgg['names'].ndim > 1 else rgg['names'][:n_genes]

        for i, gene in enumerate(genes):
            if pd.isna(gene) or gene == '':
                continue

            gene = str(gene)

            # Calculate pct.1 and pct.2
            pct1, pct2 = _calculate_pct_expressing(
                adata, gene, cluster_col, cluster, min_expression
            )

            row = {
                'cluster': str(cluster),
                'gene': gene,
                'pct.1': pct1,
                'pct.2': pct2
            }

            # Add additional stats if requested
            if include_stats:
                # Log fold change
                if 'logfoldchanges' in rgg:
                    try:
                        if hasattr(rgg['logfoldchanges'], 'dtype') and rgg['logfoldchanges'].dtype.names:
                            row['avg_log2FC'] = float(rgg['logfoldchanges'][cluster][i])
                        elif isinstance(rgg['logfoldchanges'], dict):
                            row['avg_log2FC'] = float(rgg['logfoldchanges'][cluster][i])
                        else:
                            cluster_idx = clusters.index(cluster)
                            row['avg_log2FC'] = float(rgg['logfoldchanges'][i, cluster_idx])
                    except (IndexError, KeyError, TypeError):
                        row['avg_log2FC'] = np.nan

                # Adjusted p-value
                if 'pvals_adj' in rgg:
                    try:
                        if hasattr(rgg['pvals_adj'], 'dtype') and rgg['pvals_adj'].dtype.names:
                            row['p_val_adj'] = float(rgg['pvals_adj'][cluster][i])
                        elif isinstance(rgg['pvals_adj'], dict):
                            row['p_val_adj'] = float(rgg['pvals_adj'][cluster][i])
                        else:
                            cluster_idx = clusters.index(cluster)
                            row['p_val_adj'] = float(rgg['pvals_adj'][i, cluster_idx])
                    except (IndexError, KeyError, TypeError):
                        row['p_val_adj'] = np.nan

                # Scores
                if 'scores' in rgg:
                    try:
                        if hasattr(rgg['scores'], 'dtype') and rgg['scores'].dtype.names:
                            row['scores'] = float(rgg['scores'][cluster][i])
                        elif isinstance(rgg['scores'], dict):
                            row['scores'] = float(rgg['scores'][cluster][i])
                        else:
                            cluster_idx = clusters.index(cluster)
                            row['scores'] = float(rgg['scores'][i, cluster_idx])
                    except (IndexError, KeyError, TypeError):
                        row['scores'] = np.nan

            results.append(row)

    df = pd.DataFrame(results)

    # Ensure proper column order
    col_order = ['cluster', 'gene', 'pct.1', 'pct.2']
    if include_stats:
        for col in ['avg_log2FC', 'p_val_adj', 'scores']:
            if col in df.columns:
                col_order.append(col)

    df = df[[c for c in col_order if c in df.columns]]

    print(f"Enhanced {len(df)} marker genes across {len(clusters)} clusters with pct.1/pct.2 values")

    return df
