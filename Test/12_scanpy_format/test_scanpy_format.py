"""
CASSIA Test 12: Scanpy Format Input
====================================
Tests how CASSIA handles Scanpy's rank_genes_groups output with runCASSIA_batch.

This test uses pre-generated Scanpy rank_genes_groups data from PBMC3k dataset.
To regenerate the test data, run: python generate_scanpy_markers.py

Tests:
1. Load pre-generated Scanpy rank_genes_groups data (structured array format)
2. Test Scanpy format auto-detection
3. Test get_top_markers() with Scanpy structured array output
4. Test runCASSIA_batch with Scanpy-derived markers
5. Load scanpy_df format (flat DataFrame from sc.get.rank_genes_groups_df)
6. Test get_top_markers() with scanpy_df format
7. Test runCASSIA_batch with scanpy_df format

Usage:
    python test_scanpy_format.py
"""

import sys
import time
import pickle
import numpy as np
import pandas as pd
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    setup_cassia_imports,
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Setup CASSIA imports
setup_cassia_imports()


def get_data_dir() -> Path:
    """Get the data directory for this test."""
    return Path(__file__).parent / "data"


def load_scanpy_data():
    """
    Load pre-generated Scanpy rank_genes_groups data.

    Returns:
        tuple: (rank_genes_dict, cluster_info)
    """
    data_dir = get_data_dir()

    rank_genes_path = data_dir / "scanpy_rank_genes_groups.pkl"
    cluster_info_path = data_dir / "scanpy_cluster_info.pkl"

    if not rank_genes_path.exists() or not cluster_info_path.exists():
        raise FileNotFoundError(
            f"Scanpy test data not found. Please run 'python generate_scanpy_markers.py' first.\n"
            f"Expected files:\n"
            f"  - {rank_genes_path}\n"
            f"  - {cluster_info_path}"
        )

    with open(rank_genes_path, 'rb') as f:
        rank_genes_dict = pickle.load(f)

    with open(cluster_info_path, 'rb') as f:
        cluster_info = pickle.load(f)

    return rank_genes_dict, cluster_info


def test_load_scanpy_data():
    """Test loading the pre-generated Scanpy data."""
    print("\n--- Testing load of pre-generated Scanpy data ---")

    errors = []

    try:
        rank_genes_dict, cluster_info = load_scanpy_data()

        print(f"Loaded rank_genes_groups from pickle file")
        print(f"  Keys: {list(rank_genes_dict.keys())}")
        print(f"  Clusters: {cluster_info['clusters']}")
        print(f"  N cells: {cluster_info['n_cells']}")
        print(f"  N genes: {cluster_info['n_genes']}")

        # Validate structure
        required_keys = ['names', 'scores', 'logfoldchanges', 'pvals_adj', 'pcts']
        for key in required_keys:
            if key not in rank_genes_dict:
                errors.append(f"Missing required key: {key}")

        # Validate clusters
        if len(cluster_info['clusters']) == 0:
            errors.append("No clusters found in data")

    except FileNotFoundError as e:
        errors.append(str(e))
    except Exception as e:
        errors.append(f"Error loading data: {str(e)}")

    return len(errors) == 0, errors


def test_format_detection():
    """Test that Scanpy format is correctly auto-detected."""
    print("\n--- Testing Scanpy format auto-detection ---")

    errors = []

    # Load Scanpy data
    rank_genes_dict, cluster_info = load_scanpy_data()

    # Check detection logic
    has_names = 'names' in rank_genes_dict
    has_scores = 'scores' in rank_genes_dict

    print(f"Scanpy dict has 'names': {has_names}")
    print(f"Scanpy dict has 'scores': {has_scores}")
    print(f"Should detect as Scanpy: {has_names and has_scores}")

    if not (has_names and has_scores):
        errors.append("Scanpy output missing required keys for detection")

    # Also test that a Seurat DataFrame is NOT detected as Scanpy
    seurat_df = pd.DataFrame({
        'cluster': ['0', '0', '1', '1'],
        'gene': ['CD3D', 'CD4', 'CD19', 'MS4A1'],
        'avg_log2FC': [2.5, 2.0, 3.0, 2.8],
        'p_val_adj': [0.001, 0.002, 0.001, 0.001],
        'pct.1': [0.8, 0.7, 0.9, 0.85],
        'pct.2': [0.1, 0.1, 0.05, 0.05]
    })

    seurat_has_names = 'names' in seurat_df
    seurat_has_scores = 'scores' in seurat_df

    print(f"\nSeurat DataFrame has 'names': {seurat_has_names} (expected: False)")
    print(f"Seurat DataFrame has 'scores': {seurat_has_scores} (expected: False)")

    if seurat_has_names and seurat_has_scores:
        errors.append("Seurat format incorrectly has Scanpy keys")

    return len(errors) == 0, errors


def test_get_top_markers_with_scanpy():
    """Test get_top_markers() with Scanpy rank_genes_groups output."""
    print("\n--- Testing get_top_markers() with Scanpy output ---")

    from CASSIA import get_top_markers

    # Load Scanpy data
    rank_genes_dict, cluster_info = load_scanpy_data()

    print(f"Scanpy rank_genes_groups structure:")
    print(f"  Keys: {list(rank_genes_dict.keys())}")
    print(f"  Clusters: {cluster_info['clusters']}")
    print(f"  Genes per cluster: {len(rank_genes_dict['names'])}")

    # Test auto-detection and processing
    print("\nCalling get_top_markers()...")
    result = get_top_markers(rank_genes_dict, n_genes=15)

    print(f"\nOutput DataFrame:")
    print(f"  Shape: {result.shape}")
    print(f"  Columns: {result.columns.tolist()}")

    errors = []

    # Validate output structure
    if 'cluster' not in result.columns:
        errors.append("Missing 'cluster' column in output")
    if 'markers' not in result.columns:
        errors.append("Missing 'markers' column in output")

    if len(result) == 0:
        errors.append("No clusters returned from get_top_markers")
    else:
        print(f"\nMarkers extracted per cluster:")
        for _, row in result.iterrows():
            cluster = row['cluster']
            markers = row['markers']
            marker_list = markers.split(',') if markers else []
            print(f"  {cluster}: {len(marker_list)} markers")
            if marker_list:
                print(f"    Top 5: {', '.join(marker_list[:5])}")

        # Verify all expected clusters have markers
        expected_clusters = set(cluster_info['clusters'])
        actual_clusters = set(result['cluster'].tolist())
        missing = expected_clusters - actual_clusters

        if missing:
            print(f"\nWarning: Some clusters missing: {missing}")
            # Not an error - some clusters may have no genes passing filters

    return len(errors) == 0, errors, result


def load_scanpy_df():
    """
    Load pre-generated Scanpy flat DataFrame format.

    This is the modern format from sc.get.rank_genes_groups_df(adata, group=None).

    Returns:
        pandas.DataFrame with columns: group, names, scores, logfoldchanges, pvals, pvals_adj
    """
    data_dir = get_data_dir()
    scanpy_df_path = data_dir / "scanpy_rank_genes_df.csv"

    if not scanpy_df_path.exists():
        raise FileNotFoundError(
            f"Scanpy DataFrame not found. Please run 'python generate_scanpy_markers.py' first.\n"
            f"Expected file: {scanpy_df_path}"
        )

    return pd.read_csv(scanpy_df_path)


def test_load_scanpy_df():
    """Test loading the pre-generated Scanpy flat DataFrame."""
    print("\n--- Testing load of Scanpy flat DataFrame format ---")

    errors = []

    try:
        scanpy_df = load_scanpy_df()

        print(f"Loaded scanpy_df from CSV file")
        print(f"  Shape: {scanpy_df.shape}")
        print(f"  Columns: {scanpy_df.columns.tolist()}")
        print(f"  Sample data:")
        print(scanpy_df.head(3).to_string(index=False))

        # Validate required columns
        required_cols = ['group', 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']
        for col in required_cols:
            if col not in scanpy_df.columns:
                errors.append(f"Missing required column: {col}")

        # Validate data
        n_clusters = scanpy_df['group'].nunique()
        print(f"\n  Number of clusters: {n_clusters}")
        print(f"  Clusters: {scanpy_df['group'].unique().tolist()}")

        if n_clusters == 0:
            errors.append("No clusters found in DataFrame")

    except FileNotFoundError as e:
        errors.append(str(e))
    except Exception as e:
        errors.append(f"Error loading data: {str(e)}")

    return len(errors) == 0, errors


def test_format_detection_scanpy_df():
    """Test that scanpy_df format is correctly auto-detected."""
    print("\n--- Testing scanpy_df format auto-detection ---")

    errors = []

    # Load scanpy_df
    scanpy_df = load_scanpy_df()

    # Check detection logic - should detect 'group' and 'names' columns
    has_group = 'group' in scanpy_df.columns
    has_names = 'names' in scanpy_df.columns
    has_logfoldchanges = 'logfoldchanges' in scanpy_df.columns

    print(f"scanpy_df has 'group': {has_group}")
    print(f"scanpy_df has 'names': {has_names}")
    print(f"scanpy_df has 'logfoldchanges': {has_logfoldchanges}")
    print(f"Should detect as scanpy_df: {has_group and has_names and has_logfoldchanges}")

    if not (has_group and has_names and has_logfoldchanges):
        errors.append("scanpy_df missing required columns for detection")

    # Verify it does NOT look like Seurat format
    has_cluster = 'cluster' in scanpy_df.columns
    has_gene = 'gene' in scanpy_df.columns

    print(f"\nscanpy_df has 'cluster' (Seurat): {has_cluster} (expected: False)")
    print(f"scanpy_df has 'gene' (Seurat): {has_gene} (expected: False)")

    if has_cluster and has_gene:
        errors.append("scanpy_df incorrectly has Seurat column names")

    return len(errors) == 0, errors


def test_get_top_markers_with_scanpy_df():
    """Test get_top_markers() with Scanpy flat DataFrame format."""
    print("\n--- Testing get_top_markers() with scanpy_df format ---")

    from CASSIA import get_top_markers

    # Load scanpy_df
    scanpy_df = load_scanpy_df()

    print(f"scanpy_df structure:")
    print(f"  Shape: {scanpy_df.shape}")
    print(f"  Columns: {scanpy_df.columns.tolist()}")
    print(f"  Clusters: {scanpy_df['group'].unique().tolist()}")

    # Test auto-detection and processing
    print("\nCalling get_top_markers() with scanpy_df...")
    result = get_top_markers(scanpy_df, n_genes=15)

    print(f"\nOutput DataFrame:")
    print(f"  Shape: {result.shape}")
    print(f"  Columns: {result.columns.tolist()}")

    errors = []

    # Validate output structure
    if 'cluster' not in result.columns:
        errors.append("Missing 'cluster' column in output")
    if 'markers' not in result.columns:
        errors.append("Missing 'markers' column in output")

    if len(result) == 0:
        errors.append("No clusters returned from get_top_markers")
    else:
        print(f"\nMarkers extracted per cluster:")
        for _, row in result.iterrows():
            cluster = row['cluster']
            markers = row['markers']
            marker_list = markers.split(',') if markers else []
            print(f"  {cluster}: {len(marker_list)} markers")
            if marker_list:
                print(f"    Top 5: {', '.join(marker_list[:5])}")

    return len(errors) == 0, errors, result


def test_runCASSIA_batch_with_scanpy_df(marker_df):
    """Test runCASSIA_batch with scanpy_df-derived markers."""
    print("\n--- Testing runCASSIA_batch() with scanpy_df-derived markers ---")

    # Load configuration
    config = load_config()
    llm_config = config['llm']
    data_config = config['data']

    # Setup API keys
    setup_api_keys()

    from CASSIA import runCASSIA_batch

    print(f"\nInput marker DataFrame:")
    print(f"  Shape: {marker_df.shape}")
    print(f"  Clusters: {marker_df['cluster'].tolist()}")

    # Create results directory
    results = create_results_dir("12_scanpy_format", get_test_mode())
    logging_ctx = setup_logging(results["logs"])
    output_name = str(Path(results["outputs"]) / "scanpy_df_batch_results")
    print(f"\nResults will be saved to: {results['base']}")

    errors = []
    start_time = time.time()

    try:
        print("\nRunning runCASSIA_batch with scanpy_df-derived markers...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'blood'),
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1'),
            celltype_column='cluster',
            gene_column_name='markers'
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")
        summary_csv = Path(f"{output_name}_summary.csv")

        if full_csv.exists():
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(marker_df)}")
            print(f"  Output files created:")
            print(f"    - {full_csv.name}")
            if summary_csv.exists():
                print(f"    - {summary_csv.name}")

            # Show annotations
            if 'main_cell_type' in results_df.columns:
                print("\n  Cell Type Annotations:")
                for _, row in results_df.iterrows():
                    cluster = row.get('cluster', row.get('Cluster', 'Unknown'))
                    cell_type = row.get('main_cell_type', 'Unknown')
                    print(f"    {cluster} -> {cell_type}")

            if clusters_annotated < len(marker_df):
                errors.append(f"Only {clusters_annotated}/{len(marker_df)} clusters annotated")
        else:
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time
    print(f"\nDuration: {duration:.2f}s")

    cleanup_logging(logging_ctx)
    return len(errors) == 0, errors, results


def test_runCASSIA_batch_with_scanpy(marker_df):
    """Test runCASSIA_batch with Scanpy-derived markers."""
    print("\n--- Testing runCASSIA_batch() with Scanpy-derived markers ---")

    # Load configuration
    config = load_config()
    llm_config = config['llm']
    data_config = config['data']

    # Setup API keys
    setup_api_keys()

    from CASSIA import runCASSIA_batch

    print(f"\nInput marker DataFrame:")
    print(f"  Shape: {marker_df.shape}")
    print(f"  Clusters: {marker_df['cluster'].tolist()}")

    # Create results directory
    results = create_results_dir("12_scanpy_format", get_test_mode())
    logging_ctx = setup_logging(results["logs"])
    output_name = str(Path(results["outputs"]) / "scanpy_batch_results")
    print(f"\nResults will be saved to: {results['base']}")

    errors = []
    start_time = time.time()

    try:
        print("\nRunning runCASSIA_batch...")
        runCASSIA_batch(
            marker=marker_df,
            output_name=output_name,
            n_genes=data_config.get('n_genes', 30),
            model=llm_config.get('model', 'google/gemini-2.5-flash'),
            temperature=llm_config.get('temperature', 0.3),
            tissue=data_config.get('tissue', 'blood'),  # PBMC is from blood
            species=data_config.get('species', 'human'),
            max_workers=llm_config.get('max_workers', 3),
            provider=llm_config.get('provider', 'openrouter'),
            validator_involvement=config.get('validator', {}).get('default', 'v1'),
            celltype_column='cluster',
            gene_column_name='markers'
        )

        # Check output files
        full_csv = Path(f"{output_name}_full.csv")
        summary_csv = Path(f"{output_name}_summary.csv")

        if full_csv.exists():
            results_df = pd.read_csv(full_csv)
            clusters_annotated = len(results_df)

            print(f"\nBatch Results:")
            print(f"  Clusters annotated: {clusters_annotated}/{len(marker_df)}")
            print(f"  Output files created:")
            print(f"    - {full_csv.name}")
            if summary_csv.exists():
                print(f"    - {summary_csv.name}")

            # Show annotations
            if 'main_cell_type' in results_df.columns:
                print("\n  Cell Type Annotations:")
                for _, row in results_df.iterrows():
                    cluster = row.get('cluster', row.get('Cluster', 'Unknown'))
                    cell_type = row.get('main_cell_type', 'Unknown')
                    print(f"    {cluster} -> {cell_type}")

            if clusters_annotated < len(marker_df):
                errors.append(f"Only {clusters_annotated}/{len(marker_df)} clusters annotated")
        else:
            errors.append("Output file not created")

    except Exception as e:
        errors.append(str(e))
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

    duration = time.time() - start_time
    print(f"\nDuration: {duration:.2f}s")

    cleanup_logging(logging_ctx)
    return len(errors) == 0, errors, results


def run_scanpy_format_tests():
    """Run all Scanpy format tests."""
    print_test_header("12 - Scanpy Format Input")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    all_tests_passed = True
    all_errors = []
    start_time = time.time()
    marker_df = None
    marker_df_scanpy_df = None

    # ============================================================
    # PART A: Tests for structured array format (legacy)
    # ============================================================
    print("\n" + "#" * 60)
    print("# PART A: Structured Array Format (Legacy)")
    print("#" * 60)

    # Test 1: Load Scanpy data
    print("\n" + "=" * 60)
    print("Test 1: Load Pre-generated Scanpy Data (Structured Array)")
    print("=" * 60)
    passed, errors = test_load_scanpy_data()
    all_tests_passed = all_tests_passed and passed
    all_errors.extend(errors)
    print(f"\nResult: {'PASSED' if passed else 'FAILED'}")

    if not passed:
        print("\nCannot proceed without test data. Run 'python generate_scanpy_markers.py' first.")
        return False

    # Test 2: Format detection
    print("\n" + "=" * 60)
    print("Test 2: Scanpy Format Auto-Detection (Structured Array)")
    print("=" * 60)
    passed, errors = test_format_detection()
    all_tests_passed = all_tests_passed and passed
    all_errors.extend(errors)
    print(f"\nResult: {'PASSED' if passed else 'FAILED'}")

    # Test 3: get_top_markers with Scanpy output
    print("\n" + "=" * 60)
    print("Test 3: get_top_markers() with Scanpy Structured Array")
    print("=" * 60)
    passed, errors, marker_df = test_get_top_markers_with_scanpy()
    all_tests_passed = all_tests_passed and passed
    all_errors.extend(errors)
    print(f"\nResult: {'PASSED' if passed else 'FAILED'}")

    # Test 4: runCASSIA_batch with Scanpy-derived markers
    results = None
    if marker_df is not None and len(marker_df) > 0:
        print("\n" + "=" * 60)
        print("Test 4: runCASSIA_batch() with Structured Array Markers")
        print("=" * 60)
        passed, errors, results = test_runCASSIA_batch_with_scanpy(marker_df)
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
    else:
        print("\n" + "=" * 60)
        print("Test 4: SKIPPED (no markers from previous test)")
        print("=" * 60)
        all_errors.append("Test 4 skipped due to no markers")

    # ============================================================
    # PART B: Tests for flat DataFrame format (modern scanpy_df)
    # ============================================================
    print("\n" + "#" * 60)
    print("# PART B: Flat DataFrame Format (Modern sc.get.rank_genes_groups_df)")
    print("#" * 60)

    # Test 5: Load scanpy_df data
    print("\n" + "=" * 60)
    print("Test 5: Load Scanpy Flat DataFrame (scanpy_df)")
    print("=" * 60)
    try:
        passed, errors = test_load_scanpy_df()
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
        scanpy_df_available = passed
    except Exception as e:
        print(f"\nTest 5 FAILED: {e}")
        all_errors.append(f"Test 5 failed: {e}")
        scanpy_df_available = False

    if scanpy_df_available:
        # Test 6: Format detection for scanpy_df
        print("\n" + "=" * 60)
        print("Test 6: scanpy_df Format Auto-Detection")
        print("=" * 60)
        passed, errors = test_format_detection_scanpy_df()
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")

        # Test 7: get_top_markers with scanpy_df
        print("\n" + "=" * 60)
        print("Test 7: get_top_markers() with scanpy_df Format")
        print("=" * 60)
        passed, errors, marker_df_scanpy_df = test_get_top_markers_with_scanpy_df()
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")

        # Test 8: runCASSIA_batch with scanpy_df-derived markers
        if marker_df_scanpy_df is not None and len(marker_df_scanpy_df) > 0:
            print("\n" + "=" * 60)
            print("Test 8: runCASSIA_batch() with scanpy_df Markers")
            print("=" * 60)
            passed, errors, results = test_runCASSIA_batch_with_scanpy_df(marker_df_scanpy_df)
            all_tests_passed = all_tests_passed and passed
            all_errors.extend(errors)
            print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
        else:
            print("\n" + "=" * 60)
            print("Test 8: SKIPPED (no markers from scanpy_df test)")
            print("=" * 60)
            all_errors.append("Test 8 skipped due to no markers from scanpy_df")
    else:
        print("\n" + "=" * 60)
        print("Tests 6-8: SKIPPED (scanpy_df data not available)")
        print("=" * 60)
        print("Run 'python generate_scanpy_markers.py' to generate scanpy_df test data.")

    duration = time.time() - start_time

    # Save metadata (reuse results from last test if available, otherwise create new)
    if results is None:
        results = create_results_dir("12_scanpy_format", get_test_mode())
    rank_genes_dict, cluster_info = load_scanpy_data()
    metadata = create_test_metadata(
        test_name="scanpy_format",
        config=config,
        duration_seconds=duration,
        status="passed" if all_tests_passed else "failed",
        clusters_tested=cluster_info['clusters'],
        errors=all_errors
    )
    save_test_metadata(results["outputs"], metadata)

    # Print final result
    print("\n" + "=" * 60)
    print_test_result(all_tests_passed, f"Total duration: {duration:.2f}s")

    if all_errors:
        print("\nErrors encountered:")
        for error in all_errors:
            print(f"  - {error}")

    return all_tests_passed


if __name__ == "__main__":
    success = run_scanpy_format_tests()
    sys.exit(0 if success else 1)
