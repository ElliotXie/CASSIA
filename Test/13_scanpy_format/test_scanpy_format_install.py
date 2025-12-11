"""
CASSIA Test 13: Scanpy Format Input (PIP INSTALL MODE)
=======================================================
Tests how CASSIA handles Scanpy's rank_genes_groups output with runCASSIA_batch.
Uses pip-installed CASSIA package.

This test uses pre-generated Scanpy rank_genes_groups data from PBMC3k dataset.
To regenerate the test data, run: python generate_scanpy_markers.py

Usage:
    python test_scanpy_format_install.py
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
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    print_config_summary,
    verify_cassia_pip_install,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)

# Import CASSIA directly from pip-installed package (no setup_cassia_imports)
import CASSIA


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


def load_scanpy_df():
    """
    Load pre-generated Scanpy flat DataFrame format.

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


def test_get_top_markers_with_scanpy():
    """Test get_top_markers() with Scanpy rank_genes_groups output."""
    print("\n--- Testing get_top_markers() with Scanpy output ---")

    # Load Scanpy data
    rank_genes_dict, cluster_info = load_scanpy_data()

    print(f"Scanpy rank_genes_groups structure:")
    print(f"  Keys: {list(rank_genes_dict.keys())}")
    print(f"  Clusters: {cluster_info['clusters']}")

    # Test auto-detection and processing
    print("\nCalling get_top_markers()...")
    result = CASSIA.get_top_markers(rank_genes_dict, n_genes=15)

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

    return len(errors) == 0, errors, result


def test_get_top_markers_with_scanpy_df():
    """Test get_top_markers() with Scanpy flat DataFrame format."""
    print("\n--- Testing get_top_markers() with scanpy_df format ---")

    # Load scanpy_df
    scanpy_df = load_scanpy_df()

    print(f"scanpy_df structure:")
    print(f"  Shape: {scanpy_df.shape}")
    print(f"  Columns: {scanpy_df.columns.tolist()}")
    print(f"  Clusters: {scanpy_df['group'].unique().tolist()}")

    # Test auto-detection and processing
    print("\nCalling get_top_markers() with scanpy_df...")
    result = CASSIA.get_top_markers(scanpy_df, n_genes=15)

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

    return len(errors) == 0, errors, result


def test_runCASSIA_batch_with_scanpy(marker_df):
    """Test runCASSIA_batch with Scanpy-derived markers."""
    print("\n--- Testing runCASSIA_batch() with Scanpy-derived markers ---")

    # Load configuration
    config = load_config()
    llm_config = config['llm']
    data_config = config['data']

    # Setup API keys
    setup_api_keys()

    print(f"\nInput marker DataFrame:")
    print(f"  Shape: {marker_df.shape}")
    print(f"  Clusters: {marker_df['cluster'].tolist()}")

    # Create results directory
    results = create_results_dir("13_scanpy_format", get_test_mode())
    logging_ctx = setup_logging(results["logs"])
    output_name = str(Path(results["outputs"]) / "scanpy_batch_results")
    print(f"\nResults will be saved to: {results['base']}")

    errors = []
    start_time = time.time()

    try:
        print("\nRunning runCASSIA_batch...")
        CASSIA.runCASSIA_batch(
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
    print_test_header("13 - Scanpy Format Input (PIP INSTALL MODE)")

    # Verify pip installation
    pip_info = verify_cassia_pip_install()
    print(f"\nCASSIA Installation Info:")
    print(f"  Version: {pip_info['version']}")
    print(f"  Location: {pip_info['location']}")
    print(f"  Is pip install: {pip_info['is_pip_install']}")

    if not pip_info['is_pip_install']:
        print("\nWARNING: CASSIA may not be from pip install!")

    # Load configuration
    config = load_config()
    print_config_summary(config)

    # Setup API keys
    setup_api_keys()

    all_tests_passed = True
    all_errors = []
    start_time = time.time()
    marker_df = None

    # Test 1: get_top_markers with Scanpy structured array
    print("\n" + "=" * 60)
    print("Test 1: get_top_markers() with Scanpy Structured Array")
    print("=" * 60)
    try:
        passed, errors, marker_df = test_get_top_markers_with_scanpy()
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
    except FileNotFoundError as e:
        print(f"\nSKIPPED: {e}")
        all_errors.append(str(e))

    # Test 2: get_top_markers with scanpy_df format
    print("\n" + "=" * 60)
    print("Test 2: get_top_markers() with scanpy_df Format")
    print("=" * 60)
    try:
        passed, errors, marker_df_scanpy = test_get_top_markers_with_scanpy_df()
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
    except FileNotFoundError as e:
        print(f"\nSKIPPED: {e}")
        all_errors.append(str(e))
        marker_df_scanpy = None

    # Test 3: runCASSIA_batch with Scanpy-derived markers
    results = None
    if marker_df is not None and len(marker_df) > 0:
        print("\n" + "=" * 60)
        print("Test 3: runCASSIA_batch() with Scanpy Markers")
        print("=" * 60)
        passed, errors, results = test_runCASSIA_batch_with_scanpy(marker_df)
        all_tests_passed = all_tests_passed and passed
        all_errors.extend(errors)
        print(f"\nResult: {'PASSED' if passed else 'FAILED'}")
    else:
        print("\n" + "=" * 60)
        print("Test 3: SKIPPED (no markers from previous test)")
        print("=" * 60)

    duration = time.time() - start_time

    # Save metadata
    if results is None:
        results = create_results_dir("13_scanpy_format", get_test_mode())

    try:
        rank_genes_dict, cluster_info = load_scanpy_data()
        clusters_tested = cluster_info['clusters']
    except:
        clusters_tested = []

    metadata = create_test_metadata(
        test_name="scanpy_format_install",
        config=config,
        duration_seconds=duration,
        status="passed" if all_tests_passed else "failed",
        clusters_tested=clusters_tested,
        errors=all_errors
    )
    metadata['pip_install_info'] = pip_info
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
