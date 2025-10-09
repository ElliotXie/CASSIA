#!/usr/bin/env python3
"""
Test: runCASSIA_batch - Both Processed and Unprocessed Data
Description: Test both processed and unprocessed data formats
Expected Runtime: 6-10 minutes

This test validates CASSIA batch annotation with:
1. Processed data (2-column format: cluster, markers)
2. Unprocessed data (FindAllMarkers format with statistics)

Usage:
    python test_both.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv and unprocessed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd

# Add parent directories to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import CASSIA

shared_dir = Path(__file__).parent.parent / "shared"
sys.path.insert(0, str(shared_dir))

from test_config import load_config, validate_api_key
from test_utils import (
    setup_logging,
    log_test_start,
    log_test_end,
    save_results,
    Timer,
    get_timestamp
)
from sample_data import SampleDataLoader


def run_test_with_config(config_file, logger, test_dir, results_dir):
    """Run a single test with given config."""

    config = load_config(test_dir / config_file)
    logger.info("\n" + "=" * 80)
    logger.info(f"TESTING: {config['description']}")
    logger.info("=" * 80)
    logger.info(f"Config: {config_file}")
    logger.info(f"Data: {config['data_file']}")

    # Load data
    loader = SampleDataLoader()
    dataset_name = config['data_file'].replace('.csv', '')

    if dataset_name == 'processed':
        marker_data = loader.load_processed()
    elif dataset_name == 'unprocessed':
        marker_data = loader.load_unprocessed()
    else:
        raise ValueError(f"Unknown data file: {config['data_file']}")

    logger.info(f"✓ Loaded {len(marker_data)} rows from {config['data_file']}")
    logger.info(f"  Columns: {list(marker_data.columns)}")

    # Check for ranking columns if unprocessed
    if dataset_name == 'unprocessed':
        if 'avg_log2FC' in marker_data.columns:
            logger.info(f"  ✓ Found avg_log2FC column (unprocessed data)")
        else:
            logger.warning(f"  ⚠ avg_log2FC column not found!")

    # Get unique clusters
    if 'cluster' in marker_data.columns:
        clusters = marker_data['cluster'].unique()
        logger.info(f"  - Clusters: {len(clusters)}")
        logger.info(f"  - Names: {clusters[:5].tolist()}{' ...' if len(clusters) > 5 else ''}")

    # Run batch annotation
    start_time = log_test_start(logger, config.to_dict())

    try:
        with Timer(f"CASSIA batch ({dataset_name})", logger):
            result = CASSIA.runCASSIA_batch(
                marker=marker_data,
                output_name=config['method_params']['output_name'],
                model=config['model'],
                provider=config['provider'],
                temperature=config['temperature'],
                tissue=config['method_params']['tissue'],
                species=config['method_params']['species'],
                max_workers=config['method_params']['max_workers'],
                n_genes=config['method_params']['n_genes'],
                ranking_method=config['method_params'].get('ranking_method'),
                ascending=config['method_params'].get('ascending'),
                max_retries=config['method_params'].get('max_retries', 1)
            )

        # Read results
        result_csv = f"{config['method_params']['output_name']}_full.csv"
        if not os.path.exists(result_csv):
            raise FileNotFoundError(f"Results not found: {result_csv}")

        result_df = pd.read_csv(result_csv)
        logger.info(f"\n✓ Annotation completed")
        logger.info(f"  - Results shape: {result_df.shape}")
        logger.info(f"  - Output file: {result_csv}")

        # Validate
        logger.info("\nValidating results:")
        expected_cols = config['validation']['expected_columns']
        missing_cols = set(expected_cols) - set(result_df.columns)

        if missing_cols:
            raise ValueError(f"Missing columns: {missing_cols}")

        logger.info(f"  ✓ All expected columns present")
        logger.info(f"  ✓ Row count: {len(result_df)}/{config['validation']['min_rows']} (min)")

        # Save results
        timestamp = get_timestamp()
        archived_path = save_results(
            result_df,
            config['test_name'],
            results_dir,
            prefix=dataset_name,
            suffix="results"
        )
        logger.info(f"  ✓ Archived: {archived_path.name}")

        # Clean up
        try:
            os.remove(result_csv)
            summary_csv = f"{config['method_params']['output_name']}_summary.csv"
            if os.path.exists(summary_csv):
                os.remove(summary_csv)
            logger.info(f"  ✓ Cleaned up working files")
        except:
            pass

        log_test_end(logger, start_time, success=True)

        return {
            'test': config['test_name'],
            'data': dataset_name,
            'status': 'success',
            'rows': len(result_df),
            'file': archived_path.name
        }

    except Exception as e:
        logger.error(f"\n❌ Test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return {
            'test': config['test_name'],
            'data': dataset_name,
            'status': 'failed',
            'error': str(e)
        }


def main():
    """Main test function."""
    test_dir = Path(__file__).parent

    # Setup logging
    timestamp = get_timestamp()
    log_file = test_dir / "results" / f"{timestamp}_test_both_log.txt"
    logger = setup_logging("test_both_formats", log_file=log_file)

    logger.info("=" * 80)
    logger.info("CASSIA Test: runCASSIA_batch - BOTH DATA FORMATS")
    logger.info("=" * 80)
    logger.info("Testing both processed and unprocessed data formats")

    # Validate API key
    try:
        config = load_config(test_dir / "config.json")
        api_key = validate_api_key(config['provider'])
        CASSIA.set_api_key(api_key, provider=config['provider'])
        logger.info(f"✓ API key validated for provider: {config['provider']}")
    except EnvironmentError as e:
        logger.error(str(e))
        return 1

    # Create results directory
    results_dir = test_dir / "results"
    results_dir.mkdir(exist_ok=True)

    # Run both tests
    results = []

    # Test 1: Processed data
    logger.info("\n" + "=" * 80)
    logger.info("TEST 1/2: PROCESSED DATA (2-column format)")
    logger.info("=" * 80)
    result1 = run_test_with_config("config.json", logger, test_dir, results_dir)
    results.append(result1)

    # Test 2: Unprocessed data
    logger.info("\n" + "=" * 80)
    logger.info("TEST 2/2: UNPROCESSED DATA (FindAllMarkers format)")
    logger.info("=" * 80)
    result2 = run_test_with_config("config_unprocessed.json", logger, test_dir, results_dir)
    results.append(result2)

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("FINAL SUMMARY")
    logger.info("=" * 80)

    passed = sum(1 for r in results if r['status'] == 'success')
    failed = sum(1 for r in results if r['status'] == 'failed')

    logger.info(f"Total tests: {len(results)}")
    logger.info(f"  ✓ Passed: {passed}")
    logger.info(f"  ✗ Failed: {failed}")

    logger.info("\nDetailed results:")
    for r in results:
        status_icon = "✓" if r['status'] == 'success' else "✗"
        logger.info(f"  {status_icon} {r['data']:12} - {r['status']}")
        if r['status'] == 'success':
            logger.info(f"      Rows: {r['rows']}, File: {r['file']}")
        else:
            logger.info(f"      Error: {r.get('error', 'Unknown')}")

    # Save summary
    summary = {
        'timestamp': timestamp,
        'total_tests': len(results),
        'passed': passed,
        'failed': failed,
        'results': results
    }

    summary_path = save_results(
        summary,
        "test_both_formats",
        results_dir,
        prefix="both",
        suffix="summary"
    )
    logger.info(f"\n✓ Summary saved: {summary_path.name}")

    if failed == 0:
        logger.info("\n" + "=" * 80)
        logger.info("✅ ALL TESTS PASSED")
        logger.info(f"📁 Results: {results_dir}")
        logger.info(f"📄 Log: {log_file}")
        logger.info("=" * 80)
        return 0
    else:
        logger.info("\n" + "=" * 80)
        logger.info(f"❌ {failed} TEST(S) FAILED")
        logger.info(f"📁 Results: {results_dir}")
        logger.info(f"📄 Log: {log_file}")
        logger.info("=" * 80)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
