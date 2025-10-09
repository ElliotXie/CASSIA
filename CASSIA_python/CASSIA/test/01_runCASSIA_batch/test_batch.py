#!/usr/bin/env python3
"""
Test: runCASSIA_batch
Description: Test basic batch annotation functionality
Expected Runtime: 3-5 minutes

This test validates the core batch annotation functionality of CASSIA
using clean sample data from the processed.csv dataset.

Usage:
    python test_batch.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd

# Add parent directories to path for imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Import CASSIA
import CASSIA

# Import shared utilities
shared_dir = Path(__file__).parent.parent / "shared"
sys.path.insert(0, str(shared_dir))

from test_config import load_config, validate_api_key
from test_utils import (
    setup_logging,
    log_test_start,
    log_test_end,
    save_results,
    validate_output,
    Timer,
    get_timestamp
)
from sample_data import load_sample_data


def main():
    """Main test function."""
    # Setup
    test_dir = Path(__file__).parent
    config = load_config(test_dir / "config.json")

    # Setup logging with file output
    timestamp = get_timestamp()
    log_file = test_dir / "results" / f"{timestamp}_test_log.txt"
    logger = setup_logging(config['test_name'], log_file=log_file)

    logger.info("=" * 80)
    logger.info("CASSIA Test: runCASSIA_batch")
    logger.info("=" * 80)

    # Validate API key
    try:
        api_key = validate_api_key(config['provider'])
        CASSIA.set_api_key(api_key, provider=config['provider'])
        logger.info(f"✓ API key validated for provider: {config['provider']}")
    except EnvironmentError as e:
        logger.error(str(e))
        return 1

    # Load data
    logger.info(f"Loading data: {config['data_file']}")
    try:
        # Extract base name without extension
        dataset_name = config['data_file'].replace('.csv', '')
        marker_data = load_sample_data(dataset_name)
        logger.info(f"✓ Loaded {len(marker_data)} rows from {config['data_file']}")

        # Log dataset info
        if len(marker_data.columns) >= 2:
            cluster_col = marker_data.columns[1]
            n_clusters = marker_data[cluster_col].nunique()
            logger.info(f"  - Number of clusters: {n_clusters}")
            logger.info(f"  - Cluster names: {marker_data[cluster_col].unique().tolist()}")
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}", exc_info=True)
        return 1

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    try:
        # Run CASSIA batch annotation
        logger.info("\nStarting CASSIA batch annotation...")

        with Timer("CASSIA batch annotation", logger):
            result = CASSIA.runCASSIA_batch(
                marker=marker_data,
                model=config['model'],
                provider=config['provider'],
                temperature=config['temperature'],
                **config['method_params']
            )

        # The runCASSIA_batch function saves files automatically
        # We need to read them back to validate
        logger.info("\nReading generated results...")

        # Construct expected output filenames
        output_name = config['method_params']['output_name']
        full_csv_path = f"{output_name}_full.csv"
        summary_csv_path = f"{output_name}_summary.csv"

        # Check if files were created
        if not os.path.exists(full_csv_path):
            raise FileNotFoundError(f"Expected output file not found: {full_csv_path}")

        # Read the results
        result_df = pd.read_csv(full_csv_path)
        summary_df = pd.read_csv(summary_csv_path)

        logger.info(f"✓ Read full results: {result_df.shape}")
        logger.info(f"✓ Read summary results: {summary_df.shape}")

        # Validate results
        if config['validation']['check_output_format']:
            logger.info("\nValidating output format...")

            validate_output(
                result_df,
                expected_columns=config['validation']['expected_columns'],
                min_rows=config['validation'].get('min_rows', 1),
                max_nulls=config['validation'].get('max_nulls')
            )
            logger.info("✓ Validation passed")

        # Save timestamped copies to results directory
        logger.info("\nSaving timestamped results...")
        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Save both full and summary results
        full_path = save_results(
            result_df,
            config['test_name'],
            results_dir,
            prefix="batch",
            suffix="full"
        )
        logger.info(f"✓ Saved full results: {full_path}")

        summary_path = save_results(
            summary_df,
            config['test_name'],
            results_dir,
            prefix="batch",
            suffix="summary"
        )
        logger.info(f"✓ Saved summary results: {summary_path}")

        # Clean up temporary files
        try:
            os.remove(full_csv_path)
            os.remove(summary_csv_path)
            logger.info("✓ Cleaned up temporary files")
        except:
            pass

        # Log summary statistics
        logger.info("\n" + "=" * 80)
        logger.info("TEST RESULTS SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Total clusters processed: {len(result_df)}")
        logger.info(f"Successful annotations: {result_df['Predicted Main Cell Type'].notna().sum()}")

        if 'Iterations' in result_df.columns:
            avg_iterations = result_df['Iterations'].mean()
            logger.info(f"Average iterations: {avg_iterations:.2f}")

        # Show sample results
        logger.info("\nSample annotations:")
        for idx, row in result_df.head(3).iterrows():
            cluster = row['True Cell Type']
            prediction = row['Predicted Main Cell Type']
            logger.info(f"  {cluster} → {prediction}")

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Results saved to: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ Test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
