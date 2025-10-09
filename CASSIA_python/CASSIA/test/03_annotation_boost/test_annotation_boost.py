#!/usr/bin/env python3
"""
Test: runCASSIA_annotationboost
Description: Test iterative marker analysis for deep-dive annotation
Expected Runtime: 5-10 minutes

This test validates the annotation boost functionality which performs
iterative marker analysis for ambiguous or uncertain clusters. The LLM
can request additional markers iteratively to refine its annotation.

Usage:
    python test_annotation_boost.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import json
import glob

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
    Timer,
    get_timestamp
)
from sample_data import SampleDataLoader


def create_prerequisite_batch_results(logger, config, marker_data):
    """
    Create the prerequisite batch annotation results needed for boost.

    Returns:
        str: Path to the batch results CSV
    """
    logger.info("\n" + "=" * 80)
    logger.info("CREATING PREREQUISITE BATCH RESULTS")
    logger.info("=" * 80)

    batch_output = config['prerequisites']['batch_results_file'].replace('_full.csv', '')

    logger.info(f"Running batch annotation to create prerequisite data...")
    logger.info(f"Output name: {batch_output}")

    with Timer("Prerequisite batch annotation", logger):
        CASSIA.runCASSIA_batch(
            marker=marker_data,
            output_name=batch_output,
            model=config['model'],
            provider=config['provider'],
            tissue=config['method_params']['major_cluster_info'].split()[-1],  # Extract tissue
            species=config['method_params']['major_cluster_info'].split()[0],  # Extract species
            max_workers=4,
            n_genes=50
        )

    batch_csv = f"{batch_output}_full.csv"

    if not os.path.exists(batch_csv):
        raise FileNotFoundError(f"Batch results not created: {batch_csv}")

    logger.info(f"✓ Created prerequisite batch results: {batch_csv}")
    return batch_csv


def validate_boost_outputs(output_name, logger, config):
    """Validate that all expected boost outputs were created."""

    # Expected output files
    expected_files = {
        'conversation_json': f"{output_name}_conversation.json",
        'html_summary': f"{output_name}_summary.html",
        'raw_text': f"{output_name}_raw_conversation.txt"
    }

    found_files = {}
    missing_files = []

    for file_type, filename in expected_files.items():
        if os.path.exists(filename):
            found_files[file_type] = filename
            logger.info(f"  ✓ Found {file_type}: {filename}")
        else:
            missing_files.append(filename)
            logger.warning(f"  ✗ Missing {file_type}: {filename}")

    if missing_files and config['validation']['check_html_report']:
        raise FileNotFoundError(f"Missing output files: {missing_files}")

    return found_files


def validate_conversation_history(conv_json_path, logger, config):
    """Validate the conversation history JSON structure."""

    with open(conv_json_path, 'r') as f:
        conversation = json.load(f)

    logger.info(f"\nConversation history structure:")
    logger.info(f"  - Type: {type(conversation)}")

    if isinstance(conversation, list):
        logger.info(f"  - Messages: {len(conversation)}")

        # Validate iterations (look for check_genes requests)
        iterations = 0
        for msg in conversation:
            if isinstance(msg, dict) and 'content' in msg:
                if '<check_genes>' in str(msg['content']):
                    iterations += 1

        logger.info(f"  - Iterations with gene checks: {iterations}")

        if config['validation']['min_iterations']:
            if iterations < config['validation']['min_iterations']:
                logger.warning(f"  ⚠ Fewer iterations than minimum ({iterations} < {config['validation']['min_iterations']})")

        if config['validation']['max_iterations']:
            if iterations > config['validation']['max_iterations']:
                logger.warning(f"  ⚠ More iterations than maximum ({iterations} > {config['validation']['max_iterations']})")

        # Show sample messages
        logger.info(f"\nSample conversation messages:")
        for i, msg in enumerate(conversation[:3]):
            if isinstance(msg, dict):
                role = msg.get('role', 'unknown')
                content_preview = str(msg.get('content', ''))[:100] + "..."
                logger.info(f"  [{i}] {role}: {content_preview}")

        return True

    return False


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
    logger.info("CASSIA Test: runCASSIA_annotationboost (Iterative Deep-Dive)")
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
    logger.info(f"\nLoading data: {config['data_file']}")
    try:
        loader = SampleDataLoader()
        dataset_name = config['data_file'].replace('.csv', '')
        marker_data = loader.load_processed() if dataset_name == 'processed' else loader.load_unprocessed()

        logger.info(f"✓ Loaded {len(marker_data)} rows from {config['data_file']}")

        # Get cluster info
        cluster_col = marker_data.columns[1]
        clusters = marker_data[cluster_col].unique()
        logger.info(f"  - Available clusters: {clusters.tolist()}")

        # Validate test cluster exists
        test_cluster = config['test_cluster']
        if test_cluster not in clusters.tolist():
            raise ValueError(f"Test cluster '{test_cluster}' not found in data")

        logger.info(f"  - Test cluster: {test_cluster}")
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}", exc_info=True)
        return 1

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    try:
        # Step 1: Create prerequisite batch results
        if config['prerequisites']['requires_batch_results']:
            batch_csv = create_prerequisite_batch_results(logger, config, marker_data)
        else:
            batch_csv = config['prerequisites']['batch_results_file']
            if not os.path.exists(batch_csv):
                raise FileNotFoundError(f"Batch results file not found: {batch_csv}")

        # Step 2: Run annotation boost
        logger.info("\n" + "=" * 80)
        logger.info("RUNNING ANNOTATION BOOST")
        logger.info("=" * 80)
        logger.info(f"\nBoost configuration:")
        logger.info(f"  - Cluster: {config['method_params']['cluster_name']}")
        logger.info(f"  - Context: {config['method_params']['major_cluster_info']}")
        logger.info(f"  - Max iterations: {config['method_params']['num_iterations']}")
        logger.info(f"  - Conversation mode: {config['method_params']['conversation_history_mode']}")
        logger.info(f"  - Report style: {config['method_params']['report_style']}")

        with Timer("Annotation boost analysis", logger):
            result, messages = CASSIA.runCASSIA_annotationboost(
                full_result_path=batch_csv,
                marker=marker_data,
                cluster_name=config['method_params']['cluster_name'],
                major_cluster_info=config['method_params']['major_cluster_info'],
                output_name=config['method_params']['output_name'],
                num_iterations=config['method_params']['num_iterations'],
                model=config['model'],
                provider=config['provider'],
                temperature=config['temperature'],
                conversation_history_mode=config['method_params']['conversation_history_mode'],
                report_style=config['method_params']['report_style']
            )

        logger.info(f"\n✓ Annotation boost completed")
        logger.info(f"  - Result length: {len(str(result))} characters")
        logger.info(f"  - Messages count: {len(messages) if isinstance(messages, list) else 'N/A'}")

        # Step 3: Validate outputs
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATING BOOST OUTPUTS")
        logger.info("=" * 80)

        found_files = validate_boost_outputs(
            config['method_params']['output_name'],
            logger,
            config
        )

        # Step 4: Validate conversation history
        if config['validation']['check_conversation_history'] and 'conversation_json' in found_files:
            logger.info("\nValidating conversation history:")
            validate_conversation_history(
                found_files['conversation_json'],
                logger,
                config
            )

        # Step 5: Archive results
        logger.info("\n" + "=" * 80)
        logger.info("ARCHIVING RESULTS")
        logger.info("=" * 80)

        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Copy all boost outputs to results directory
        for file_type, filepath in found_files.items():
            if os.path.exists(filepath):
                dest_path = results_dir / f"{timestamp}_{file_type}{Path(filepath).suffix}"
                import shutil
                shutil.copy2(filepath, dest_path)
                logger.info(f"  ✓ Archived {file_type}: {dest_path.name}")

                # Clean up original
                try:
                    os.remove(filepath)
                except:
                    pass

        # Clean up batch results
        try:
            if os.path.exists(batch_csv):
                os.remove(batch_csv)
            batch_summary = batch_csv.replace('_full.csv', '_summary.csv')
            if os.path.exists(batch_summary):
                os.remove(batch_summary)
            logger.info(f"  ✓ Cleaned up prerequisite files")
        except:
            pass

        # Create test summary
        summary = {
            'test_cluster': config['test_cluster'],
            'boost_outputs': list(found_files.keys()),
            'conversation_messages': len(messages) if isinstance(messages, list) else 0,
            'result_length': len(str(result))
        }

        summary_path = save_results(
            summary,
            config['test_name'],
            results_dir,
            prefix="boost",
            suffix="summary"
        )
        logger.info(f"  ✓ Saved test summary: {summary_path.name}")

        # Log final summary
        logger.info("\n" + "=" * 80)
        logger.info("BOOST ANALYSIS SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Cluster analyzed: {config['test_cluster']}")
        logger.info(f"Boost outputs created: {len(found_files)}")
        logger.info(f"Files archived: {len(found_files) + 1}")  # +1 for summary

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ ANNOTATION BOOST TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Results directory: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ Boost test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
