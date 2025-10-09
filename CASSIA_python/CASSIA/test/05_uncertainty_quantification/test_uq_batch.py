#!/usr/bin/env python3
"""
Test: Uncertainty Quantification
Description: Test annotation stability through multiple stochastic runs
Expected Runtime: 15-20 minutes

This test validates annotation stability and confidence by running batch
annotation multiple times with stochastic settings (temperature > 0) and
analyzing the consistency of results across runs.

Usage:
    python test_uq_batch.py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/processed.csv
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter
import json

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


def run_single_batch_iteration(marker_data, config, logger, run_number):
    """
    Run a single batch annotation iteration.

    Args:
        marker_data: Input marker DataFrame
        config: Test configuration
        logger: Logger instance
        run_number: Iteration number (1-indexed)

    Returns:
        pd.DataFrame: Batch annotation results
    """
    logger.info(f"\n  Running iteration {run_number}/{config['method_params']['n_iterations']}...")

    output_name = f"{config['method_params']['output_name']}_run_{run_number}"

    # Run batch annotation
    CASSIA.runCASSIA_batch(
        marker=marker_data,
        output_name=output_name,
        model=config['model'],
        provider=config['provider'],
        temperature=config['temperature'],
        tissue=config['method_params']['tissue'],
        species=config['method_params']['species'],
        max_workers=config['method_params']['max_workers'],
        n_genes=config['method_params']['n_genes']
    )

    # Read results
    result_file = f"{output_name}_full.csv"
    if not os.path.exists(result_file):
        raise FileNotFoundError(f"Run {run_number} results not found: {result_file}")

    result_df = pd.read_csv(result_file)
    logger.info(f"    ✓ Iteration {run_number} completed: {len(result_df)} clusters annotated")

    # Clean up intermediate files
    try:
        os.remove(result_file)
        summary_file = f"{output_name}_summary.csv"
        if os.path.exists(summary_file):
            os.remove(summary_file)
    except:
        pass

    return result_df


def calculate_pairwise_similarity(results_list, logger):
    """
    Calculate pairwise similarity between annotation runs.

    Args:
        results_list: List of DataFrames from each run
        logger: Logger instance

    Returns:
        np.ndarray: Similarity matrix
    """
    n_runs = len(results_list)
    similarity_matrix = np.zeros((n_runs, n_runs))

    logger.info("\nCalculating pairwise similarity scores:")

    for i in range(n_runs):
        for j in range(n_runs):
            if i == j:
                similarity_matrix[i, j] = 1.0
            elif i < j:  # Only calculate upper triangle
                # Calculate exact match similarity
                annotations_i = results_list[i]['Predicted Main Cell Type'].tolist()
                annotations_j = results_list[j]['Predicted Main Cell Type'].tolist()

                matches = sum(1 for a, b in zip(annotations_i, annotations_j) if a == b)
                similarity = matches / len(annotations_i)

                similarity_matrix[i, j] = similarity
                similarity_matrix[j, i] = similarity  # Mirror to lower triangle

                logger.info(f"  Run {i+1} vs Run {j+1}: {similarity:.3f} ({matches}/{len(annotations_i)} matches)")

    return similarity_matrix


def calculate_stability_metrics(results_list, similarity_matrix, logger):
    """
    Calculate stability and confidence metrics.

    Args:
        results_list: List of DataFrames from each run
        similarity_matrix: Pairwise similarity matrix
        logger: Logger instance

    Returns:
        dict: Stability metrics
    """
    logger.info("\nCalculating stability metrics:")

    # Overall average similarity (excluding diagonal)
    n_runs = len(results_list)
    upper_triangle_indices = np.triu_indices(n_runs, k=1)
    avg_similarity = similarity_matrix[upper_triangle_indices].mean()
    logger.info(f"  - Average pairwise similarity: {avg_similarity:.3f}")

    # Per-cluster consistency
    cluster_consistency = {}
    clusters = results_list[0]['True Cell Type'].tolist()

    for cluster_idx, cluster_name in enumerate(clusters):
        # Get all annotations for this cluster across runs
        annotations = [df.iloc[cluster_idx]['Predicted Main Cell Type'] for df in results_list]

        # Count most common annotation
        annotation_counts = Counter(annotations)
        most_common, count = annotation_counts.most_common(1)[0]

        consistency = count / n_runs
        cluster_consistency[cluster_name] = {
            'consistency': consistency,
            'most_common': most_common,
            'count': count,
            'total_runs': n_runs,
            'unique_annotations': len(annotation_counts)
        }

        logger.info(f"  - {cluster_name}: {consistency:.3f} consistency ({count}/{n_runs} agreement on '{most_common}')")

    # Overall stability score
    stability_score = np.mean([v['consistency'] for v in cluster_consistency.values()])
    logger.info(f"\n  Overall stability score: {stability_score:.3f}")

    metrics = {
        'n_runs': n_runs,
        'avg_similarity': float(avg_similarity),
        'stability_score': float(stability_score),
        'cluster_consistency': cluster_consistency,
        'similarity_range': {
            'min': float(similarity_matrix[upper_triangle_indices].min()),
            'max': float(similarity_matrix[upper_triangle_indices].max()),
            'std': float(similarity_matrix[upper_triangle_indices].std())
        }
    }

    return metrics


def create_consensus_annotations(results_list, logger):
    """
    Create consensus annotations from multiple runs.

    Args:
        results_list: List of DataFrames from each run
        logger: Logger instance

    Returns:
        pd.DataFrame: Consensus annotations with confidence scores
    """
    logger.info("\nCreating consensus annotations:")

    # Start with the first run's structure
    consensus_df = results_list[0][['True Cell Type']].copy()

    clusters = consensus_df['True Cell Type'].tolist()
    consensus_annotations = []
    confidence_scores = []

    for cluster_idx, cluster_name in enumerate(clusters):
        # Get all annotations for this cluster
        annotations = [df.iloc[cluster_idx]['Predicted Main Cell Type'] for df in results_list]

        # Find most common annotation
        annotation_counts = Counter(annotations)
        most_common, count = annotation_counts.most_common(1)[0]

        # Calculate confidence as proportion of runs agreeing
        confidence = count / len(results_list)

        consensus_annotations.append(most_common)
        confidence_scores.append(confidence)

        logger.info(f"  - {cluster_name}: '{most_common}' (confidence: {confidence:.3f})")

    consensus_df['Consensus_Annotation'] = consensus_annotations
    consensus_df['Confidence_Score'] = confidence_scores

    return consensus_df


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
    logger.info("CASSIA Test: Uncertainty Quantification (Multiple Stochastic Runs)")
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

        # Load full dataset first
        full_data = loader.load_processed() if dataset_name == 'processed' else loader.load_unprocessed()

        # Use subset for UQ testing (faster)
        n_subset = config['method_params']['n_clusters_subset']
        marker_data = loader.load_subset(
            dataset=dataset_name,
            n_clusters=n_subset,
            random_state=42
        )

        logger.info(f"✓ Loaded subset: {len(marker_data)} rows ({n_subset} clusters)")

        cluster_col = marker_data.columns[1]
        clusters = marker_data[cluster_col].unique()
        logger.info(f"  - Clusters: {clusters.tolist()}")
    except Exception as e:
        logger.error(f"Failed to load data: {str(e)}", exc_info=True)
        return 1

    # Start test
    start_time = log_test_start(logger, config.to_dict())

    try:
        # Run multiple iterations
        logger.info("\n" + "=" * 80)
        logger.info("RUNNING MULTIPLE STOCHASTIC ITERATIONS")
        logger.info("=" * 80)
        logger.info(f"\nUQ configuration:")
        logger.info(f"  - Number of runs: {config['method_params']['n_iterations']}")
        logger.info(f"  - Temperature: {config['temperature']} (stochastic)")
        logger.info(f"  - Clusters: {n_subset}")
        logger.info(f"  - Model: {config['model']}")

        results_list = []

        with Timer(f"All {config['method_params']['n_iterations']} UQ iterations", logger):
            for run_num in range(1, config['method_params']['n_iterations'] + 1):
                result_df = run_single_batch_iteration(marker_data, config, logger, run_num)
                results_list.append(result_df)

        logger.info(f"\n✓ All {len(results_list)} iterations completed")

        # Calculate similarity matrix
        logger.info("\n" + "=" * 80)
        logger.info("ANALYZING ANNOTATION STABILITY")
        logger.info("=" * 80)

        similarity_matrix = calculate_pairwise_similarity(results_list, logger)

        # Calculate stability metrics
        stability_metrics = calculate_stability_metrics(results_list, similarity_matrix, logger)

        # Create consensus annotations
        consensus_df = create_consensus_annotations(results_list, logger)

        # Validate results
        logger.info("\n" + "=" * 80)
        logger.info("VALIDATING UQ RESULTS")
        logger.info("=" * 80)

        if config['validation']['check_all_runs_complete']:
            if len(results_list) == config['method_params']['n_iterations']:
                logger.info(f"  ✓ All {len(results_list)} runs completed")
            else:
                raise ValueError(f"Expected {config['method_params']['n_iterations']} runs, got {len(results_list)}")

        if config['validation']['check_similarity_scores']:
            avg_sim = stability_metrics['avg_similarity']
            min_sim = config['validation']['min_similarity']
            max_sim = config['validation']['max_similarity']

            if min_sim <= avg_sim <= max_sim:
                logger.info(f"  ✓ Average similarity within expected range: {avg_sim:.3f} ∈ [{min_sim}, {max_sim}]")
            else:
                logger.warning(f"  ⚠ Average similarity outside expected range: {avg_sim:.3f} ∉ [{min_sim}, {max_sim}]")

        if config['validation']['check_stability_metrics']:
            stability_score = stability_metrics['stability_score']
            logger.info(f"  ✓ Stability score calculated: {stability_score:.3f}")

        # Archive results
        logger.info("\n" + "=" * 80)
        logger.info("ARCHIVING RESULTS")
        logger.info("=" * 80)

        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        # Save individual run results
        if config['reporting']['save_individual_runs']:
            for i, result_df in enumerate(results_list, 1):
                run_path = save_results(
                    result_df,
                    config['test_name'],
                    results_dir,
                    prefix=f"run_{i}",
                    suffix="annotations"
                )
                logger.info(f"  ✓ Archived run {i}: {run_path.name}")

        # Save similarity matrix
        if config['reporting']['save_similarity_matrix']:
            similarity_df = pd.DataFrame(
                similarity_matrix,
                columns=[f"Run_{i+1}" for i in range(len(results_list))],
                index=[f"Run_{i+1}" for i in range(len(results_list))]
            )
            sim_path = save_results(
                similarity_df,
                config['test_name'],
                results_dir,
                prefix="similarity",
                suffix="matrix"
            )
            logger.info(f"  ✓ Archived similarity matrix: {sim_path.name}")

        # Save consensus annotations
        consensus_path = save_results(
            consensus_df,
            config['test_name'],
            results_dir,
            prefix="consensus",
            suffix="annotations"
        )
        logger.info(f"  ✓ Archived consensus annotations: {consensus_path.name}")

        # Save stability metrics
        metrics_path = save_results(
            stability_metrics,
            config['test_name'],
            results_dir,
            prefix="stability",
            suffix="metrics"
        )
        logger.info(f"  ✓ Archived stability metrics: {metrics_path.name}")

        # Log final summary
        logger.info("\n" + "=" * 80)
        logger.info("UNCERTAINTY QUANTIFICATION SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Clusters analyzed: {n_subset}")
        logger.info(f"Stochastic runs: {len(results_list)}")
        logger.info(f"Average similarity: {stability_metrics['avg_similarity']:.3f}")
        logger.info(f"Overall stability: {stability_metrics['stability_score']:.3f}")
        logger.info(f"Similarity range: {stability_metrics['similarity_range']['min']:.3f} - {stability_metrics['similarity_range']['max']:.3f}")

        logger.info("\nPer-cluster confidence:")
        for cluster, metrics in stability_metrics['cluster_consistency'].items():
            logger.info(f"  - {cluster}: {metrics['consistency']:.3f} ({metrics['count']}/{metrics['total_runs']} consensus)")

        log_test_end(logger, start_time, success=True)

        logger.info("\n" + "=" * 80)
        logger.info("✅ UNCERTAINTY QUANTIFICATION TEST COMPLETED SUCCESSFULLY")
        logger.info(f"📁 Results directory: {results_dir}")
        logger.info(f"📄 Log file: {log_file}")
        logger.info("=" * 80)

        return 0

    except Exception as e:
        logger.error(f"\n❌ UQ test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
