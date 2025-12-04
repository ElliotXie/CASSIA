"""
Scoring functions for CASSIA cell type annotations.

This module provides functions for evaluating and scoring cell type
annotation results using LLM-based quality assessment.
"""

import re
import threading
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import dependencies
try:
    from CASSIA.core.llm_utils import call_llm
except ImportError:
    try:
        from ..core.llm_utils import call_llm
    except ImportError:
        from llm_utils import call_llm

try:
    from CASSIA.core.progress_tracker import BatchProgressTracker
except ImportError:
    try:
        from ..core.progress_tracker import BatchProgressTracker
    except ImportError:
        from progress_tracker import BatchProgressTracker

try:
    from CASSIA.core.utils import get_column_value, MARKER_COLUMN_OPTIONS, HISTORY_COLUMN_OPTIONS, CLUSTER_ID_COLUMN_OPTIONS
except ImportError:
    try:
        from ..core.utils import get_column_value, MARKER_COLUMN_OPTIONS, HISTORY_COLUMN_OPTIONS, CLUSTER_ID_COLUMN_OPTIONS
    except ImportError:
        from utils import get_column_value, MARKER_COLUMN_OPTIONS, HISTORY_COLUMN_OPTIONS, CLUSTER_ID_COLUMN_OPTIONS


def prompt_creator_score(major_cluster_info, marker, annotation_history):
    """
    Create a prompt for scoring cell type annotations.

    Args:
        major_cluster_info (str): Information about species and tissue
        marker (str): Comma-separated list of marker genes
        annotation_history (str): History of annotation conversation

    Returns:
        str: Formatted prompt for the scoring LLM
    """
    prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and rate single-cell annotation results, focusing on their correctness and ability to capture the overall picture of the data. You will provide a score from 0 to 100 and justify your rating.

Here are the single-cell annotation results to evaluate:



<marker>
{marker}
</marker>

<Cluster Origin>
{major_cluster_info}
</Cluster Origin>

<annotation_history>
{annotation_history}
</annotation_history>

Carefully analyze these results, paying particular attention to the following aspects:
1. Correctness of the annotations
2. Balanced consideration of multiple markers rather than over-focusing on a specific one
3. Ability to capture the general picture of the cell populations

When evaluating, consider:
- Are the annotations scientifically accurate?
- Is there a good balance in the use of different markers?
- Does the annotation provide a comprehensive view of the cell types present?
- Are there any obvious misclassifications or oversights?
- Did it consider the rank of the marker? marker appear first is more important.

Provide your analysis in the following format:
1. Start with a <reasoning> tag, where you explain your evaluation of the annotation results. Discuss the strengths and weaknesses you've identified, referring to specific examples from the results where possible.
2. After your reasoning, use a <score> tag to provide a numerical score from 0 to 100, where 0 represents completely incorrect or unusable results, and 100 represents perfect annotation that captures all aspects of the data correctly.

Your response should look like this:

<reasoning>
[Your detailed analysis and justification here]
</reasoning>

<score>[Your numerical score between 0 and 100]</score>

Remember, the focus is on correctness and the ability to see the general picture, rather than the structure of the results. Be critical but fair in your assessment.
    """
    return prompt


def extract_score_and_reasoning(text):
    """
    Extract both score and reasoning from annotation text.

    Args:
        text (str): Text containing score and reasoning between XML-like tags

    Returns:
        tuple: (score, reasoning_text) where score is int or None and reasoning_text is str or None

    Example:
        >>> score, reasoning = extract_score_and_reasoning("<reasoning>Good analysis</reasoning><score>85</score>")
        >>> print(f"Score: {score}, Reasoning: {reasoning[:20]}...")
        Score: 85, Reasoning: Good analysis...
    """
    try:
        # Initialize results
        score = None
        reasoning = None

        # Extract score - try multiple patterns
        score_patterns = [
            r'<score>(\d+)</score>',  # Original format
            r'Score:\s*(\d+)',        # "Score: 85"
            r'score:\s*(\d+)',        # "score: 85"
            r'(\d+)/100',             # "85/100"
            r'(\d+)\s*out\s*of\s*100', # "85 out of 100"
            r'rating.*?(\d+)',        # "rating of 85"
            r'(\d+)%'                 # "85%"
        ]

        for pattern in score_patterns:
            score_match = re.search(pattern, text, re.IGNORECASE)
            if score_match:
                score = int(score_match.group(1))
                break

        # Extract reasoning - try multiple patterns
        reasoning_patterns = [
            r'<reasoning>(.*?)</reasoning>',  # Original format
            r'Reasoning:\s*(.*?)(?=Score:|$)',  # "Reasoning: ..." until "Score:" or end
            r'reasoning:\s*(.*?)(?=score:|$)',  # lowercase version
            r'Analysis:\s*(.*?)(?=Score:|$)',   # "Analysis: ..."
            r'Evaluation:\s*(.*?)(?=Score:|$)' # "Evaluation: ..."
        ]

        for pattern in reasoning_patterns:
            reasoning_match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
            if reasoning_match:
                reasoning = reasoning_match.group(1).strip()
                break

        # If no specific reasoning found, use the entire text as reasoning
        if reasoning is None and text.strip():
            reasoning = text.strip()

        return score, reasoning

    except Exception as e:
        print(f"Error extracting data: {str(e)}")
        return None, None


def score_single_analysis(major_cluster_info, marker, annotation_history, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Score a single cell type annotation analysis.

    Args:
        major_cluster_info (str): Information about species and tissue
        marker (str): Comma-separated list of marker genes
        annotation_history (str): History of annotation conversation
        model (str): Model to use (e.g., "gpt-4" for OpenAI or "claude-3-5-sonnet-20241022" for Anthropic)
        provider (str): AI provider to use ('openai', 'anthropic', or 'openrouter')

    Returns:
        tuple: (score, reasoning) where score is int and reasoning is str
    """
    prompt = prompt_creator_score(major_cluster_info, marker, annotation_history)

    # Add explicit max_tokens to ensure responses aren't truncated
    response = call_llm(
        prompt=prompt,
        provider=provider,
        model=model,
        max_tokens=4096  # Maximum tokens allowed for most models
    )

    score, reasoning = extract_score_and_reasoning(response)
    return score, reasoning


def _get_task_name(row):
    """
    Get a task name from a row for progress tracking.

    Args:
        row: pandas Series row

    Returns:
        str: Task name (cluster ID or row number)
    """
    for col in CLUSTER_ID_COLUMN_OPTIONS:
        if col in row.index:
            return str(row[col])
    return f"Row {row.name + 1}"


def process_single_row(row_data, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Process a single row of data for scoring.

    Args:
        row_data (tuple): (idx, row) containing index and row data
        model (str): Model to use
        provider (str): AI provider to use ('openai', 'anthropic', or 'openrouter')

    Returns:
        tuple: (idx, score, reasoning)
    """
    idx, row = row_data

    try:
        major_cluster_info = f"{row['Species']} {row['Tissue']}"

        # Get marker using helper function
        marker = get_column_value(row, MARKER_COLUMN_OPTIONS)

        # Get annotation history using helper function
        annotation_history = get_column_value(row, HISTORY_COLUMN_OPTIONS)

        # Try up to 3 times for a valid score if we get None
        score, reasoning = None, None
        max_retries_for_none = 3
        retry_count = 0

        while score is None and retry_count < max_retries_for_none:
            score, reasoning = score_single_analysis(
                major_cluster_info,
                marker,
                annotation_history,
                model=model,
                provider=provider
            )

            if score is not None:
                break

            retry_count += 1

        return (idx, score, reasoning)

    except Exception as e:
        return (idx, None, f"Error: {str(e)}")


def runCASSIA_score_batch(input_file, output_file=None, max_workers=4, model="deepseek/deepseek-chat-v3-0324", provider="openrouter", max_retries=1, generate_report=True):
    """
    Run scoring on a batch of cell type annotations with progress tracking.

    Args:
        input_file (str): Path to input CSV file (with or without .csv extension)
        output_file (str, optional): Path to output CSV file (with or without .csv extension)
        max_workers (int): Maximum number of parallel workers
        model (str): Model to use
        provider (str): AI provider to use ('openai', 'anthropic', or 'openrouter')
        max_retries (int): Maximum number of retries for failed analyses
        generate_report (bool): Whether to generate HTML report (default True, set False when called from pipeline)

    Returns:
        pd.DataFrame: Results DataFrame with scores
    """
    # Add .csv extension if not present
    if not input_file.lower().endswith('.csv'):
        input_file = input_file + '.csv'

    if output_file and not output_file.lower().endswith('.csv'):
        output_file = output_file + '.csv'

    print(f"Starting scoring process with {max_workers} workers using {provider} ({model})...")

    try:
        # Read the input file
        results = pd.read_csv(input_file)

        # Initialize new columns if they don't exist
        if 'Score' not in results.columns:
            results['Score'] = None
        if 'Scoring_Reasoning' not in results.columns:
            results['Scoring_Reasoning'] = None

        # Determine output file path (before early return check)
        actual_output_file = output_file if output_file else input_file.replace('.csv', '_scored.csv')

        # Create a list of unscored rows to process
        rows_to_process = [
            (idx, row) for idx, row in results.iterrows()
            if pd.isna(row['Score'])
        ]

        if not rows_to_process:
            print("All rows already scored!")
            # Save the results to output file before returning
            results.to_csv(actual_output_file, index=False)
            return results

        # Initialize progress tracker
        total_to_score = len(rows_to_process)
        tracker = BatchProgressTracker(total_to_score)
        print(f"\nStarting scoring of {total_to_score} rows with {max_workers} parallel workers...\n")

        # Set up a lock for DataFrame updates
        df_lock = threading.Lock()
        failed_analyses = []  # Track failed rows for reporting

        # Define a function that includes retry logic and progress tracking
        def process_with_retry(row_data):
            idx, row = row_data
            task_name = _get_task_name(row)
            tracker.start_task(task_name)

            for attempt in range(max_retries + 1):
                try:
                    result = process_single_row(row_data, model=model, provider=provider)
                    tracker.complete_task(task_name)
                    return result
                except Exception as exc:
                    # Don't retry authentication errors
                    if "401" in str(exc) or "API key" in str(exc) or "authentication" in str(exc).lower():
                        tracker.complete_task(task_name)
                        return idx, None, f"API error: {str(exc)}"

                    # For other errors, retry if attempts remain
                    if attempt >= max_retries:
                        tracker.complete_task(task_name)
                        return idx, None, f"Failed after {max_retries + 1} attempts: {str(exc)}"

        # Process rows in parallel
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_row = {
                executor.submit(process_with_retry, row_data): row_data
                for row_data in rows_to_process
            }

            # Process completed jobs
            for future in as_completed(future_to_row):
                row_data = future_to_row[future]
                idx, row = row_data
                task_name = _get_task_name(row)
                try:
                    idx, score, reasoning = future.result()

                    # Safely update DataFrame
                    with df_lock:
                        results.loc[idx, 'Score'] = score
                        results.loc[idx, 'Scoring_Reasoning'] = reasoning

                        # Track failures
                        if score is None:
                            failed_analyses.append((task_name, reasoning))

                        # Save intermediate results
                        results.to_csv(actual_output_file, index=False)
                except Exception as exc:
                    failed_analyses.append((task_name, str(exc)))

        # Finalize progress display
        tracker.finish()

        # Report any failures
        if failed_analyses:
            print(f"\nWarning: {len(failed_analyses)} scoring(s) failed:")
            for task_name, error in failed_analyses:
                print(f"  - {task_name}: {error[:100]}{'...' if len(error) > 100 else ''}")
            print()

        # Print summary statistics
        total_rows = len(results)
        scored_rows = results['Score'].notna().sum()
        print(f"Scoring completed. Results saved.")
        print(f"\nSummary:")
        print(f"Total rows: {total_rows}")
        print(f"Successfully scored: {scored_rows}")
        print(f"Failed/Skipped: {total_rows - scored_rows}")

        # Generate HTML report with scores (unless called from pipeline which handles its own)
        if generate_report:
            try:
                try:
                    from CASSIA.reports.generate_batch_report import generate_batch_html_report_from_data
                except ImportError:
                    try:
                        from ..reports.generate_batch_report import generate_batch_html_report_from_data
                    except ImportError:
                        from generate_batch_report import generate_batch_html_report_from_data

                # Convert DataFrame to list of dicts for report generation
                rows_data = results.to_dict('records')

                # Generate report path based on output file
                report_output_path = actual_output_file.replace('.csv', '_report.html')
                generate_batch_html_report_from_data(
                    rows=rows_data,
                    output_path=report_output_path,
                    report_title="CASSIA Scored Analysis Report"
                )
                print(f"HTML report generated: {report_output_path}")
            except Exception as e:
                print(f"Warning: Could not generate HTML report: {e}")

        return results

    except Exception as e:
        print(f"Error in runCASSIA_score_batch: {str(e)}")
        raise


def score_annotation_batch(results_file_path, output_file_path=None, max_workers=4, model="deepseek/deepseek-chat-v3-0324", provider="openrouter"):
    """
    Process and score all rows in a results CSV file in parallel.

    This is a simplified wrapper around runCASSIA_score_batch for backward compatibility.

    Args:
        results_file_path (str): Path to the results CSV file
        output_file_path (str, optional): Path to save the updated results
        max_workers (int): Maximum number of parallel threads
        model (str): Model to use
        provider (str): AI provider to use ('openai', 'anthropic', or 'openrouter')

    Returns:
        pd.DataFrame: Original results with added score and reasoning columns
    """
    return runCASSIA_score_batch(
        input_file=results_file_path,
        output_file=output_file_path,
        max_workers=max_workers,
        model=model,
        provider=provider,
        max_retries=0  # Original behavior had no retry logic
    )
