import pandas as pd
import json
import re
import csv
import os
import sys
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from openai import OpenAI
try:
    from CASSIA.engine.main_function_code import *
    from CASSIA.core.model_settings import ModelSettings
except ImportError:
    try:
        from .main_function_code import *
        from ..core.model_settings import ModelSettings
    except ImportError:
        from main_function_code import *
        from model_settings import ModelSettings

import requests
import threading
import numpy as np
from importlib import resources
import datetime
import shutil
import atexit
from collections import Counter

# Suppress httpx and API client logs to reduce noise during batch operations
logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("openai").setLevel(logging.WARNING)
logging.getLogger("anthropic").setLevel(logging.WARNING)

# Import CASSIA logger for actionable error messages
try:
    from CASSIA.core.logging_config import get_logger
except ImportError:
    try:
        from ..core.logging_config import get_logger
    except ImportError:
        from logging_config import get_logger

logger = get_logger(__name__)

try:
    from CASSIA.core.llm_utils import *
except ImportError:
    try:
        from ..core.llm_utils import *
    except ImportError:
        from llm_utils import *

try:
    from CASSIA.core.model_settings import resolve_model_name, get_recommended_model
except ImportError:
    try:
        from ..core.model_settings import resolve_model_name, get_recommended_model
    except ImportError:
        from model_settings import resolve_model_name, get_recommended_model

try:
    from CASSIA.core.validation import (
        validate_runCASSIA_inputs,
        validate_runCASSIA_batch_inputs,
        validate_runCASSIA_with_reference_inputs
    )
    from CASSIA.core.exceptions import CASSIAValidationError
except ImportError:
    try:
        from ..core.validation import (
            validate_runCASSIA_inputs,
            validate_runCASSIA_batch_inputs,
            validate_runCASSIA_with_reference_inputs
        )
        from ..core.exceptions import CASSIAValidationError
    except ImportError:
        from validation import (
            validate_runCASSIA_inputs,
            validate_runCASSIA_batch_inputs,
            validate_runCASSIA_with_reference_inputs
        )
        from exceptions import CASSIAValidationError

try:
    from CASSIA.agents.merging.merging_annotation import *
except ImportError:
    try:
        from ..agents.merging.merging_annotation import *
    except ImportError:
        from merging_annotation import *

try:
    from CASSIA.evaluation.cell_type_comparison import compareCelltypes
except ImportError:
    try:
        from ..evaluation.cell_type_comparison import compareCelltypes
    except ImportError:
        from cell_type_comparison import compareCelltypes

# Reference Agent for intelligent reference retrieval
try:
    from CASSIA.agents.reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt
except ImportError:
    try:
        from ..agents.reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt
    except ImportError:
        try:
            from reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt
        except ImportError:
            # Reference agent not available - provide stub
            ReferenceAgent = None
            get_reference_content = None
            format_reference_for_prompt = None

# Import from extracted modules (these were previously defined locally)
try:
    from CASSIA.core.progress_tracker import BatchProgressTracker
except ImportError:
    try:
        from ..core.progress_tracker import BatchProgressTracker
    except ImportError:
        from progress_tracker import BatchProgressTracker

try:
    from CASSIA.core.utils import safe_get, natural_sort_key, clean_conversation_history, write_csv
except ImportError:
    try:
        from ..core.utils import safe_get, natural_sort_key, clean_conversation_history, write_csv
    except ImportError:
        from utils import safe_get, natural_sort_key, clean_conversation_history, write_csv

try:
    from CASSIA.core.marker_utils import split_markers, get_top_markers, _validate_ranking_parameters, _prepare_ranking_column, _get_sort_direction
except ImportError:
    try:
        from ..core.marker_utils import split_markers, get_top_markers, _validate_ranking_parameters, _prepare_ranking_column, _get_sort_direction
    except ImportError:
        from marker_utils import split_markers, get_top_markers, _validate_ranking_parameters, _prepare_ranking_column, _get_sort_direction


def _normalize_reasoning(reasoning):
    """
    Normalize reasoning parameter to dict format.

    Accepts either a string ("high", "medium", "low") or a dict.
    Converts string to dict format for internal use.

    Args:
        reasoning: String like "high" or dict like {"effort": "high"}

    Returns:
        None, or dict like {"effort": "high"}
    """
    if reasoning is None:
        return None
    if isinstance(reasoning, str):
        return {"effort": reasoning.lower()}
    return reasoning  # Already a dict


def set_openai_api_key(api_key):
    os.environ["OPENAI_API_KEY"] = api_key

def set_anthropic_api_key(api_key):
    """Set the Anthropic API key in environment variables."""
    os.environ["ANTHROPIC_API_KEY"] = api_key

def set_openrouter_api_key(api_key):
    os.environ["OPENROUTER_API_KEY"] = api_key

def set_api_key(api_key, provider="openai"):
    """
    Set the API key for the specified provider in environment variables.

    Args:
        api_key (str): The API key to set
        provider (str): The provider to set the key for ('openai', 'anthropic', 'openrouter', or a custom base URL)

    Raises:
        ValueError: If provider is not recognized
    """
    if provider.lower() == "openai":
        os.environ["OPENAI_API_KEY"] = api_key
    elif provider.lower() == "anthropic":
        os.environ["ANTHROPIC_API_KEY"] = api_key
    elif provider.lower() == "openrouter":
        os.environ["OPENROUTER_API_KEY"] = api_key
    elif provider.lower().startswith("http"):
        os.environ["CUSTOMIZED_API_KEY"] = api_key
    else:
        raise ValueError("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)")


def check_formatted_output(structured_output):
    return 'main_cell_type' in structured_output and 'sub_cell_types' in structured_output


def rerun_formatting_agent(agent, full_conversation_history):
    full_text = "\n\n".join([f"{role}: {message}" for role, message in full_conversation_history])
    formatted_result = agent(full_text, "user")
    return extract_json_from_reply(formatted_result)


def _run_core_analysis(model, temperature, marker_list, tissue, species, additional_info, provider, validator_involvement, reasoning=None):
    """
    Internal helper that runs the actual LLM analysis.

    Args:
        reasoning: Optional reasoning configuration for models that support it.

    Returns:
        tuple: (analysis_result, conversation_history)
    """
    if provider.lower() == "openai":
        return run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info, validator_involvement, reasoning=reasoning)
    elif provider.lower() == "anthropic":
        return run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info, validator_involvement, reasoning=reasoning)
    elif provider.lower() == "openrouter":
        return run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info, validator_involvement, reasoning=reasoning)
    elif provider.lower().startswith("http"):
        api_key = os.environ.get("CUSTOMIZED_API_KEY")
        if not api_key:
            raise ValueError("CUSTOMIZED_API_KEY environment variable is not set. Please call set_api_key with your API key and provider (base URL).")
        return run_cell_type_analysis_custom(
            base_url=provider,
            api_key=api_key,
            model=model,
            temperature=temperature,
            marker_list=marker_list,
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            validator_involvement=validator_involvement,
            reasoning=reasoning
        )
    else:
        raise ValueError("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)")


def runCASSIA(
    model="google/gemini-2.5-flash-preview",
    temperature=0,
    marker_list=None,
    tissue="lung",
    species="human",
    additional_info=None,
    provider="openrouter",
    validator_involvement="v1",
    reasoning=None,
    # Reference parameters (optional, default off for backward compatibility)
    use_reference=False,
    reference_threshold=40,
    reference_provider=None,
    reference_model=None,
    skip_reference_llm=False,
    verbose=False
):
    """
    Run cell type analysis using OpenAI, Anthropic, OpenRouter, or a custom OpenAI-compatible provider.

    Optionally enhances analysis with intelligent reference retrieval based on marker complexity.

    Args:
        model (str): Model name to use
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for the analysis
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a base URL)
        validator_involvement (str): Validator involvement level ('v1', 'v0', etc.)
        reasoning (str or dict): Optional reasoning effort for models that support it.
            Controls how much the model "thinks" before responding.
            Simple: "high", "medium", or "low"
            Dict: {"effort": "high|medium|low"} (for advanced use)
            Supported models: OpenAI GPT-5 series, Anthropic Claude Opus 4.5, compatible via OpenRouter.
        use_reference (bool): Whether to use intelligent reference retrieval (default: False)
        reference_threshold (float): Complexity score threshold for triggering reference (0-100)
        reference_provider (str): Provider for reference complexity assessment (default: same as provider)
        reference_model (str): Model for reference complexity assessment (default: fast model)
        skip_reference_llm (bool): Skip LLM complexity assessment, use rules only
        verbose (bool): Print reference retrieval info (default: False)

    Returns:
        tuple: (analysis_result, conversation_history, reference_info)
            - reference_info is a dict with keys: reference_used, complexity_score,
              preliminary_cell_type, references_used, reason

    Raises:
        CASSIAValidationError: If input validation fails
    """
    # Normalize reasoning parameter (accept string or dict)
    reasoning = _normalize_reasoning(reasoning)

    # Validate all inputs early (fail-fast)
    validated = validate_runCASSIA_inputs(
        model=model,
        temperature=temperature,
        marker_list=marker_list,
        tissue=tissue,
        species=species,
        provider=provider,
        additional_info=additional_info,
        validator_involvement=validator_involvement
    )

    # Use validated values (some may be normalized)
    marker_list = validated['marker_list']
    temperature = validated['temperature']
    tissue = validated['tissue']
    species = validated['species']
    provider = validated['provider']
    model = validated['model']
    validator_involvement = validated['validator_involvement']

    # Resolve fuzzy model names to full model names (e.g., "gpt" -> "gpt-5.1")
    settings = ModelSettings()
    model, provider = settings.resolve_model_name(model, provider, verbose=False)

    # Initialize reference_info (always returned)
    reference_info = {
        "reference_used": False,
        "complexity_score": None,
        "preliminary_cell_type": None,
        "references_used": [],
        "reason": ""
    }

    # If reference is disabled, run standard analysis
    if not use_reference:
        reference_info["reason"] = "Reference disabled (use_reference=False)"
        result, history = _run_core_analysis(
            model, temperature, marker_list, tissue, species,
            additional_info, provider, validator_involvement, reasoning=reasoning
        )
        return result, history, reference_info

    # Reference is enabled - check if ReferenceAgent is available
    if ReferenceAgent is None:
        if verbose:
            logger.warning("Reference agent not available. Running without reference.")
        reference_info["reason"] = "Reference agent not available"
        result, history = _run_core_analysis(
            model, temperature, marker_list, tissue, species,
            additional_info, provider, validator_involvement, reasoning=reasoning
        )
        return result, history, reference_info

    # Initialize reference agent
    ref_provider = reference_provider or provider
    ref_agent = ReferenceAgent(provider=ref_provider, model=reference_model)

    # Get reference content
    if verbose:
        print("Analyzing markers for reference retrieval...")

    ref_result = ref_agent.get_reference_for_markers(
        markers=marker_list[:20] if marker_list else [],
        tissue=tissue,
        species=species,
        threshold=reference_threshold,
        skip_llm=skip_reference_llm
    )

    reference_info["complexity_score"] = ref_result.get("complexity_score")
    reference_info["preliminary_cell_type"] = ref_result.get("preliminary_cell_type")
    reference_info["cell_type_range"] = ref_result.get("cell_type_range", [])

    # Determine if reference should be used
    if ref_result.get("should_use_reference") and ref_result.get("content"):
        reference_info["reference_used"] = True
        reference_info["references_used"] = ref_result.get("references_used", [])
        reference_info["reason"] = ref_result.get("reasoning", "Reference retrieved")

        if verbose:
            print(f"  Complexity score: {reference_info['complexity_score']}/100")
            print(f"  Preliminary cell type: {reference_info['preliminary_cell_type']}")
            print(f"  References used: {', '.join(reference_info['references_used'])}")

        # Format reference for injection
        reference_content = format_reference_for_prompt(ref_result)

        # Combine with existing additional_info
        if additional_info:
            combined_info = f"{additional_info}\n\n{reference_content}"
        else:
            combined_info = reference_content
    else:
        reference_info["reason"] = ref_result.get("reasoning", "Reference not needed")
        combined_info = additional_info

        if verbose:
            print(f"  Reference not needed. {reference_info['reason']}")

    # Run analysis with (possibly enhanced) additional_info
    result, history = _run_core_analysis(
        model, temperature, marker_list, tissue, species,
        combined_info, provider, validator_involvement, reasoning=reasoning
    )

    return result, history, reference_info


def runCASSIA_with_reference(*args, **kwargs):
    """
    DEPRECATED: Use runCASSIA(use_reference=True) instead.

    This function is kept for backward compatibility and will be removed in a future version.
    """
    import warnings
    warnings.warn(
        "runCASSIA_with_reference is deprecated. Use runCASSIA(use_reference=True) instead.",
        DeprecationWarning,
        stacklevel=2
    )
    # Set use_reference=True and verbose=True (original defaults)
    kwargs.setdefault('use_reference', True)
    kwargs.setdefault('verbose', True)
    return runCASSIA(*args, **kwargs)


def runCASSIA_batch(
    marker,
    output_name="cell_type_analysis_results.json",
    n_genes=50,
    model="google/gemini-2.5-flash-preview",
    temperature=0,
    tissue="lung",
    species="human",
    additional_info=None,
    celltype_column=None,
    gene_column_name=None,
    max_workers=10,
    provider="openrouter",
    max_retries=1,
    ranking_method="avg_log2FC",
    ascending=None,
    validator_involvement="v1",
    reasoning=None,
    # Reference parameters (NEW)
    use_reference=False,
    reference_model=None,
    verbose=True
):
    """
    Run cell type analysis on multiple clusters in parallel.

    Optionally uses intelligent per-cluster reference retrieval based on marker complexity.

    Args:
        marker: Input marker data (pandas DataFrame or path to CSV file)
        output_name (str): Base name for output files
        n_genes (int): Number of top genes to extract per cluster
        model (str): Model name to use for analysis
        temperature (float): Temperature parameter for the model
        tissue (str): Tissue type being analyzed
        species (str): Species being analyzed
        additional_info (str): Additional information for analysis
        celltype_column (str): Column name containing cell type names (default: first column)
        gene_column_name (str): Column name containing gene markers (default: second column)
        max_workers (int): Maximum number of parallel workers
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or base URL)
        max_retries (int): Maximum number of retries for failed analyses
        ranking_method (str): Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending (bool): Sort direction (None uses default for each method)
        validator_involvement (str): Validator involvement level ('v1', 'v0', etc.)
        reasoning (str or dict): Optional reasoning effort for models that support it.
            Controls how much the model "thinks" before responding.
            Simple: "high", "medium", or "low"
            Dict: {"effort": "high|medium|low"} (for advanced use)
            Supported: OpenAI GPT-5, Anthropic Claude Opus 4.5, compatible via OpenRouter.
        use_reference (bool): Whether to use intelligent reference retrieval per cluster (default: False)
        reference_model (str): Model for reference complexity assessment (default: fast model)
        verbose (bool): Print progress information (default: True)

    Returns:
        dict: Results dictionary containing analysis results for each cell type

    Raises:
        CASSIAValidationError: If input validation fails
    """
    # Normalize reasoning parameter (accept string or dict)
    reasoning = _normalize_reasoning(reasoning)

    # Validate all inputs early (fail-fast)
    validated = validate_runCASSIA_batch_inputs(
        marker=marker,
        output_name=output_name,
        n_genes=n_genes,
        model=model,
        temperature=temperature,
        tissue=tissue,
        species=species,
        max_workers=max_workers,
        provider=provider,
        max_retries=max_retries,
        ranking_method=ranking_method,
        validator_involvement=validator_involvement
    )

    # Use validated values
    marker = validated['marker']
    temperature = validated['temperature']
    tissue = validated['tissue']
    species = validated['species']
    provider = validated['provider']
    model = validated['model']
    n_genes = validated['n_genes']
    max_workers = validated['max_workers']
    max_retries = validated['max_retries']
    ranking_method = validated['ranking_method']
    validator_involvement = validated['validator_involvement']

    # Resolve fuzzy model names ONCE before batch starts (e.g., "gpt" -> "gpt-5.1")
    settings = ModelSettings()
    model, provider = settings.resolve_model_name(model, provider, verbose=True)

    # Load the dataframe
    if isinstance(marker, pd.DataFrame):
        df = marker.copy()
    elif isinstance(marker, str):
        df = pd.read_csv(marker)
    else:
        raise ValueError("marker must be either a pandas DataFrame or a string path to a CSV file")

    # If dataframe has only two columns, assume it's already processed
    if len(df.columns) == 2:
        print("Using input dataframe directly as it appears to be pre-processed (2 columns)")
    else:
        print("Processing input dataframe to get top markers")
        df = get_top_markers(df, n_genes=n_genes, ranking_method=ranking_method, ascending=ascending)

    # If celltype_column is not specified, use the first column
    if celltype_column is None:
        celltype_column = df.columns[0]

    # If gene_column_name is not specified, use the second column
    if gene_column_name is None:
        gene_column_name = df.columns[1]


    # Check if reference agent is available when requested
    if use_reference and ReferenceAgent is None:
        if verbose:
            print("Warning: Reference agent not available. Running without reference retrieval.")
        use_reference = False

    # Initialize progress tracker
    total_clusters = len(df)
    tracker = BatchProgressTracker(total_clusters)

    ref_status = "with reference retrieval" if use_reference else ""
    if verbose:
        print(f"\nStarting analysis of {total_clusters} clusters {ref_status} with {max_workers} parallel workers...\n")

    # Track reference usage statistics
    reference_stats = {'used': 0, 'not_needed': 0, 'errors': 0}
    stats_lock = threading.Lock()

    def analyze_cell_type(cell_type, marker_list):
        tracker.start_task(cell_type)
        for attempt in range(max_retries + 1):
            try:
                result, conversation_history, ref_info = runCASSIA(
                    model=model,
                    temperature=temperature,
                    marker_list=marker_list,
                    tissue=tissue,
                    species=species,
                    additional_info=additional_info,
                    provider=provider,
                    validator_involvement=validator_involvement,
                    reasoning=reasoning,
                    use_reference=use_reference,
                    reference_model=reference_model,
                    verbose=False  # Suppress per-cluster verbose output
                )
                # Add the number of markers and marker list to the result
                result['num_markers'] = len(marker_list)
                result['marker_list'] = marker_list
                result['reference_info'] = ref_info

                # Track reference statistics
                with stats_lock:
                    if ref_info.get('reference_used'):
                        reference_stats['used'] += 1
                    else:
                        reference_stats['not_needed'] += 1

                tracker.complete_task(cell_type)
                return cell_type, result, conversation_history, ref_info
            except Exception as exc:
                # Don't retry authentication errors
                if "401" in str(exc) or "API key" in str(exc) or "authentication" in str(exc).lower():
                    tracker.complete_task(cell_type)  # Remove from active even on failure
                    logger.error(
                        f"Authentication failed for cluster '{cell_type}'. "
                        f"Please check your API key is valid. "
                        f"Set it with: CASSIA.set_api_key(provider, 'your-key'). "
                        f"Error: {exc}"
                    )
                    raise exc

                # For other errors, retry if attempts remain
                if attempt < max_retries:
                    wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s...
                    logger.warning(
                        f"API call failed for cluster '{cell_type}' (attempt {attempt+1}/{max_retries+1}). "
                        f"Retrying in {wait_time}s... "
                        f"If this keeps happening, check your API key or try a different model. "
                        f"Error: {exc}"
                    )
                    time.sleep(wait_time)
                else:
                    tracker.complete_task(cell_type)  # Remove from active on final failure
                    with stats_lock:
                        reference_stats['errors'] += 1
                    logger.error(
                        f"All {max_retries+1} attempts failed for cluster '{cell_type}'. "
                        f"Consider using a different model or reducing max_workers. "
                        f"Error: {exc}"
                    )
                    raise exc

    results = {}
    failed_analyses = []  # Track failed cell types for reporting

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_celltype = {executor.submit(analyze_cell_type, row[celltype_column], split_markers(row[gene_column_name])): row[celltype_column] for _, row in df.iterrows()}

        # Process completed tasks
        for future in as_completed(future_to_celltype):
            cell_type = future_to_celltype[future]
            try:
                cell_type, result, conversation_history, ref_info = future.result()
                if result:
                    results[cell_type] = {
                        "analysis_result": result,
                        "conversation_history": conversation_history,
                        "iterations": result.get("iterations", 1),
                        "reference_info": ref_info
                    }
            except Exception as exc:
                failed_analyses.append((cell_type, str(exc)))

    # Finalize progress display
    tracker.finish()

    # Report any failures
    if failed_analyses:
        print(f"\nWarning: {len(failed_analyses)} analysis/analyses failed:")
        for cell_type, error in failed_analyses:
            print(f"  - {cell_type}: {error[:100]}{'...' if len(error) > 100 else ''}")
        print()

    # Report reference statistics if reference was used
    if use_reference and verbose:
        print(f"\nReference retrieval statistics:")
        print(f"  - Used: {reference_stats['used']}")
        print(f"  - Not needed: {reference_stats['not_needed']}")
        if reference_stats['errors'] > 0:
            print(f"  - Errors: {reference_stats['errors']}")
        print()

    print(f"All analyses completed. Results saved to '{output_name}'.")

    # Prepare data for both CSV files and HTML report

    full_data = []
    full_data_for_html = []  # Keep raw conversation history with newlines for HTML
    summary_data = []

    for true_cell_type, details in results.items():
        main_cell_type = safe_get(details, 'analysis_result', 'main_cell_type')
        sub_cell_types = ', '.join(safe_get(details, 'analysis_result', 'sub_cell_types') or [])
        possible_mixed_cell_types = ', '.join(safe_get(details, 'analysis_result', 'possible_mixed_cell_types') or [])
        marker_number = safe_get(details, 'analysis_result', 'num_markers')
        marker_list = ', '.join(safe_get(details, 'analysis_result', 'marker_list') or [])
        iterations = safe_get(details, 'analysis_result', 'iterations')

        # Extract reference info if available
        ref_info = details.get('reference_info', {})
        ref_used = "Yes" if ref_info.get('reference_used') else "No"
        complexity_score = ref_info.get('complexity_score', '')

        # Process conversation history - keep raw version for HTML
        raw_conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in safe_get(details, 'conversation_history') or []])
        conversation_history = clean_conversation_history(raw_conversation_history)  # Clean for CSV storage

        # Data for HTML report (with raw conversation history preserving newlines)
        html_row = {
            'Cluster ID': true_cell_type,
            'Predicted General Cell Type': main_cell_type,
            'Predicted Detailed Cell Type': sub_cell_types,
            'Possible Mixed Cell Types': possible_mixed_cell_types,
            'Marker Number': marker_number,
            'Marker List': marker_list,
            'Iterations': iterations,
            'Model': model,
            'Provider': provider,
            'Tissue': tissue,
            'Species': species,
            'Additional Info': additional_info or "N/A",
            'Conversation History': raw_conversation_history  # Raw with newlines
        }
        if use_reference:
            html_row['Reference Used'] = ref_used
            html_row['Complexity Score'] = complexity_score
        full_data_for_html.append(html_row)

        # Build full data row
        full_row = [
            true_cell_type,
            main_cell_type,
            sub_cell_types,
            possible_mixed_cell_types,
            marker_number,
            marker_list,
            iterations,
            model,
            provider,
            tissue,
            species,
            additional_info or "N/A",
            conversation_history  # Cleaned for CSV
        ]
        if use_reference:
            full_row.extend([ref_used, complexity_score])
        full_data.append(full_row)

        # Build summary data row
        summary_row = [
            true_cell_type,
            main_cell_type,
            sub_cell_types,
            possible_mixed_cell_types,
            marker_list,
            iterations,
            model,
            provider,
            tissue,
            species
        ]
        if use_reference:
            summary_row.append(ref_used)
        summary_data.append(summary_row)

    # Generate output filenames based on input JSON filename
    base_name = os.path.splitext(output_name)[0]
    full_csv_name = f"{base_name}_full.csv"
    summary_csv_name = f"{base_name}_summary.csv"

    # Make sure the output directory exists
    output_dir = os.path.dirname(full_csv_name)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Sort the data by Cluster ID with natural/numeric ordering
    # This ensures "cluster 1", "cluster 2", "cluster 10" (not "cluster 1", "cluster 10", "cluster 2")
    full_data.sort(key=lambda x: natural_sort_key(x[0]))  # Sort by first column (Cluster ID) numerically
    summary_data.sort(key=lambda x: natural_sort_key(x[0]))  # Sort by first column (Cluster ID) numerically

    # Build headers (with optional reference columns)
    full_headers = ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
                    'Possible Mixed Cell Types', 'Marker Number', 'Marker List', 'Iterations',
                    'Model', 'Provider', 'Tissue', 'Species', 'Additional Info', 'Conversation History']
    summary_headers = ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
                       'Possible Mixed Cell Types', 'Marker List', 'Iterations', 'Model', 'Provider',
                       'Tissue', 'Species']
    if use_reference:
        full_headers.extend(['Reference Used', 'Complexity Score'])
        summary_headers.append('Reference Used')

    # Write the full data CSV
    write_csv(full_csv_name, full_headers, full_data)

    # Write the summary data CSV
    write_csv(summary_csv_name, summary_headers, summary_data)

    # Generate HTML report with raw conversation history (preserving newlines)
    html_report_name = f"{base_name}_report.html"
    try:
        # Try relative import first (when used as package), fall back to absolute
        try:
            from CASSIA.reports.generate_batch_report import generate_batch_html_report_from_data
        except ImportError:
            try:
                from ..reports.generate_batch_report import generate_batch_html_report_from_data
            except ImportError:
                from generate_batch_report import generate_batch_html_report_from_data
        # Sort the HTML data the same way as CSV
        full_data_for_html.sort(key=lambda x: natural_sort_key(x['Cluster ID']))
        generate_batch_html_report_from_data(full_data_for_html, html_report_name)
        print(f"Three files have been created:")
        print(f"1. {full_csv_name} (full data CSV)")
        print(f"2. {summary_csv_name} (summary data CSV)")
        print(f"3. {html_report_name} (interactive HTML report)")
    except Exception as e:
        print(f"Warning: Could not generate HTML report: {e}")
        print(f"Two CSV files have been created:")
        print(f"1. {full_csv_name} (full data)")
        print(f"2. {summary_csv_name} (summary data)")


def runCASSIA_batch_with_reference(*args, **kwargs):
    """
    DEPRECATED: Use runCASSIA_batch(use_reference=True) instead.

    This function is kept for backward compatibility and will be removed in a future version.
    """
    import warnings
    warnings.warn(
        "runCASSIA_batch_with_reference is deprecated. Use runCASSIA_batch(use_reference=True) instead.",
        DeprecationWarning,
        stacklevel=2
    )
    # Set use_reference=True (original default behavior when this function was called)
    kwargs.setdefault('use_reference', True)
    return runCASSIA_batch(*args, **kwargs)


# =============================================================================
# RE-EXPORTS FOR BACKWARD COMPATIBILITY
# =============================================================================
# The functions below have been moved to separate modules but are re-exported
# here to maintain backward compatibility with existing code that imports from
# tools_function.py (including R package via reticulate).

# From annotation_boost.py
try:
    from CASSIA.agents.annotation_boost.annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task
except ImportError:
    try:
        from ..agents.annotation_boost.annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task
    except ImportError:
        from annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task

# From generate_reports.py (HTML report generation)
try:
    from CASSIA.reports.generate_reports import (
        generate_analysis_html_report as generate_html_report,  # Renamed for backward compat
        process_single_report,
        generate_score_index_page as generate_index_page,  # Renamed for backward compat
        runCASSIA_generate_score_report
    )
except ImportError:
    try:
        from ..reports.generate_reports import (
            generate_analysis_html_report as generate_html_report,
            process_single_report,
            generate_score_index_page as generate_index_page,
            runCASSIA_generate_score_report
        )
    except ImportError:
        from generate_reports import (
            generate_analysis_html_report as generate_html_report,
            process_single_report,
            generate_score_index_page as generate_index_page,
            runCASSIA_generate_score_report
        )

# From pipeline.py (main pipeline)
try:
    from CASSIA.pipeline.pipeline import runCASSIA_pipeline
except ImportError:
    try:
        from ..pipeline.pipeline import runCASSIA_pipeline
    except ImportError:
        from pipeline import runCASSIA_pipeline

# From marker_utils.py (marker loading utilities)
try:
    from CASSIA.core.marker_utils import loadmarker, list_available_markers
except ImportError:
    try:
        from ..core.marker_utils import loadmarker, list_available_markers
    except ImportError:
        from marker_utils import loadmarker, list_available_markers

# From scoring.py (scoring functions)
try:
    from CASSIA.evaluation.scoring import (
        prompt_creator_score,
        extract_score_and_reasoning,
        score_single_analysis,
        process_single_row,
        runCASSIA_score_batch,
        score_annotation_batch
    )
except ImportError:
    try:
        from ..evaluation.scoring import (
            prompt_creator_score,
            extract_score_and_reasoning,
            score_single_analysis,
            process_single_row,
            runCASSIA_score_batch,
            score_annotation_batch
        )
    except ImportError:
        from scoring import (
            prompt_creator_score,
            extract_score_and_reasoning,
            score_single_analysis,
            process_single_row,
            runCASSIA_score_batch,
            score_annotation_batch
        )

# From Uncertainty_quantification.py (required by R package for uncertainty functions)
try:
    from CASSIA.agents.uncertainty.Uncertainty_quantification import (
        runCASSIA_n_times,
        runCASSIA_batch_n_times,
        runCASSIA_similarity_score_batch,
        runCASSIA_n_times_similarity_score
    )
except ImportError:
    try:
        from ..agents.uncertainty.Uncertainty_quantification import (
            runCASSIA_n_times,
            runCASSIA_batch_n_times,
            runCASSIA_similarity_score_batch,
            runCASSIA_n_times_similarity_score
        )
    except ImportError:
        try:
            from Uncertainty_quantification import (
                runCASSIA_n_times,
                runCASSIA_batch_n_times,
                runCASSIA_similarity_score_batch,
                runCASSIA_n_times_similarity_score
            )
        except ImportError:
            # Uncertainty_quantification module may not be available in all deployments
            pass

# End of tools_function.py
