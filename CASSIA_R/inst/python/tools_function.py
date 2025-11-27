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
    from .main_function_code import *
    from .model_settings import ModelSettings
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
    from .logging_config import get_logger
except ImportError:
    from logging_config import get_logger

logger = get_logger(__name__)

try:
    from .llm_utils import *
except ImportError:
    from llm_utils import *

try:
    from .model_settings import resolve_model_name, get_recommended_model
except ImportError:
    from model_settings import resolve_model_name, get_recommended_model

try:
    from .validation import (
        validate_runCASSIA_inputs,
        validate_runCASSIA_batch_inputs,
        validate_runCASSIA_with_reference_inputs
    )
    from .exceptions import CASSIAValidationError
except ImportError:
    from validation import (
        validate_runCASSIA_inputs,
        validate_runCASSIA_batch_inputs,
        validate_runCASSIA_with_reference_inputs
    )
    from exceptions import CASSIAValidationError

try:
    from .merging_annotation import *
except ImportError:
    from merging_annotation import *

try:
    from .cell_type_comparison import compareCelltypes
except ImportError:
    from cell_type_comparison import compareCelltypes

# Reference Agent for intelligent reference retrieval
try:
    from .reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt
except ImportError:
    try:
        from reference_agent import ReferenceAgent, get_reference_content, format_reference_for_prompt
    except ImportError:
        # Reference agent not available - provide stub
        ReferenceAgent = None
        get_reference_content = None
        format_reference_for_prompt = None


class BatchProgressTracker:
    """Thread-safe progress tracker for batch processing with visual progress bar."""

    # Spinner animation frames
    SPINNER_FRAMES = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']

    def __init__(self, total, bar_width=40, refresh_rate=0.1):
        self.total = total
        self.completed = 0
        self.in_progress = set()
        self.lock = threading.Lock()
        self.bar_width = bar_width
        self._lines_printed = 0
        self._spinner_idx = 0
        self._running = True
        self._refresh_rate = refresh_rate

        # Start background thread for continuous spinner animation
        self._animation_thread = threading.Thread(target=self._animate, daemon=True)
        self._animation_thread.start()

        # Hide cursor during animation to prevent flashing
        sys.stdout.write('\033[?25l')
        sys.stdout.flush()

        # Register cleanup to ensure cursor is restored on unexpected exit
        atexit.register(self._restore_cursor)

    def _animate(self):
        """Background thread that continuously updates the spinner."""
        while self._running:
            time.sleep(self._refresh_rate)
            with self.lock:
                if self._running and len(self.in_progress) > 0:
                    self._render()

    def _restore_cursor(self):
        """Restore cursor visibility. Called on exit or finish."""
        sys.stdout.write('\033[?25h')
        sys.stdout.flush()

    def start_task(self, name):
        """Mark a task as started/in-progress."""
        with self.lock:
            self.in_progress.add(name)
            self._render()

    def complete_task(self, name):
        """Mark a task as completed."""
        with self.lock:
            self.in_progress.discard(name)
            self.completed += 1
            self._render()

    def _render(self):
        """Render the progress display, updating in place."""
        # Calculate progress percentage
        pct = self.completed / self.total if self.total > 0 else 0
        filled = int(self.bar_width * pct)
        bar = '█' * filled + '░' * (self.bar_width - filled)

        # Calculate counts
        processing = len(self.in_progress)
        pending = self.total - self.completed - processing

        # Get spinner character (only animate when processing)
        if processing > 0:
            spinner = self.SPINNER_FRAMES[self._spinner_idx % len(self.SPINNER_FRAMES)]
            self._spinner_idx += 1
        else:
            spinner = '✓' if self.completed == self.total else '○'

        # Truncate active task names if too many
        active_names = list(self.in_progress)[:3]
        active_str = ', '.join(str(name) for name in active_names)
        if len(self.in_progress) > 3:
            active_str += f', ... (+{len(self.in_progress)-3} more)'

        # Build display lines
        lines = [
            f"CASSIA Batch Analysis {spinner}",
            f"[{bar}] {pct*100:.0f}%",
            f"Completed: {self.completed} | Processing: {processing} | Pending: {pending}",
            f"Active: {active_str if active_str else 'None'}"
        ]

        # Move cursor up to overwrite previous output
        if self._lines_printed > 0:
            sys.stdout.write(f'\033[{self._lines_printed}A')

        # Print each line, clearing to end of line
        for line in lines:
            sys.stdout.write(f'\033[K{line}\n')

        sys.stdout.flush()
        self._lines_printed = len(lines)

    def finish(self):
        """Finalize the progress display."""
        # Stop the animation thread
        self._running = False
        self._animation_thread.join(timeout=0.5)

        with self.lock:
            # Final render to show 100% with checkmark
            self._render()
            # Add blank line after completion
            sys.stdout.write('\033[K\n')
            sys.stdout.flush()

        # Restore cursor visibility
        self._restore_cursor()

        # Unregister atexit handler since we've cleaned up normally
        try:
            atexit.unregister(self._restore_cursor)
        except Exception:
            pass  # Ignore if already unregistered


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
    
def split_markers(marker_string):
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
    """Validate ranking method and column existence"""
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
    """Prepare the ranking column, calculating if necessary"""
    df_copy = df.copy()
    
    if ranking_method == "pct_diff":
        if "pct.1" not in df_copy.columns or "pct.2" not in df_copy.columns:
            raise ValueError("Columns 'pct.1' and 'pct.2' required for pct_diff ranking")
        df_copy["pct_diff"] = df_copy["pct.1"] - df_copy["pct.2"]
        return df_copy, "pct_diff"
    
    return df_copy, ranking_method


def _get_sort_direction(ranking_method, ascending):
    """Get the sort direction for the ranking method"""
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
        df: Either a pandas DataFrame (Seurat format) or dictionary (Scanpy format)
        n_genes: Number of top genes to return per cluster
        format_type: Either 'seurat', 'scanpy', or None (auto-detect)
        ranking_method: Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending: Sort direction (None uses default for each method)
    
    Returns:
        pandas DataFrame with cluster and marker columns
    """
    # Auto-detect format if not specified
    if format_type is None:
        if 'names' in df and 'scores' in df:
            format_type = 'scanpy'
        else:
            format_type = 'seurat'
    
    if format_type == 'scanpy':
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
        # Convert string \'inf\' and \'-inf\' to numeric values first
        df['avg_log2FC'] = pd.to_numeric(df['avg_log2FC'].replace({'inf': np.inf, '-inf': -np.inf}), errors='coerce')
        
        # Replace inf and -inf values with max and min non-inf values respectively
        max_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).max()
        min_non_inf = df['avg_log2FC'].replace([np.inf, -np.inf], np.nan).min()
        df['avg_log2FC'] = df['avg_log2FC'].replace([np.inf, -np.inf], [max_non_inf, min_non_inf])
        
        # Filter by adjusted p-value, positive log2FC, and PCT
        df_filtered = df[
            (df['p_val_adj'] < 0.05) & 
            (df['avg_log2FC'] > 0.25) &
            ((df['pct.1'] >=0.1) | (df['pct.2'] >=0.1))  # Add PCT filter
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
                     .groupby('cluster',observed=True)
                     .agg({'gene': lambda x: ','.join(x)})
                     .rename(columns={'gene': 'markers'})
                     .reset_index())
            
            return result
        else:
            return pd.DataFrame(columns=['cluster', 'markers'])
        

def check_formatted_output(structured_output):
    return 'main_cell_type' in structured_output and 'sub_cell_types' in structured_output


def rerun_formatting_agent(agent, full_conversation_history):
    full_text = "\n\n".join([f"{role}: {message}" for role, message in full_conversation_history])
    formatted_result = agent(full_text, "user")
    return extract_json_from_reply(formatted_result)


def safe_get(dict_obj, *keys):
    """Safely get nested dictionary values"""
    for key in keys:
        if isinstance(dict_obj, dict) and key in dict_obj:
            dict_obj = dict_obj[key]
        else:
            return None
    return dict_obj


def natural_sort_key(cell_type):
    """
    Create a sort key that handles numeric cluster names properly.

    Handles various formats:
    - Pure numbers: "0", "1", "10" → sorted as integers 0, 1, 10 (priority 0)
    - "cluster X": "cluster 0", "cluster 1", "cluster 10" → sorted by X numerically (priority 1)
    - "Cluster X": case-insensitive
    - Other text: sorted alphabetically (priority 2)

    Args:
        cell_type (str): The cell type or cluster name

    Returns:
        tuple: (sort_priority, numeric_value, string_value) for proper sorting
    """
    import re

    if not cell_type or not isinstance(cell_type, str):
        return (3, 0, str(cell_type))  # Non-string values go last

    cell_type_str = str(cell_type).strip()

    # Try to parse as pure integer
    try:
        return (0, int(cell_type_str), "")  # Pure numbers have priority 0
    except ValueError:
        pass

    # Try to extract number from "cluster X" or "Cluster X" format (case-insensitive)
    cluster_match = re.match(r'^cluster\s+(\d+)$', cell_type_str, re.IGNORECASE)
    if cluster_match:
        cluster_num = int(cluster_match.group(1))
        return (1, cluster_num, "")  # Cluster numbers have priority 1 (after pure numbers)

    # For any other text, sort alphabetically (priority 2)
    return (2, 0, cell_type_str.lower())


def clean_conversation_history(history_text):
    """
    Clean conversation history for safe CSV storage while preserving full content.

    Args:
        history_text (str): Raw conversation history text

    Returns:
        str: Cleaned text safe for CSV storage (no truncation)
    """
    if not history_text:
        return ""

    # Replace newlines with spaces (prevents row breaks in CSV/Excel)
    cleaned = history_text.replace('\n', ' ').replace('\r', ' ')

    # Collapse multiple spaces into single spaces
    cleaned = ' '.join(cleaned.split())

    # Double quotes will be handled by csv.writer's automatic escaping
    # No need to replace quotes - csv module handles this correctly

    # Return full content without truncation
    return cleaned


def write_csv(filename, headers, row_data):
    import os
    
    # Make sure the directory exists
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(row_data)


def runCASSIA(model="google/gemini-2.5-flash-preview", temperature=0, marker_list=None, tissue="lung", species="human", additional_info=None, provider="openrouter", validator_involvement="v1"):
    """
    Wrapper function to run cell type analysis using OpenAI, Anthropic, OpenRouter, or a custom OpenAI-compatible provider.

    Args:
        model (str): Model name to use
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for the analysis
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or a base URL)

    Returns:
        tuple: (analysis_result, conversation_history)

    Raises:
        CASSIAValidationError: If input validation fails
    """
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
    # Use verbose=False since this may be called from batch (which handles its own verbose message)
    settings = ModelSettings()
    model, provider = settings.resolve_model_name(model, provider, verbose=False)

    if provider.lower() == "openai":
        return run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
    elif provider.lower() == "anthropic":
        return run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
    elif provider.lower() == "openrouter":
        return run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info, validator_involvement)
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
            validator_involvement=validator_involvement
        )
    else:
        raise ValueError("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)")


def runCASSIA_with_reference(
    model="google/gemini-2.5-flash-preview",
    temperature=0,
    marker_list=None,
    tissue="lung",
    species="human",
    additional_info=None,
    provider="openrouter",
    validator_involvement="v1",
    use_reference=True,
    reference_threshold=40,
    reference_provider=None,
    reference_model=None,
    skip_reference_llm=False,
    verbose=True
):
    """
    Run CASSIA cell type analysis with intelligent reference retrieval.

    This function enhances the standard runCASSIA by automatically retrieving
    relevant expert-curated reference documents based on marker complexity.

    Args:
        model (str): Model name for annotation
        temperature (float): Temperature parameter for the model
        marker_list (list): List of markers to analyze
        tissue (str): Tissue type
        species (str): Species type
        additional_info (str): Additional information for the analysis
        provider (str): AI provider for annotation ('openai', 'anthropic', 'openrouter', or base URL)
        validator_involvement (str): Validator involvement level
        use_reference (bool): Whether to use reference retrieval (default True)
        reference_threshold (float): Complexity score threshold for triggering reference (0-100)
        reference_provider (str): Provider for reference complexity assessment (default: same as provider)
        reference_model (str): Model for reference complexity assessment (default: fast model)
        skip_reference_llm (bool): Skip LLM complexity assessment, use rules only
        verbose (bool): Print reference retrieval info

    Returns:
        tuple: (analysis_result, conversation_history, reference_info)
            - reference_info contains details about reference retrieval

    Raises:
        CASSIAValidationError: If input validation fails
    """
    # Validate all inputs early (fail-fast)
    validated = validate_runCASSIA_with_reference_inputs(
        model=model,
        temperature=temperature,
        marker_list=marker_list,
        tissue=tissue,
        species=species,
        provider=provider,
        additional_info=additional_info,
        validator_involvement=validator_involvement,
        reference_threshold=reference_threshold
    )

    # Use validated values (some may be normalized)
    marker_list = validated['marker_list']
    temperature = validated['temperature']
    tissue = validated['tissue']
    species = validated['species']
    provider = validated['provider']
    model = validated['model']
    validator_involvement = validated['validator_involvement']
    if 'reference_threshold' in validated:
        reference_threshold = validated['reference_threshold']

    if ReferenceAgent is None:
        if verbose:
            print("Warning: Reference agent not available. Running standard CASSIA.")
        result, history = runCASSIA(
            model=model,
            temperature=temperature,
            marker_list=marker_list,
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            provider=provider,
            validator_involvement=validator_involvement
        )
        return result, history, {"reference_used": False, "reason": "Reference agent not available"}

    reference_info = {
        "reference_used": False,
        "complexity_score": None,
        "preliminary_cell_type": None,
        "references_used": [],
        "reason": ""
    }

    if not use_reference:
        reference_info["reason"] = "Reference disabled by user"
        result, history = runCASSIA(
            model=model,
            temperature=temperature,
            marker_list=marker_list,
            tissue=tissue,
            species=species,
            additional_info=additional_info,
            provider=provider,
            validator_involvement=validator_involvement
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

    # Run CASSIA with (possibly enhanced) additional_info
    result, history = runCASSIA(
        model=model,
        temperature=temperature,
        marker_list=marker_list,
        tissue=tissue,
        species=species,
        additional_info=combined_info,
        provider=provider,
        validator_involvement=validator_involvement
    )

    return result, history, reference_info



def runCASSIA_batch(marker, output_name="cell_type_analysis_results.json", n_genes=50, model="google/gemini-2.5-flash-preview", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column=None, gene_column_name=None, max_workers=10, provider="openrouter", max_retries=1, ranking_method="avg_log2FC", ascending=None, validator_involvement="v1"):
    """
    Run cell type analysis on multiple clusters in parallel.

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

    Returns:
        dict: Results dictionary containing analysis results for each cell type

    Raises:
        CASSIAValidationError: If input validation fails
    """
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
    
    
    # Initialize progress tracker
    total_clusters = len(df)
    tracker = BatchProgressTracker(total_clusters)
    print(f"\nStarting analysis of {total_clusters} clusters with {max_workers} parallel workers...\n")

    def analyze_cell_type(cell_type, marker_list):
        tracker.start_task(cell_type)
        for attempt in range(max_retries + 1):
            try:
                result, conversation_history = runCASSIA(
                    model=model,
                    temperature=temperature,
                    marker_list=marker_list,
                    tissue=tissue,
                    species=species,
                    additional_info=additional_info,
                    provider=provider,
                    validator_involvement=validator_involvement
                )
                # Add the number of markers and marker list to the result
                result['num_markers'] = len(marker_list)
                result['marker_list'] = marker_list
                tracker.complete_task(cell_type)
                return cell_type, result, conversation_history
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
                cell_type, result, conversation_history = future.result()
                if result:
                    results[cell_type] = {
                        "analysis_result": result,
                        "conversation_history": conversation_history,
                        "iterations": result.get("iterations", 1)
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

        # Process conversation history - keep raw version for HTML
        raw_conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in safe_get(details, 'conversation_history') or []])
        conversation_history = clean_conversation_history(raw_conversation_history)  # Clean for CSV storage

        # Data for HTML report (with raw conversation history preserving newlines)
        full_data_for_html.append({
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
        })

        full_data.append([
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
        ])
        summary_data.append([
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
        ])

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

    # Write the full data CSV with updated headers
    write_csv(full_csv_name,
              ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
               'Possible Mixed Cell Types', 'Marker Number', 'Marker List', 'Iterations',
               'Model', 'Provider', 'Tissue', 'Species', 'Additional Info', 'Conversation History'],
              full_data)

    # Write the summary data CSV with updated headers
    write_csv(summary_csv_name,
              ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
               'Possible Mixed Cell Types', 'Marker List', 'Iterations', 'Model', 'Provider',
               'Tissue', 'Species'],
              summary_data)

    # Generate HTML report with raw conversation history (preserving newlines)
    html_report_name = f"{base_name}_report.html"
    try:
        # Try relative import first (when used as package), fall back to absolute
        try:
            from .generate_batch_report import generate_batch_html_report_from_data
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


def runCASSIA_batch_with_reference(
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
    use_reference=False,
    reference_model=None,
    verbose=True
):
    """
    Run cell type analysis on multiple clusters with intelligent per-cluster reference retrieval.

    Uses a TWO-STEP ReAct workflow for each cluster:
        - Step 1: LLM assesses marker complexity and decides if reference is needed
        - Step 2: If needed, LLM sees router and selects specific references
        - Selected references are added to that cluster's additional_info

    Args:
        marker: Input marker data (pandas DataFrame or path to CSV file)
        output_name (str): Base name for output files
        n_genes (int): Number of top genes to extract per cluster
        model (str): Model name to use for analysis
        temperature (float): Temperature parameter for the model
        tissue (str): Tissue type being analyzed
        species (str): Species being analyzed
        additional_info (str): Base additional information (shared across all clusters)
        celltype_column (str): Column name containing cell type names (default: first column)
        gene_column_name (str): Column name containing gene markers (default: second column)
        max_workers (int): Maximum number of parallel workers
        provider (str): AI provider to use ('openai', 'anthropic', 'openrouter', or base URL)
        max_retries (int): Maximum number of retries for failed analyses
        ranking_method (str): Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score')
        ascending (bool): Sort direction (None uses default for each method)
        validator_involvement (str): Validator mode ('v1', 'v0', etc.)
        use_reference (bool): If True, use two-step reference retrieval per cluster
        reference_model (str): Model to use for reference agent LLM calls (defaults to fast model)
        verbose (bool): Print reference retrieval progress

    Returns:
        dict: Results dictionary containing analysis results for each cell type
    """
    # Check if reference agent is available
    if use_reference and ReferenceAgent is None:
        print("Warning: Reference agent not available. Running without reference retrieval.")
        use_reference = False

    # Initialize reference agent if needed
    ref_agent = None
    if use_reference:
        ref_agent = ReferenceAgent(
            provider=provider,
            model=reference_model
        )
        if verbose:
            print("Reference agent initialized (two-step ReAct workflow enabled)")

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

    # Initialize progress tracker
    total_clusters = len(df)
    tracker = BatchProgressTracker(total_clusters)

    ref_status = "with reference retrieval" if use_reference else "without reference retrieval"
    print(f"\nStarting analysis of {total_clusters} clusters {ref_status} with {max_workers} parallel workers...\n")

    # Track reference usage statistics
    reference_stats = {'used': 0, 'not_needed': 0, 'errors': 0}

    def analyze_cell_type_with_reference(cell_type, marker_list):
        """Inner function to analyze a single cluster with optional reference retrieval."""
        tracker.start_task(cell_type)

        # Start with base additional_info
        cluster_additional_info = additional_info or ""
        ref_info = None

        # Get reference for this cluster if enabled
        if use_reference and ref_agent is not None:
            try:
                ref_result = ref_agent.get_reference_for_markers(
                    markers=marker_list[:20],
                    tissue=tissue,
                    species=species
                )

                if ref_result.get('should_use_reference') and ref_result.get('content'):
                    # Combine base additional_info with reference content
                    ref_content = format_reference_for_prompt(ref_result)
                    if cluster_additional_info:
                        cluster_additional_info = f"{cluster_additional_info}\n\n{ref_content}"
                    else:
                        cluster_additional_info = ref_content

                    reference_stats['used'] += 1
                    ref_info = {
                        'preliminary_cell_type': ref_result.get('preliminary_cell_type'),
                        'complexity_score': ref_result.get('complexity_score'),
                        'references_used': ref_result.get('references_used', [])
                    }
                else:
                    reference_stats['not_needed'] += 1
                    ref_info = {
                        'preliminary_cell_type': ref_result.get('preliminary_cell_type'),
                        'complexity_score': ref_result.get('complexity_score'),
                        'references_used': [],
                        'reason': 'LLM determined reference not needed'
                    }
            except Exception as ref_exc:
                reference_stats['errors'] += 1
                if verbose:
                    print(f"  Warning: Reference retrieval failed for {cell_type}: {str(ref_exc)[:50]}")

        # Run CASSIA with per-cluster additional_info
        for attempt in range(max_retries + 1):
            try:
                result, conversation_history = runCASSIA(
                    model=model,
                    temperature=temperature,
                    marker_list=marker_list,
                    tissue=tissue,
                    species=species,
                    additional_info=cluster_additional_info if cluster_additional_info else None,
                    provider=provider,
                    validator_involvement=validator_involvement
                )
                # Add metadata to result
                result['num_markers'] = len(marker_list)
                result['marker_list'] = marker_list
                if ref_info:
                    result['reference_info'] = ref_info

                tracker.complete_task(cell_type)
                return cell_type, result, conversation_history
            except Exception as exc:
                # Don't retry authentication errors
                if "401" in str(exc) or "API key" in str(exc) or "authentication" in str(exc).lower():
                    tracker.complete_task(cell_type)
                    raise exc

                # For other errors, retry if attempts remain
                if attempt < max_retries:
                    pass  # Continue to next attempt
                else:
                    tracker.complete_task(cell_type)
                    raise exc

    results = {}
    failed_analyses = []

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_celltype = {
            executor.submit(
                analyze_cell_type_with_reference,
                row[celltype_column],
                split_markers(row[gene_column_name])
            ): row[celltype_column] for _, row in df.iterrows()
        }

        # Process completed tasks
        for future in as_completed(future_to_celltype):
            cell_type = future_to_celltype[future]
            try:
                cell_type, result, conversation_history = future.result()
                if result:
                    results[cell_type] = {
                        "analysis_result": result,
                        "conversation_history": conversation_history,
                        "iterations": result.get("iterations", 1)
                    }
            except Exception as exc:
                failed_analyses.append((cell_type, str(exc)))

    # Finalize progress display
    tracker.finish()

    # Report reference usage statistics
    if use_reference and verbose:
        print(f"\nReference retrieval statistics:")
        print(f"  - References used: {reference_stats['used']}")
        print(f"  - Not needed (LLM decision): {reference_stats['not_needed']}")
        if reference_stats['errors'] > 0:
            print(f"  - Errors: {reference_stats['errors']}")

    # Report any failures
    if failed_analyses:
        print(f"\nWarning: {len(failed_analyses)} analysis/analyses failed:")
        for cell_type, error in failed_analyses:
            print(f"  - {cell_type}: {error[:100]}{'...' if len(error) > 100 else ''}")
        print()

    print(f"All analyses completed. Results saved to '{output_name}'.")

    # Prepare data for both CSV files and HTML report
    full_data = []
    full_data_for_html = []
    summary_data = []

    for true_cell_type, details in results.items():
        main_cell_type = safe_get(details, 'analysis_result', 'main_cell_type')
        sub_cell_types = ', '.join(safe_get(details, 'analysis_result', 'sub_cell_types') or [])
        possible_mixed_cell_types = ', '.join(safe_get(details, 'analysis_result', 'possible_mixed_cell_types') or [])
        marker_number = safe_get(details, 'analysis_result', 'num_markers')
        marker_list_str = ', '.join(safe_get(details, 'analysis_result', 'marker_list') or [])
        iterations = safe_get(details, 'analysis_result', 'iterations')

        # Reference info
        ref_info = safe_get(details, 'analysis_result', 'reference_info')
        ref_used = 'Yes' if ref_info and ref_info.get('references_used') else 'No'
        ref_complexity = ref_info.get('complexity_score', 'N/A') if ref_info else 'N/A'

        # Process conversation history
        raw_conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in safe_get(details, 'conversation_history') or []])
        conversation_history = clean_conversation_history(raw_conversation_history)

        # Data for HTML report
        full_data_for_html.append({
            'Cluster ID': true_cell_type,
            'Predicted General Cell Type': main_cell_type,
            'Predicted Detailed Cell Type': sub_cell_types,
            'Possible Mixed Cell Types': possible_mixed_cell_types,
            'Marker Number': marker_number,
            'Marker List': marker_list_str,
            'Iterations': iterations,
            'Model': model,
            'Provider': provider,
            'Tissue': tissue,
            'Species': species,
            'Additional Info': additional_info or "N/A",
            'Reference Used': ref_used,
            'Complexity Score': ref_complexity,
            'Conversation History': raw_conversation_history
        })

        full_data.append([
            true_cell_type,
            main_cell_type,
            sub_cell_types,
            possible_mixed_cell_types,
            marker_number,
            marker_list_str,
            iterations,
            model,
            provider,
            tissue,
            species,
            additional_info or "N/A",
            ref_used,
            ref_complexity,
            conversation_history
        ])
        summary_data.append([
            true_cell_type,
            main_cell_type,
            sub_cell_types,
            possible_mixed_cell_types,
            marker_list_str,
            iterations,
            model,
            provider,
            tissue,
            species,
            ref_used
        ])

    # Generate output filenames
    base_name = os.path.splitext(output_name)[0]
    full_csv_name = f"{base_name}_full.csv"
    summary_csv_name = f"{base_name}_summary.csv"

    # Make sure the output directory exists
    output_dir = os.path.dirname(full_csv_name)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Sort the data by Cluster ID with natural ordering
    full_data.sort(key=lambda x: natural_sort_key(x[0]))
    summary_data.sort(key=lambda x: natural_sort_key(x[0]))

    # Write CSV files with reference columns
    write_csv(full_csv_name,
              ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
               'Possible Mixed Cell Types', 'Marker Number', 'Marker List', 'Iterations',
               'Model', 'Provider', 'Tissue', 'Species', 'Additional Info',
               'Reference Used', 'Complexity Score', 'Conversation History'],
              full_data)

    write_csv(summary_csv_name,
              ['Cluster ID', 'Predicted General Cell Type', 'Predicted Detailed Cell Type',
               'Possible Mixed Cell Types', 'Marker List', 'Iterations', 'Model', 'Provider',
               'Tissue', 'Species', 'Reference Used'],
              summary_data)

    # Generate HTML report
    html_report_name = f"{base_name}_report.html"
    try:
        try:
            from .generate_batch_report import generate_batch_html_report_from_data
        except ImportError:
            from generate_batch_report import generate_batch_html_report_from_data
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

    return results


# =============================================================================
# RE-EXPORTS FOR BACKWARD COMPATIBILITY
# =============================================================================
# The functions below have been moved to separate modules but are re-exported
# here to maintain backward compatibility with existing code that imports from
# tools_function.py (including R package via reticulate).

# From annotation_boost.py
try:
    from .annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task
except ImportError:
    from annotation_boost import runCASSIA_annotationboost, runCASSIA_annotationboost_additional_task

# From generate_reports.py (HTML report generation)
try:
    from .generate_reports import (
        generate_analysis_html_report as generate_html_report,  # Renamed for backward compat
        process_single_report,
        generate_score_index_page as generate_index_page,  # Renamed for backward compat
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
    from .pipeline import runCASSIA_pipeline
except ImportError:
    from pipeline import runCASSIA_pipeline

# From marker_utils.py (marker loading utilities)
try:
    from .marker_utils import loadmarker, list_available_markers
except ImportError:
    from marker_utils import loadmarker, list_available_markers

# From scoring.py (scoring functions)
try:
    from .scoring import (
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
    from .Uncertainty_quantification import (
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
