import os
import json
import requests
import pandas as pd
import csv
import concurrent.futures
from functools import partial
from typing import Dict, Any, Optional, Union, List
try:
    from CASSIA.core.llm_utils import call_llm
    from CASSIA.core.progress_tracker import BatchProgressTracker
except ImportError:
    try:
        from .llm_utils import call_llm
        from .progress_tracker import BatchProgressTracker
    except ImportError:
        from llm_utils import call_llm
        from progress_tracker import BatchProgressTracker

def merge_annotations(
    csv_path: str,
    output_path: Optional[str] = None,
    provider: str = "openai",
    model: Optional[str] = None,
    api_key: Optional[str] = None,
    additional_context: Optional[str] = None,
    batch_size: int = 20,
    detail_level: str = "broad"  # New parameter for controlling grouping granularity
) -> pd.DataFrame:
    """
    Agent function that reads a CSV file with cell cluster annotations and merges/groups them.
    
    Args:
        csv_path: Path to the CSV file containing cluster annotations
        output_path: Path to save the results (if None, returns DataFrame without saving)
        provider: LLM provider to use ("openai", "anthropic", or "openrouter")
        model: Specific model to use (if None, uses default for provider)
        api_key: API key for the provider (if None, gets from environment)
        additional_context: Optional domain-specific context to help with annotation
        batch_size: Number of clusters to process in each LLM call (for efficiency)
        detail_level: Level of detail for the groupings:
                     - "broad": More general cell categories (e.g., "Myeloid cells" for macrophages and dendritic cells)
                     - "detailed": More specific groupings that still consolidate very specific clusters
                     - "very_detailed": Most specific groupings with normalized and consistent naming
        
    Returns:
        DataFrame with original annotations and suggested cell groupings
    """
    # Validate detail_level parameter
    if detail_level not in ["broad", "detailed", "very_detailed"]:
        raise ValueError("detail_level must be one of: 'broad', 'detailed', or 'very_detailed'")
    
    # Set column name for results based on detail level
    result_column_map = {
        "broad": "Merged_Grouping_1",
        "detailed": "Merged_Grouping_2",
        "very_detailed": "Merged_Grouping_3"
    }
    result_column = result_column_map[detail_level]
    
    # Read the CSV file
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        raise ValueError(f"Error reading CSV file: {str(e)}")
    
    # Map expected column names to actual column names
    # Check for the expected column names based on the new information
    column_mapping = {
        "cluster": "Cluster ID",  # The cluster identifier column
        "general_annotation": "Predicted General Cell Type",
        "subtype_annotation": "Predicted Detailed Cell Type"
    }
    
    # Verify that we found the necessary columns
    missing_columns = [col for col in column_mapping.values() if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Required columns not found: {', '.join(missing_columns)}. Available columns: {', '.join(df.columns)}")
    
    # Create a working copy of the DataFrame
    working_df = df.copy()
    
    # Extract first subtype if the subtype column contains comma-separated values
    subtype_col = column_mapping["subtype_annotation"]
    if working_df[subtype_col].dtype == 'object':  # Check if the column contains strings
        # Extract the first subtype from comma-separated list
        working_df["processed_subtype"] = working_df[subtype_col].apply(
            lambda x: x.split(',')[0].strip() if isinstance(x, str) and ',' in x else x
        )
    else:
        working_df["processed_subtype"] = working_df[subtype_col]
    
    # Create a new column for the merged annotations
    working_df[result_column] = ""
    
    # Process in batches for efficiency
    total_rows = len(working_df)
    for i in range(0, total_rows, batch_size):
        batch_end = min(i + batch_size, total_rows)
        batch = working_df.iloc[i:batch_end]
        
        # Prepare prompt for LLM based on detail level
        prompt = _create_annotation_prompt(
            batch, 
            cluster_col=column_mapping["cluster"],
            general_col=column_mapping["general_annotation"],
            subtype_col="processed_subtype",
            additional_context=additional_context,
            detail_level=detail_level
        )
        
        # Call LLM to get suggested groupings
        try:
            response = call_llm(
                prompt=prompt,
                provider=provider,
                model=model,
                api_key=api_key,
                temperature=0.3,  # Lower temperature for more consistent results
                max_tokens=4096,
                system_prompt="You are an expert cell biologist specializing in single-cell analysis. Your task is to analyze cluster annotations and suggest general cell groupings."
            )
            
            # Parse LLM response and update DataFrame
            groupings = _parse_llm_response(response, batch.index)
            for idx, grouping in groupings.items():
                working_df.at[idx, result_column] = grouping
        except Exception as e:
            # Silently continue - errors will be visible from incomplete results
            pass
    
    # Add the result column to the original DataFrame
    df[result_column] = working_df[result_column]
    
    # Save results if output path is provided
    if output_path:
        df.to_csv(output_path, index=False)
    
    return df

def _create_annotation_prompt(
    batch: pd.DataFrame, 
    cluster_col: str,
    general_col: str, 
    subtype_col: str,
    additional_context: Optional[str],
    detail_level: str = "broad"
) -> str:
    """
    Create a prompt for the LLM to suggest groupings based on cluster annotations.
    
    Args:
        batch: DataFrame batch containing clusters to process
        cluster_col: Name of the cluster ID column
        general_col: Name of the general annotation column
        subtype_col: Name of the subtype annotation column
        additional_context: Optional domain-specific context
        detail_level: Level of detail for the groupings ("broad", "detailed", or "very_detailed")
        
    Returns:
        Formatted prompt string
    """
    # Format clusters data
    clusters_text = "\n".join([
        f"Cluster {row[cluster_col]}: General annotation: {row[general_col]}, Subtype: {row[subtype_col]}"
        for _, row in batch.iterrows()
    ])
    
    # Create the prompt based on detail level
    if detail_level == "broad":
        # Broad groupings prompt (original)
        prompt = f"""I have single-cell RNA-seq cluster annotations and need to suggest broader cell groupings.
For each cluster, I'll provide the general annotation and subtype annotation.
Based on these annotations, suggest an appropriate broader cell grouping category.

For example:
- "macrophage, inflammatory macrophage" → "Myeloid cells"
- "CD4 T cell, naive CD4 T cell" → "T cells"
- "B cell, memory B cell" → "B cells"

Use general cell lineage categories when possible, combining related cell types into a single group.
Prioritize creating broader categories that span multiple specific cell types.

Annotations to process:
{clusters_text}

Please respond with a JSON object where keys are cluster identifiers and values are the suggested groupings. 
For example:
{{
  "1": "Myeloid cells",
  "2": "T cells"
}}
"""
    elif detail_level == "detailed":
        # Detailed groupings prompt
        prompt = f"""I have single-cell RNA-seq cluster annotations and need to suggest intermediate-level cell groupings.
For each cluster, I'll provide the general annotation and subtype annotation.
Based on these annotations, suggest a moderately specific cell grouping that balances detail and generality.

For example:
- "macrophage, inflammatory macrophage" → "Macrophages" (not as broad as "Myeloid cells")
- "CD4 T cell, naive CD4 T cell" → "CD4 T cells" (more specific than just "T cells")
- "CD8 T cell, cytotoxic CD8 T cell" → "CD8 T cells" (more specific than just "T cells")
- "B cell, memory B cell" → "B cells" (specific cell type)

Maintain biological specificity when important, but still group very similar subtypes together.
Aim for a middle ground - not too general, but also not too specific.
The grouping should be more detailed than broad categories like "Myeloid cells" or "Lymphoid cells", 
but should still consolidate highly specific annotations.

Annotations to process:
{clusters_text}

Please respond with a JSON object where keys are cluster identifiers and values are the suggested groupings. 
For example:
{{
  "1": "Macrophages",
  "2": "CD4 T cells",
  "3": "CD8 T cells"
}}
"""
    else:  # very_detailed
        # Very detailed groupings prompt
        prompt = f"""I have single-cell RNA-seq cluster annotations and need to normalize and standardize cell type names 
while preserving the most specific and detailed biological information.
For each cluster, I'll provide the general annotation and subtype annotation.

Your task is to create a consistent and standardized cell type label that:
1. Maintains the highest level of biological specificity from the annotations
2. Uses consistent nomenclature across similar cell types
3. Follows standard cell type naming conventions
4. Preserves functional or activation state information when present
5. Normalizes naming variations (e.g., "inflammatory macrophage" vs "M1 macrophage" should use one consistent term)

Examples:
- "macrophage, inflammatory macrophage" → "Inflammatory macrophages" (preserve activation state)
- "CD4 T cell, naive CD4 T cell" → "Naive CD4+ T cells" (preserve naive state, standardize CD4+)
- "CD8 T cell, cytotoxic CD8 T cell" → "Cytotoxic CD8+ T cells" (preserve function, standardize CD8+)
- "dendritic cell, plasmacytoid dendritic cell" → "Plasmacytoid dendritic cells" (preserve specific subtype)
- "B cell, memory B cell" → "Memory B cells" (preserve memory state)
- "NK cell, CD56bright NK cell" → "CD56bright NK cells" (preserve specific marker)

Annotations to process:
{clusters_text}

Please respond with a JSON object where keys are cluster identifiers and values are the normalized, specific cell type labels.
For example:
{{
  "1": "Inflammatory macrophages",
  "2": "Naive CD4+ T cells",
  "3": "Memory B cells"
}}
"""
    
    # Add additional context if provided
    if additional_context:
        prompt += f"\n\nAdditional context that may help with the analysis:\n{additional_context}"
    
    return prompt

def _parse_llm_response(response: str, indices: pd.Index) -> Dict[int, str]:
    """
    Parse the LLM response to extract suggested groupings.
    
    Args:
        response: LLM response text
        indices: DataFrame indices for the batch
        
    Returns:
        Dictionary mapping DataFrame indices to suggested groupings
    """
    groupings = {}
    
    # Try to find and parse JSON in the response
    try:
        # Extract JSON if it's embedded in text
        import re
        json_match = re.search(r'({[\s\S]*})', response)
        if json_match:
            json_str = json_match.group(1)
            parsed = json.loads(json_str)
            
            # Map the parsed results to DataFrame indices
            for i, (cluster_id, grouping) in enumerate(parsed.items()):
                if i < len(indices):
                    groupings[indices[i]] = grouping
        else:
            # Fallback: Try parsing line by line
            lines = [line.strip() for line in response.split('\n') if line.strip()]
            for i, line in enumerate(lines):
                if i < len(indices):
                    if ':' in line:
                        groupings[indices[i]] = line.split(':', 1)[1].strip()
    except Exception as e:
        # Fallback: Use the raw response on parse error
        for i, idx in enumerate(indices):
            groupings[idx] = "Error parsing response"
    
    return groupings

def merge_annotations_all(
    csv_path: str,
    output_path: Optional[str] = None,
    provider: str = "openai",
    model: Optional[str] = None,
    api_key: Optional[str] = None,
    additional_context: Optional[str] = None,
    batch_size: int = 20
) -> pd.DataFrame:
    """
    Process all three detail levels in parallel and return a combined DataFrame.
    
    Args:
        csv_path: Path to the CSV file containing cluster annotations
        output_path: Path to save the results (if None, returns DataFrame without saving)
        provider: LLM provider to use ("openai", "anthropic", or "openrouter")
        model: Specific model to use (if None, uses default for provider)
        api_key: API key for the provider (if None, gets from environment)
        additional_context: Optional domain-specific context to help with annotation
        batch_size: Number of clusters to process in each LLM call (for efficiency)
        
    Returns:
        DataFrame with original annotations and all three levels of suggested cell groupings
    """
    # Define the detail levels to process
    detail_levels = ["broad", "detailed", "very_detailed"]

    # Initialize progress tracker
    tracker = BatchProgressTracker(total=len(detail_levels), title="Annotation Merging")

    # Create a partial function with common arguments
    merge_func = partial(
        merge_annotations,
        csv_path=csv_path,
        output_path=None,  # We'll save the combined result at the end
        provider=provider,
        model=model,
        api_key=api_key,
        additional_context=additional_context,
        batch_size=batch_size
    )

    # Use ThreadPoolExecutor to run in parallel with 3 threads (instead of ProcessPoolExecutor)
    results = {}
    errors = []
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            # Submit all tasks and mark them as started
            future_to_level = {}
            for level in detail_levels:
                future = executor.submit(merge_func, detail_level=level)
                future_to_level[future] = level
                tracker.start_task(level)

            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_level):
                level = future_to_level[future]
                try:
                    df = future.result()
                    results[level] = df
                    tracker.complete_task(level)
                except Exception as e:
                    tracker.complete_task(level)
                    errors.append(f"Error processing {level}: {str(e)}")
    finally:
        tracker.finish()
    
    # If we have results, combine them
    if results:
        # Start with the first detail level's DataFrame
        combined_df = None
        for level in detail_levels:
            if level in results:
                if combined_df is None:
                    combined_df = results[level].copy()
                else:
                    # Get the column name for this level
                    result_column_map = {
                        "broad": "Merged_Grouping_1",
                        "detailed": "Merged_Grouping_2",
                        "very_detailed": "Merged_Grouping_3"
                    }
                    result_column = result_column_map[level]
                    
                    # Add this level's results to the combined DataFrame
                    combined_df[result_column] = results[level][result_column]
        
        # Save combined results if output path is provided
        if output_path and combined_df is not None:
            combined_df.to_csv(output_path, index=False)

        # Print any errors that occurred
        if errors:
            for error in errors:
                print(error)

        return combined_df
    else:
        raise ValueError("All parallel processing tasks failed. Check logs for details.")








