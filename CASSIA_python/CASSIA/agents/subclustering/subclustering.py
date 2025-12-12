try:
    from CASSIA.engine.tools_function import *
    from CASSIA.core.model_settings import get_agent_default
except ImportError:
    try:
        from .tools_function import *
        from ..core.model_settings import get_agent_default
    except ImportError:
        from tools_function import *
        from model_settings import get_agent_default

try:
    from CASSIA.core.llm_utils import *
except ImportError:
    try:
        from .llm_utils import *
    except ImportError:
        from llm_utils import *

import pandas as pd


def _get_get_top_markers():
    """Lazy import of get_top_markers to avoid circular imports."""
    try:
        from CASSIA.core.marker_utils import get_top_markers
        return get_top_markers
    except ImportError:
        try:
            from ...core.marker_utils import get_top_markers
            return get_top_markers
        except ImportError:
            from marker_utils import get_top_markers
            return get_top_markers


def subcluster_agent_annotate_subcluster(user_message, model=None, temperature=None, provider="openrouter"):
    """
    Unified function to call LLM for subcluster annotation.

    Args:
        user_message: The prompt message for subcluster annotation
        model: Model to use (defaults to provider's default if None)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)

    Returns:
        The generated annotation as a string
    """
    # Set default model and temperature based on provider if not specified
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]

    # Use the unified call_llm function
    result = call_llm(
        prompt=user_message,
        provider=provider,
        model=model,
        temperature=temperature,
        max_tokens=4096
    )

    return result if result else ''



def construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=50):
    # Process DataFrame if it has more than 2 columns
    if len(marker.columns) > 2:
        print(f"Processing input dataframe to get top {n_genes} markers")
        get_top_markers = _get_get_top_markers()
        marker = get_top_markers(marker, n_genes=n_genes)
    else:
        print("Using input dataframe directly as it appears to be pre-processed (2 columns)")
        marker = marker.copy()
    
    # Initialize the prompt with the major cluster information
    prompt = f"""

You are an expert biologist specializing in cell type annotation, with deep expertise in immunology, cancer biology, and developmental biology.You will be given sets of highly expressed markers ranked by significance for some subclusters from the {major_cluster_info} cluster, identify what is the most likely top2 cell type each marker set implies.

Take a deep breath and work step by step. You'd better do a really good job or 1000 grandma are going to be in danger.
You will be tipped $10,000 if you do a good job.

For each output, provide:
1.Key marker:
2.Explanation:
3.Most likely top2 cell types:

Remember these subclusters are from a {major_cluster_info} big cluster. You must include all clusters mentioned in the analysis.
"""

    # Iterate over each row in the DataFrame
    for i, (index, row) in enumerate(marker.iterrows(), start=1):
        cluster_name = row.iloc[0]  # Use iloc for positional indexing
        markers = row.iloc[1]       # Use iloc for positional indexing
        prompt += f"{i}.{markers}\n"

    return prompt



def annotate_subclusters(marker, major_cluster_info, model=None, temperature=None, provider="openrouter", n_genes=50):
    """
    Annotate subclusters using an LLM.

    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        model: Model to use (defaults to provider's subclustering default)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use

    Returns:
        The generated annotation as a string
    """
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt, model=model, temperature=temperature, provider=provider)
    return output_text



def extract_subcluster_results_with_llm_multiple_output(analysis_text, provider="openrouter", model=None, temperature=None):
    """
    Extract multiple output results from subcluster analysis text using XML tags.
    Uses brief reasoning to save tokens in batch runs.

    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use
        temperature: Temperature for generation (0-1)

    Returns:
        Extracted results as string with XML-tagged format
    """
    # Define the prompt to instruct the LLM with XML output format (brief reasoning for batch)
    prompt = f"""You are an expert in analyzing celltype annotation for subclusters.

Extract the cell type annotations from the following analysis. For each cluster, output in this exact XML format:

<cluster id="1">
<celltype1>first cell type</celltype1>
<celltype2>second cell type</celltype2>
<reason>brief reason</reason>
</cluster>

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    result = subcluster_agent_annotate_subcluster(prompt, provider=provider, model=model, temperature=temperature)
    return result




def extract_subcluster_results_with_llm(analysis_text, provider="openrouter", model=None, temperature=None):
    """
    Extract results with reasons from subcluster analysis text using XML tags.

    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use
        temperature: Temperature for generation (0-1)

    Returns:
        Extracted results as string with XML-tagged format
    """
    # Define the prompt to instruct the LLM with XML output format
    prompt = f"""You are an expert in analyzing celltype annotation for subclusters.

Extract the cell type annotations from the following analysis. For each cluster, output in this exact XML format:

<cluster id="1">
<celltype1>first cell type</celltype1>
<celltype2>second cell type</celltype2>
<reason>the complete explanation from the analysis</reason>
</cluster>

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    result = subcluster_agent_annotate_subcluster(prompt, provider=provider, model=model, temperature=temperature)
    return result



def write_results_to_csv(results, output_name='subcluster_results'):
    """
    Extract cell type results from LLM output (XML format) and write to CSV file

    Args:
        results (str): LLM analysis results with XML-tagged clusters
        output_name (str): Base name for output file (will add .csv if not present)

    Returns:
        pandas.DataFrame: DataFrame containing the extracted results
    """
    # Add .csv suffix if not present
    if not output_name.lower().endswith('.csv'):
        output_name = output_name + '.csv'

    results_str = str(results)

    # Parse XML-tagged clusters: <cluster id="1">...</cluster>
    cluster_pattern = r'<cluster[^>]*id=["\']?(\d+)["\']?[^>]*>(.*?)</cluster>'
    cluster_matches = re.findall(cluster_pattern, results_str, re.DOTALL | re.IGNORECASE)

    rows = []
    if cluster_matches:
        for cluster_id, content in cluster_matches:
            # Extract celltype1
            ct1_match = re.search(r'<celltype1>(.*?)</celltype1>', content, re.DOTALL | re.IGNORECASE)
            celltype1 = ct1_match.group(1).strip() if ct1_match else 'Unknown'

            # Extract celltype2
            ct2_match = re.search(r'<celltype2>(.*?)</celltype2>', content, re.DOTALL | re.IGNORECASE)
            celltype2 = ct2_match.group(1).strip() if ct2_match else 'Unknown'

            # Extract reason
            reason_match = re.search(r'<reason>(.*?)</reason>', content, re.DOTALL | re.IGNORECASE)
            reason = reason_match.group(1).strip() if reason_match else ''

            rows.append([cluster_id, celltype1, celltype2, '', reason])

    if rows:
        df = pd.DataFrame(rows, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'])
        df.to_csv(output_name, index=False)
        print(f"Results have been written to {output_name}")
        return df

    # If XML parsing failed, save raw results for debugging
    raw_output_file = f"{output_name}.txt"
    with open(raw_output_file, "w") as f:
        f.write(results_str)

    raise RuntimeError(
        f"\n{'='*60}\n"
        f"SUBCLUSTERING FAILED - Could not parse LLM results\n"
        f"{'='*60}\n"
        f"Raw output saved to: {raw_output_file}\n"
        f"Expected XML format: <cluster id=\"1\"><celltype1>...</celltype1>...</cluster>\n"
        f"{'='*60}"
    )



def runCASSIA_subclusters(marker, major_cluster_info, output_name,
                       model=None, temperature=None, provider="openrouter", n_genes=50):
    """
    Process subclusters from marker data and generate annotated results.

    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        output_name: Base name for output file (will add .csv if not present)
        model: Model name to use (defaults to provider's subclustering default)
        temperature: Temperature parameter for API calls (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use for analysis

    Returns:
        None: Results are saved to a CSV file
    """
    # Apply agent defaults if model or temperature not specified
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]

    # Construct prompt and get analysis from LLM
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt, model=model, temperature=temperature, provider=provider)

    # Extract structured results from the analysis text
    results = extract_subcluster_results_with_llm(output_text, provider=provider, model=model, temperature=temperature)
    # print(results)  # Remove or comment out this line to avoid showing conversation history
    
    # Save results to CSV
    write_results_to_csv(results, output_name)

    # --- Generate HTML report for the single run CSV ---
    try:
        from CASSIA.reports.generate_reports import process_evaluation_csv
    except ImportError:
        try:
            from ...reports.generate_reports import process_evaluation_csv
        except ImportError:
            from reports.generate_reports import process_evaluation_csv
    import os
    csv_file = output_name if output_name.lower().endswith('.csv') else output_name + '.csv'
    if os.path.exists(csv_file):
        process_evaluation_csv(csv_file, overwrite=True, model_name=model)
    
    return None



def runCASSIA_n_subcluster(n, marker, major_cluster_info, base_output_name,
                          model=None, temperature=None,
                          provider="openrouter", max_workers=5, n_genes=50):
    """
    Run multiple subcluster analyses in parallel and save results.

    Args:
        n: Number of analyses to run
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        base_output_name: Base name for output files
        model: Model name to use (defaults to provider's subclustering_n default)
        temperature: Temperature parameter for API calls (defaults to subclustering_n default: 0.3)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        max_workers: Maximum number of parallel workers
        n_genes: Number of top genes to use for analysis

    Returns:
        None: Results are saved to CSV files
    """
    # Apply agent defaults for n-times subclustering (uses subclustering_n for variability)
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering_n", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]
    
    def run_single_analysis(i):
        # Run the annotation process
        output_text = annotate_subclusters(marker, major_cluster_info,
                                         model=model, temperature=temperature, provider=provider, n_genes=n_genes)

        # Extract results using XML format
        results = extract_subcluster_results_with_llm_multiple_output(output_text, provider=provider, model=model, temperature=temperature)
        results_str = str(results)

        # Parse XML-tagged clusters: <cluster id="1">...</cluster>
        cluster_pattern = r'<cluster[^>]*id=["\']?(\d+)["\']?[^>]*>(.*?)</cluster>'
        cluster_matches = re.findall(cluster_pattern, results_str, re.DOTALL | re.IGNORECASE)

        rows = []
        if cluster_matches:
            for cluster_id, content in cluster_matches:
                # Extract celltype1
                ct1_match = re.search(r'<celltype1>(.*?)</celltype1>', content, re.DOTALL | re.IGNORECASE)
                celltype1 = ct1_match.group(1).strip() if ct1_match else 'Unknown'

                # Extract celltype2
                ct2_match = re.search(r'<celltype2>(.*?)</celltype2>', content, re.DOTALL | re.IGNORECASE)
                celltype2 = ct2_match.group(1).strip() if ct2_match else 'Unknown'

                # Extract reason
                reason_match = re.search(r'<reason>(.*?)</reason>', content, re.DOTALL | re.IGNORECASE)
                reason = reason_match.group(1).strip() if reason_match else ''

                rows.append([cluster_id, celltype1, celltype2, '', reason])

        if rows:
            df = pd.DataFrame(rows, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'])
        else:
            raise RuntimeError(
                f"Iteration {i+1}: Could not parse LLM response. "
                f"Expected XML format: <cluster id=\"1\"><celltype1>...</celltype1>...</cluster>. "
                f"Response preview: {results_str[:200]}"
            )

        try:
            # Try to get top markers, but handle the case where required columns are missing
            try:
                get_top_markers = _get_get_top_markers()
                marker_df = get_top_markers(marker, n_genes=n_genes)
            except KeyError as e:
                # If get_top_markers fails due to missing columns, use the original marker dataframe
                print(f"Warning: {str(e)}. Using original marker dataframe.")
                marker_df = marker.copy()
            
            # Convert types to ensure compatibility
            df['Result ID'] = df['Result ID'].astype(str)

            # Make a copy of the original values before swapping
            original_cluster_ids = df['Result ID'].copy()

            # Check if marker_df has at least one row and one column
            if marker_df.shape[0] > 0 and marker_df.shape[1] > 0:
                original_marker_first_col = marker_df.iloc[:, 0].copy()

                # Perform the swap safely - only if there are enough rows in both dataframes
                min_rows = min(len(df), len(marker_df))
                if min_rows > 0:
                    df.loc[:min_rows-1, 'Result ID'] = original_marker_first_col[:min_rows].values
                    # Only update marker_df if we're actually going to use it later
                    # marker_df.iloc[:min_rows, 0] = original_true_cell_types[:min_rows].values
            else:
                print("Warning: Marker dataframe is empty or has no columns. Skipping column swap.")
        except Exception as e:
            print(f"Warning: Error during column swap: {str(e)}. Continuing without swapping.")

        # Write the DataFrame to a CSV file with an index
        indexed_csv_file_path = f'{base_output_name}_{i+1}.csv'
        df.to_csv(indexed_csv_file_path, index=False)
        
        return indexed_csv_file_path

    failed_iterations = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_analysis, i): i for i in range(n)}
        result_files = []
        for future in as_completed(futures):
            i = futures[future]
            try:
                result_file = future.result()
                print(f"Results for iteration {i+1} have been written to {result_file}")
                result_files.append(result_file)
            except Exception as exc:
                print(f"Iteration {i+1} generated an exception: {exc}")
                failed_iterations.append((i+1, str(exc)))

    # Check if ALL iterations failed - this is a total failure
    if len(result_files) == 0 and len(failed_iterations) > 0:
        error_sample = failed_iterations[0][1][:200] if failed_iterations else "Unknown"
        raise RuntimeError(
            f"\n{'='*60}\n"
            f"SUBCLUSTERING BATCH FAILED - All {n} iterations failed\n"
            f"{'='*60}\n"
            f"Sample error: {error_sample}\n"
            f"{'='*60}"
        )

    # Warn about partial failures
    if failed_iterations:
        print(f"\nWarning: {len(failed_iterations)} of {n} iterations failed:")
        for iter_num, err in failed_iterations[:5]:
            print(f"  - Iteration {iter_num}: {err[:80]}...")
        if len(failed_iterations) > 5:
            print(f"  ... and {len(failed_iterations) - 5} more")

    # --- Generate HTML reports for all batch CSVs ---
    try:
        from CASSIA.reports.generate_reports import process_evaluation_csv, create_index_html
    except ImportError:
        try:
            from ...reports.generate_reports import process_evaluation_csv, create_index_html
        except ImportError:
            from reports.generate_reports import process_evaluation_csv, create_index_html
    import os

    # Generate HTML report for each CSV
    for csv_file in result_files:
        if os.path.exists(csv_file):
            process_evaluation_csv(csv_file, overwrite=True)

    # Create an index.html summary in the same directory as the first result file
    if result_files:
        output_dir = os.path.dirname(result_files[0]) or '.'
        create_index_html(result_files, output_dir)
        print(f"Batch HTML reports and index generated in {output_dir}")


def test_custom_api_parsing():
    """
    Test function to simulate a response from a custom API provider and test the parsing functionality.
    This is useful for debugging the parsing logic without making actual API calls.
    """
    import pandas as pd
    import os
    try:
        from CASSIA.engine.tools_function import get_top_markers
    except ImportError:
        try:
            from .tools_function import get_top_markers
        except ImportError:
            from tools_function import get_top_markers

    # Sample structured response that mimics what DeepSeek or other custom APIs might return
    sample_response = [
        {
            'cluster': 1,
            'key_markers': 'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'explanation': 'The presence of IL7R, CD8A, and CD8B suggests a CD8+ T cell identity. CCL4 is associated with effector functions, while KLRB1 (CD161) and ITK indicate a memory-like or tissue-resident phenotype.',
            'most_likely_top2_cell_types': ['CD8+ memory T cells', 'Tissue-resident memory CD8+ T cells (TRM)']
        },
        {
            'cluster': 2,
            'key_markers': 'LAYN, HAVCR2 (TIM-3), TIGIT, IKZF2, KLRC2, KLRC3',
            'explanation': 'LAYN (Lag-3) and HAVCR2 (TIM-3) are markers of exhausted or chronically stimulated CD8+ T cells.',
            'most_likely_top2_cell_types': ['Exhausted CD8+ T cells', 'NK-like CD8+ T cells']
        },
        {
            'cluster': 3,
            'key_markers': 'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'explanation': 'GZMK, GZMH, PRF1, and NKG7 are markers of cytotoxic activity, typical of effector CD8+ T cells.',
            'most_likely_top2_cell_types': ['Effector CD8+ T cells', 'Central memory CD8+ T cells']
        },
        {
            'cluster': 4,
            'key_markers': 'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'explanation': 'WFDC2, CEACAM7, and CLDN8 are markers associated with epithelial or secretory cells.',
            'most_likely_top2_cell_types': ['Regulatory CD8+ T cells', 'Epithelial-like CD8+ T cells (rare subset)']
        }
    ]
    
    # Test the write_results_to_csv function
    print("Testing write_results_to_csv with structured data...")
    df = write_results_to_csv(sample_response, output_name='test_parsing_result')
    print(f"Generated DataFrame:\n{df}")
    
    # Test the extract_subcluster_results_with_llm function
    print("\nTesting extract_subcluster_results_with_llm with structured data...")
    result = extract_subcluster_results_with_llm(sample_response)
    print(f"Result type: {type(result)}")
    
    # Test the extract_subcluster_results_with_llm_multiple_output function
    print("\nTesting extract_subcluster_results_with_llm_multiple_output with structured data...")
    result = extract_subcluster_results_with_llm_multiple_output(sample_response)
    print(f"Result type: {type(result)}")
    
    # Create a simple marker dataframe for testing
    print("\nCreating test marker dataframe...")
    test_marker_df = pd.DataFrame({
        'cluster': ['cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5', 'cluster6'],
        'markers': [
            'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3',
            'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'GNLY, KLRF1, FCER1G, TYROBP, CD38, KIR2DL4',
            'LPL, SNAI2, HAND2, SOX2, NES, PDGFRA'
        ]
    })
    print(f"Created test marker dataframe with shape: {test_marker_df.shape}")
    
    # Save the test dataframe to a temporary file
    temp_file = 'test_markers_temp.csv'
    test_marker_df.to_csv(temp_file, index=False)
    print(f"Saved test marker dataframe to {temp_file}")
    
    try:
        # Create a mock function to simulate annotate_subclusters
        def mock_annotate(*args, **kwargs):
            return sample_response
            
        # Save the original function
        original_annotate = globals()['annotate_subclusters']
        
        # Replace with mock function
        globals()['annotate_subclusters'] = mock_annotate
        
        print("\nTesting runCASSIA_n_subcluster with mock data...")
        # Run with n=1 to test a single iteration
        runCASSIA_n_subcluster(
            n=1,
            marker=test_marker_df,
            major_cluster_info="cd8 t cell",
            base_output_name="test_n_subcluster",
            model="dummy-model",
            provider="dummy-provider"
        )
        
        # Check if the output file was created
        output_file = "test_n_subcluster_1.csv"
        if os.path.exists(output_file):
            print(f"Successfully created output file: {output_file}")
            result_df = pd.read_csv(output_file)
            print(f"Output file contents:\n{result_df}")
        else:
            print(f"Error: Output file {output_file} was not created")
        
        # Restore the original function
        globals()['annotate_subclusters'] = original_annotate
        
        print("\nAll tests completed.")
        
    except Exception as e:
        print(f"Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()
    finally:
        # Clean up temporary files
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"Removed temporary file: {temp_file}")

# Uncomment to run the test function when this file is executed directly
# if __name__ == "__main__":
#     test_custom_api_parsing()

