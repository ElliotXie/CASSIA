try:
    from CASSIA.engine.tools_function import *
except ImportError:
    try:
        from .tools_function import *
    except ImportError:
        from tools_function import *

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


def subcluster_agent_annotate_subcluster(user_message, model=None, temperature=0, provider="anthropic"):
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
    # Set default model based on provider if not specified
    if model is None:
        if provider == "openai":
            model = "gpt-4o"
        elif provider == "anthropic":
            model = "claude-3-5-sonnet-20241022"
        elif provider == "openrouter":
            model = "anthropic/claude-3.5-sonnet"
        elif provider.startswith("http"):
            # For custom API endpoints, use a default model if none specified
            model = model or "deepseek-chat"
            print(f"Using model: {model} with custom provider: {provider}")
    
    # Add JSON tags for providers that need them
    modified_message = user_message
    if provider == "anthropic" or provider == "openrouter" or provider.startswith("http"):
        # Add JSON tags to help with structured output for these providers
        if not "<json>" in modified_message.lower():
            modified_message = f"{modified_message}\n\nPlease format your response as JSON:\n<json>\n{{\"response\": \"Your detailed analysis here\"}}\n</json>"
    
    # Use the unified call_llm function
    result = call_llm(
        prompt=modified_message,
        provider=provider,
        model=model,
        temperature=temperature,
        max_tokens=4096
    )
    
    # Process the result for providers that use JSON tags
    if provider == "anthropic" or provider == "openrouter" or provider.startswith("http"):
        # Extract content from JSON tags if present
        import re
        import json
        
        # First check for JSON format with ```json or <json> tags
        json_match = re.search(r'(?:```json|<json>)(.*?)(?:```|</json>)', result, re.DOTALL)
        if json_match:
            try:
                json_content = json_match.group(1).strip()
                # Try to parse as JSON
                parsed = json.loads(json_content)
                if isinstance(parsed, dict) and "response" in parsed:
                    # If it's a dict with a response key, return the response
                    if isinstance(parsed["response"], list):
                        return parsed["response"]
                    return parsed["response"]
                return parsed
            except json.JSONDecodeError:
                pass  # If JSON parsing fails, continue with other methods
        
        # Check if the entire response is a JSON string (common with DeepSeek API)
        if result.strip().startswith('{') and result.strip().endswith('}'):
            try:
                parsed = json.loads(result)
                if isinstance(parsed, dict) and "response" in parsed:
                    # If it's a dict with a response key, return the response
                    if isinstance(parsed["response"], list):
                        return parsed["response"]
                    return parsed["response"]
                return parsed
            except json.JSONDecodeError:
                pass  # If JSON parsing fails, return the original result
    
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



def annotate_subclusters(marker, major_cluster_info, model="claude-3-5-sonnet-20241022", temperature=0, provider="anthropic", n_genes=50):
    """
    Annotate subclusters using an LLM.
    
    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        model: Model to use (defaults to Claude 3.5 Sonnet)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use
        
    Returns:
        The generated annotation as a string
    """
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt, model=model, temperature=temperature, provider=provider)
    return output_text



def extract_subcluster_results_with_llm_multiple_output(analysis_text, provider="anthropic", model="claude-3-5-sonnet-20241022", temperature=0):
    """
    Extract multiple output results from subcluster analysis text.
    
    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use (defaults to Claude 3.5 Sonnet)
        temperature: Temperature for generation (0-1)
        
    Returns:
        Extracted results in the format: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.
        For custom APIs, may return a list of dictionaries directly.
    """
    # If the analysis_text is already in a structured format (list of dictionaries),
    # we can return it directly for processing
    if isinstance(analysis_text, list) and all(isinstance(item, dict) for item in analysis_text):
        print(f"Analysis text is already in structured format (list of {len(analysis_text)} dictionaries)")
        return analysis_text
    
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    result = subcluster_agent_annotate_subcluster(prompt, provider=provider, model=model, temperature=temperature)

    # Check if the result is already in structured format (e.g., from DeepSeek API)
    if isinstance(result, list) and all(isinstance(item, dict) for item in result):
        return result
    
    return result




def extract_subcluster_results_with_llm(analysis_text, provider="anthropic", model="claude-3-5-sonnet-20241022", temperature=0):
    """
    Extract results with reasons from subcluster analysis text.
    
    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use (defaults to Claude 3.5 Sonnet)
        temperature: Temperature for generation (0-1)
        
    Returns:
        Extracted results in the format: results1(celltype1, celltype2, reason), results2(celltype1, celltype2, reason), etc.
        For custom APIs, may return a list of dictionaries directly.
    """
    # If the analysis_text is already in a structured format (list of dictionaries),
    # we can return it directly for processing by write_results_to_csv
    if isinstance(analysis_text, list) and all(isinstance(item, dict) for item in analysis_text):
        print(f"Analysis text is already in structured format (list of {len(analysis_text)} dictionaries)")
        return analysis_text
    
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2,reason), results2(celltype1, celltype2,reason), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    result = subcluster_agent_annotate_subcluster(prompt, provider=provider, model=model, temperature=temperature)
    
    # Check if the result is already in structured format (e.g., from DeepSeek API)
    if isinstance(result, list) and all(isinstance(item, dict) for item in result):
        return result
    
    return result



def write_results_to_csv(results, output_name='subcluster_results'):
    """
    Extract cell type results from LLM output and write to CSV file
    
    Args:
        results (str or list): LLM analysis results as string or structured data
        output_name (str): Base name for output file (will add .csv if not present)
        
    Returns:
        pandas.DataFrame: DataFrame containing the extracted results
    """
    # Add .csv suffix if not present
    if not output_name.lower().endswith('.csv'):
        output_name = output_name + '.csv'
    
    # Handle different result formats
    if isinstance(results, list):
        try:
            rows = []
            for i, item in enumerate(results, 1):
                if isinstance(item, dict):
                    # Always check both capitalized and lowercase keys
                    cluster_id = str(item.get('cluster', i))
                    # Main type
                    main_type = (
                        item.get('most_likely_top2_cell_types', ['Unknown'])[0]
                        if isinstance(item.get('most_likely_top2_cell_types'), list) and len(item.get('most_likely_top2_cell_types')) > 0
                        else item.get('main_cell_type',
                            item.get('Most likely top2 cell types', ['Unknown'])[0]
                            if isinstance(item.get('Most likely top2 cell types', ['Unknown']), list) and len(item.get('Most likely top2 cell types', ['Unknown'])) > 0
                            else item.get('Main cell type', 'Unknown')
                        )
                    )
                    # Sub type
                    sub_type = (
                        item.get('most_likely_top2_cell_types', ['Unknown', 'Unknown'])[1]
                        if isinstance(item.get('most_likely_top2_cell_types'), list) and len(item.get('most_likely_top2_cell_types')) > 1
                        else item.get('sub_cell_type',
                            item.get('Most likely top2 cell types', ['Unknown', 'Unknown'])[1]
                            if isinstance(item.get('Most likely top2 cell types', ['Unknown', 'Unknown']), list) and len(item.get('Most likely top2 cell types', ['Unknown', 'Unknown'])) > 1
                            else item.get('Sub cell type', 'Unknown')
                        )
                    )
                    # Key markers
                    key_markers = item.get('key_markers', item.get('Key marker', item.get('Key Marker', '')))
                    # Reason/explanation
                    reason = item.get('explanation', item.get('Explanation', item.get('reason', '')))
                    rows.append([cluster_id, main_type, sub_type, key_markers, reason])
            if rows:
                df = pd.DataFrame(rows, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'])
                df.to_csv(output_name, index=False)
                print(f"Results have been written to {output_name}")
                return df
        except Exception as e:
            print(f"Error processing list results: {str(e)}")
            print("Attempting to convert to string and process...")
            # Fall back to string processing
            results = str(results)
    
    # Process as string (original method)
    if isinstance(results, str):
        # Updated regex pattern to capture the reason
        pattern = r"results(\\d+)\\(([^,]+),\\s*([^)]+)\\)"
        matches = re.findall(pattern, results)

        if matches:
            # Convert matches to a DataFrame with the reason column
            df = pd.DataFrame(matches, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'])
            df.to_csv(output_name, index=False)
            print(f"Results have been written to {output_name}")
            return df
        else:
            # Try alternative pattern without reason
            alt_pattern = r"results(\\d+)\\(([^,]+),\\s*([^)]+)\\)"
            alt_matches = re.findall(alt_pattern, results)
            if alt_matches:
                # Convert matches to a DataFrame without reason column
                df = pd.DataFrame(alt_matches, columns=['Result ID', 'main_cell_type', 'sub_cell_type'])
                # Add empty reason column
                df['reason'] = ""
                df.to_csv(output_name, index=False)
                print(f"Results have been written to {output_name} (without reasons)")
                return df
    
    # If we get here, we couldn't process the results
    print(f"Warning: Could not extract results in expected format. Results type: {type(results)}")
    print(f"Saving raw results to {output_name}.txt for inspection")
    
    # Save raw results for debugging
    with open(f"{output_name}.txt", "w") as f:
        f.write(str(results))
    
    # Create a minimal DataFrame to avoid errors
    df = pd.DataFrame([["1", "Unknown", "Unknown", "Could not parse results"]], 
                     columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'])
    df.to_csv(output_name, index=False)
    
    return df



def runCASSIA_subclusters(marker, major_cluster_info, output_name, 
                       model="google/gemini-2.5-flash-preview", temperature=0, provider="openrouter", n_genes=50):
    """
    Process subclusters from marker data and generate annotated results.
    
    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        output_name: Base name for output file (will add .csv if not present)
        model: Model name to use
        temperature: Temperature parameter for API calls (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use for analysis
        
    Returns:
        None: Results are saved to a CSV file
    """
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
                                         model="google/gemini-2.5-flash-preview", temperature=0, 
                          provider="openrouter", max_workers=5, n_genes=50):
    """
    Run multiple subcluster analyses in parallel and save results.
    
    Args:
        n: Number of analyses to run
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        base_output_name: Base name for output files
        model: Model name to use
        temperature: Temperature parameter for API calls (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        max_workers: Maximum number of parallel workers
        n_genes: Number of top genes to use for analysis
        
    Returns:
        None: Results are saved to CSV files
    """
    # Force temperature to 0.3 for all runs
    temperature = 0.3
    
    def run_single_analysis(i):
        # Run the annotation process
        output_text = annotate_subclusters(marker, major_cluster_info, 
                                         model=model, temperature=temperature, provider=provider, n_genes=n_genes)
        
        # Extract results
        results = extract_subcluster_results_with_llm_multiple_output(output_text, provider=provider, model=model, temperature=temperature)
        
        # Create DataFrame based on result type
        if isinstance(results, list) and all(isinstance(item, dict) for item in results):
            # Handle structured data (list of dictionaries)
            rows = []
            for idx, item in enumerate(results, 1):
                # Try to extract data from different possible formats
                if isinstance(item, dict):
                    # Format 1: Dictionary with 'cluster' key and 'most_likely_top2_cell_types'
                    if 'cluster' in item and 'most_likely_top2_cell_types' in item:
                        cluster_id = str(item.get('cluster', idx))
                        cell_types = item.get('most_likely_top2_cell_types', ['Unknown', 'Unknown'])
                        main_type = cell_types[0] if isinstance(cell_types, list) and len(cell_types) > 0 else "Unknown"
                        sub_type = cell_types[1] if isinstance(cell_types, list) and len(cell_types) > 1 else "Unknown"
                        rows.append([cluster_id, main_type, sub_type])
                    # Format 2: Dictionary with key like 'results1'
                    elif any(key.startswith('result') for key in item.keys()):
                        key = next(k for k in item.keys() if k.startswith('result'))
                        value = item[key]
                        
                        if isinstance(value, list) and len(value) >= 2:
                            # Format: {'results1': ['celltype1', 'celltype2']}
                            main_type = value[0] if len(value) > 0 else ""
                            sub_type = value[1] if len(value) > 1 else ""
                            result_id = ''.join(filter(str.isdigit, key)) or str(idx)
                            rows.append([result_id, main_type, sub_type])
                        elif isinstance(value, dict) and 'celltype1' in value and 'celltype2' in value:
                            # Format: {'results1': {'celltype1': 'type1', 'celltype2': 'type2'}}
                            main_type = value.get('celltype1', "")
                            sub_type = value.get('celltype2', "")
                            result_id = ''.join(filter(str.isdigit, key)) or str(idx)
                            rows.append([result_id, main_type, sub_type])
                    # Format 3: Dictionary with 'key_markers', 'explanation', etc. (with underscores)
                    elif 'key_markers' in item or 'explanation' in item:
                        cluster_id = str(item.get('cluster', idx))
                        cell_types = item.get('most_likely_top2_cell_types', ['Unknown', 'Unknown'])
                        main_type = cell_types[0] if isinstance(cell_types, list) and len(cell_types) > 0 else "Unknown"
                        sub_type = cell_types[1] if isinstance(cell_types, list) and len(cell_types) > 1 else "Unknown"
                        rows.append([cluster_id, main_type, sub_type])
                    # Format 4: Dictionary with 'Key marker', 'Explanation', etc. (with spaces and capital letters)
                    elif 'Key marker' in item or 'Explanation' in item or 'Most likely top2 cell types' in item:
                        cluster_id = str(idx)  # Use index as cluster ID since no cluster field
                        cell_types = item.get('Most likely top2 cell types', ['Unknown', 'Unknown'])
                        main_type = cell_types[0] if isinstance(cell_types, list) and len(cell_types) > 0 else "Unknown"
                        sub_type = cell_types[1] if isinstance(cell_types, list) and len(cell_types) > 1 else "Unknown"
                        rows.append([cluster_id, main_type, sub_type])
            
            if rows:
                df = pd.DataFrame(rows, columns=['Cluster ID', 'main_cell_type', 'sub_cell_type'])
            else:
                # Fallback if we couldn't extract structured data
                df = pd.DataFrame([[str(i+1), "Unknown", "Unknown"] for i in range(len(results))], 
                                 columns=['Cluster ID', 'main_cell_type', 'sub_cell_type'])
        else:
            # Use regex to extract the results (original method)
            pattern = r"results(\\d+)\\(([^,]+),\\s*([^)]+)\\)"
            matches = re.findall(pattern, results)
        
            if matches:
                # Convert matches to a DataFrame
                df = pd.DataFrame(matches, columns=['Cluster ID', 'main_cell_type', 'sub_cell_type'])
            else:
                # Fallback if regex didn't match
                df = pd.DataFrame([[str(i+1), "Unknown", "Unknown"]], 
                                 columns=['Cluster ID', 'main_cell_type', 'sub_cell_type'])

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
            df['Cluster ID'] = df['Cluster ID'].astype(str)

            # Make a copy of the original values before swapping
            original_cluster_ids = df['Cluster ID'].copy()

            # Check if marker_df has at least one row and one column
            if marker_df.shape[0] > 0 and marker_df.shape[1] > 0:
                original_marker_first_col = marker_df.iloc[:, 0].copy()

                # Perform the swap safely - only if there are enough rows in both dataframes
                min_rows = min(len(df), len(marker_df))
                if min_rows > 0:
                    df.loc[:min_rows-1, 'Cluster ID'] = original_marker_first_col[:min_rows].values
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

