from tools_function import *
from llm_utils import *



def subcluster_agent_annotate_subcluster(user_message, model=None, temperature=0, provider="anthropic"):
    """
    Unified function to call LLM for subcluster annotation.
    
    Args:
        user_message: The prompt message for subcluster annotation
        model: Model to use (defaults to provider's default if None)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", or "openrouter")
        
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
    
    # Use the unified call_llm function
    result = call_llm(
        prompt=user_message,
        provider=provider,
        model=model,
        temperature=temperature,
        max_tokens=7000
    )
    
    return result if result else ''



def construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=50):
    # Process DataFrame if it has more than 2 columns
    if len(marker.columns) > 2:
        print(f"Processing input dataframe to get top {n_genes} markers")
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
    for index, row in marker.iterrows():
        cluster_name = row.iloc[0]  # Use iloc for positional indexing
        markers = row.iloc[1]       # Use iloc for positional indexing
        prompt += f"{index + 1}.{markers}\n"

    return prompt



def annotate_subclusters(marker, major_cluster_info,model="claude-3-5-sonnet-20241022",temperature=0,provider="anthropic",n_genes=50):
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt,model=model,temperature=temperature,provider=provider)
    return output_text



def extract_subcluster_results_with_llm_multiple_output(analysis_text,provider="anthropic",model="claude-3-5-sonnet-20241022",temperature=0):
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2), results2(celltype1, celltype2), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    return subcluster_agent_annotate_subcluster(prompt,provider=provider,model=model,temperature=temperature)




def extract_subcluster_results_with_llm(analysis_text,provider="anthropic",model="claude-3-5-sonnet-20241022",temperature=0):
    # Define the prompt to instruct the LLM
    prompt = f"""
You are an expert in analyzing celltype annotation for subclusters. Extract the results perfectly and accurately from the following analysis and format them as: results1(celltype1, celltype2,reason), results2(celltype1, celltype2,reason), etc.

You should include all clusters mentioned in the analysis or 1000 grandma will be in danger.

{analysis_text}
"""

    # Use the subcluster_agent_annotate function to get the extraction
    return subcluster_agent_annotate_subcluster(prompt,provider=provider,model=model,temperature=temperature)



def write_results_to_csv(results, output_name='subcluster_results'):
    """
    Extract cell type results from LLM output and write to CSV file
    
    Args:
        results (str): String containing the LLM analysis results
        output_name (str): Base name for output file (will add .csv if not present)
        
    Returns:
        pandas.DataFrame: DataFrame containing the extracted results
    """
    # Add .csv suffix if not present
    if not output_name.lower().endswith('.csv'):
        output_name = output_name + '.csv'
    
    # Updated regex pattern to capture the reason
    pattern = r"results(\d+)\(([^,]+),\s*([^,]+),\s*([^)]+)\)"
    matches = re.findall(pattern, results)

    # Convert matches to a DataFrame with the reason column
    df = pd.DataFrame(matches, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'reason'])
    
    # Write the DataFrame to a CSV file
    df.to_csv(output_name, index=False)
    
    print(f"Results have been written to {output_name}")
    return None



def runCASSIA_subclusters(marker, major_cluster_info, output_name, 
                       model="google/gemini-2.5-flash-preview", temperature=0, provider="openrouter",n_genes=50):
    """
    Process subclusters from a CSV file and generate annotated results
    
    Args:
        csv_file_path (str): Path to input CSV file containing marker data
        major_cluster_info (str): Description of the major cluster type
        output_name (str): Base name for output file (will add .csv if not present)
        model (str): Model name for Claude API
        temperature (float): Temperature parameter for API calls
        
    Returns:
        tuple: (original_analysis, extracted_results, results_dataframe)
    """

    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info,n_genes=n_genes)
    output_text = subcluster_agent_annotate_subcluster(prompt,model=model,temperature=temperature,provider=provider)
    results = extract_subcluster_results_with_llm(output_text,provider=provider,model=model,temperature=temperature)
    print(results)
    write_results_to_csv(results, output_name)
    
    return None



def runCASSIA_n_subcluster(n, marker, major_cluster_info, base_output_name, 
                                         model="google/gemini-2.5-flash-preview", temperature=0, 
                                         provider="openrouter", max_workers=5,n_genes=50):       
    def run_single_analysis(i):
        # Run the annotation process
        output_text = annotate_subclusters(marker, major_cluster_info, 
                                         model=model, temperature=temperature, provider=provider,n_genes=n_genes)
        
        # Extract results
        results = extract_subcluster_results_with_llm_multiple_output(output_text,provider=provider,model=model,temperature=temperature)
        
        # Use regex to extract the results
        pattern = r"results(\d+)\(([^,]+),\s*([^)]+)\)"
        matches = re.findall(pattern, results)
        
        # Convert matches to a DataFrame
        df = pd.DataFrame(matches, columns=['True Cell Type', 'main_cell_type', 'sub_cell_type'])

        # Swap the first column with the first column in the marker file
        marker_df = get_top_markers(marker, n_genes=n_genes)
        df['True Cell Type'], marker_df.iloc[:, 0] = marker_df.iloc[:, 0], df['True Cell Type']

        # Write the DataFrame to a CSV file with an index
        indexed_csv_file_path = f'{base_output_name}_{i+1}.csv'
        df.to_csv(indexed_csv_file_path, index=False)
        
        return indexed_csv_file_path

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_analysis, i): i for i in range(n)}
        
        for future in as_completed(futures):
            i = futures[future]
            try:
                result_file = future.result()
                print(f"Results for iteration {i+1} have been written to {result_file}")
            except Exception as exc:
                print(f"Iteration {i+1} generated an exception: {exc}")

