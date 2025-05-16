import os
import re
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Any, Optional, Union

from .llm_utils import call_llm

def prompt_hypothesis_generator(major_cluster_info: str, comma_separated_genes: str, annotation_history: str) -> str:
    """
    Generate a prompt for iterative marker analysis without additional tasks.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful senior computational biologist called in whenever an annotation needs deeper scrutiny, disambiguation, or simply a second opinion. Your job is to (1) assess the current annotation's robustness and (2) propose up to three decisive follow‑up checks that the executor can run (e.g., examine expression of key positive or negative markers). You should do a good job or 10 grandma are going to be in danger. You never rush to conclusions and are always careful

Context Provided to You

Cluster summary：{major_cluster_info}

Top ranked markers (high → low)：
 {comma_separated_genes}

Prior annotation dialogue：
 {annotation_history}

What to Do
1. Brief Evaluation – One concise paragraph that:

    - Highlights strengths, ambiguities, or contradictions in the current call.

    - Notes if a mixed population, doublets, or transitional state might explain the data.

2. Design up to 3 follow‑up checks (cell types or biological hypotheses):

    - Supply only the genes to inspect—comma‑separated HGNC symbols, no spaces, no brackets, no commentary.

    - Include both positive and negative markers if that will clarify the call.

    - Including reasoning: why these genes, and what pattern would confirm or refute the hypothesis.

3. Upon receiving gene expression results, refine the hypothesis or generate new ones. Continue Step 2 iteratively until the cluster is confidently and fully annotated. Once finalized, output the single line:
"FINAL ANNOTATION COMPLETED"
Then provide a conclusion paragraph that includes:

1.The final cell type
2.Confidence level (high, medium, or low)
3.Key markers supporting your conclusion
4.Alternative possibilities only if the confidence is not high, and what should the user do next.



Output Template：

Evaluation
[One short paragraph]

celltype to check 1

<check_genes>
GENE1, GENE2, GENE3
</check_genes>

<reasoning>
Why these genes and what we expect to see.
</reasoning>

celltype to check 2

<check_genes>
…
</check_genes>
<reasoning>
…
</reasoning>

hypothesis to check 3
<check_genes>
…
</check_genes>

<reasoning>
…
</reasoning>

*Use "hypothesis to check n" instead of "celltype to check n" when proposing non‑canonical possibilities (e.g., "cycling subpopulation", "doublet").
*Provide no more than three total blocks (celltype + hypothesis combined).


Tone & Style Guidelines

Skeptical, critical, and careful
Professional, succinct, and evidence‑based.
Reference established biology when helpful ("LYZ vs. CTSW distinguishes myeloid from T‑cell lineages").


"""
    return prompt

def prompt_hypothesis_generator_additional_task(major_cluster_info: str, comma_separated_genes: str, annotation_history: str, additional_task: str) -> str:
    """
    Generate a prompt for iterative marker analysis with an additional task.
    
    Args:
        major_cluster_info: Information about the cluster being analyzed
        comma_separated_genes: Comma-separated list of marker genes
        annotation_history: Previous annotation history
        additional_task: Additional analysis task to perform
        
    Returns:
        str: Generated prompt text
    """
    prompt = f"""You are a careful professional biologist, specializing in single-cell RNA-seq analysis. 
I'll provide you with some genes (comma-separated) from a cell type in {major_cluster_info} and I want you to help identify the cell type.

Here are the marker genes that are differentially expressed in this cluster:
{comma_separated_genes}

Previous annotation/analysis:
{annotation_history if annotation_history else "No previous annotation available."}

Please follow these steps carefully:
1. Based on the marker genes, generate hypotheses about what cell type this might be.
2. For each hypothesis, provide supporting evidence from the marker genes.
3. To validate your hypotheses, list additional marker genes that should be checked using this format:
<check_genes>gene1, gene2, gene3, etc.</check_genes>
4. After I provide additional gene information, refine your hypotheses.
5. Continue this process until you reach a confident conclusion.
6. {additional_task} - perform this analysis based on the marker gene expression and patterns.
7. When you are ready to make a final determination, state: "FINAL ANALYSIS COMPLETED" followed by your conclusion.

Your final analysis should include:
- General cell type
- Specific cell subtype (if applicable)
- Confidence level (high, medium, low)
- Alternative possibilities if confidence is not high
- Key markers supporting your conclusion
- Analysis of the additional task: {additional_task}

Please start by analyzing the provided markers and forming initial hypotheses.
"""
    return prompt

def get_marker_info(gene_list: List[str], marker: Union[pd.DataFrame, Any]) -> str:
    """
    Extract information about marker genes from a marker dataset.
    
    Args:
        gene_list: List of gene names to filter the marker dataset
        marker: DataFrame containing marker gene expression data
        
    Returns:
        str: Formatted string with marker information
    """
    def filter_marker(gene_names: List[str]) -> Tuple[pd.DataFrame, List[str]]:
        # Convert marker to pandas DataFrame if it's not already
        if not isinstance(marker, pd.DataFrame):
            marker_df = pd.DataFrame(marker)
        else:
            marker_df = marker.copy()
        
        # Remove any 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in marker_df.columns:
            marker_df = marker_df.drop(columns=['Unnamed: 0'])

        # Identify valid genes and NA genes
        valid_genes = []
        na_genes = []
        
        # Check for 'gene' column first to improve gene lookup
        if 'gene' in marker_df.columns:
            for gene in gene_names:
                # Try to find exact match in the gene column
                gene_rows = marker_df[marker_df['gene'].str.contains(f"^{gene}$", case=True, regex=True, na=False)]
                if not gene_rows.empty:
                    row_idx = gene_rows.index[0]
                    valid_genes.append(gene)
                else:
                    # Try in index
                    if gene in marker_df.index:
                        valid_genes.append(gene)
                    else:
                        na_genes.append(gene)
        else:
            # Original method - look up in index
            for gene in gene_names:
                if gene in marker_df.index:
                    # Check if all values for this gene are NA
                    gene_data = marker_df.loc[gene]
                    if gene_data.isna().all() or (gene_data == 'NA').all():
                        na_genes.append(gene)
                    else:
                        valid_genes.append(gene)
                else:
                    na_genes.append(gene)
        
        # Create result DataFrame with only valid genes
        if valid_genes:
            result = marker_df.loc[valid_genes].copy()
        else:
            # If no valid genes, create an empty dataframe with the same columns
            result = pd.DataFrame(columns=marker_df.columns)
            
        # Only try to format numeric columns that exist
        numeric_cols = result.select_dtypes(include=[np.number]).columns
        for col in numeric_cols:
            try:
                result[col] = result[col].apply(lambda x: f"{float(x):.2e}" if pd.notnull(x) and x != 'NA' else x)
            except:
                continue

        return result.iloc[:, 0:5], na_genes

    # Filter to rows based on gene name and get NA genes list
    marker_filtered, na_genes = filter_marker(gene_list)
    
    # Generate marker info string from valid genes only
    marker_string = marker_filtered.to_string()
    
    # If there are genes with all NA values, add a message
    if na_genes:
        na_genes_message = f"\nNote: The following genes are not in the differential expression list: {', '.join(na_genes)}"
        marker_string += na_genes_message

    return marker_string

def extract_genes_from_conversation(conversation: str) -> List[str]:
    """
    Extract gene lists from conversation using the check_genes tag.
    
    Args:
        conversation: Text containing gene lists in check_genes tags
        
    Returns:
        List[str]: List of unique gene names
    """
    # Extract gene lists and get marker info
    gene_lists = re.findall(r'<check_genes>\s*(.*?)\s*</check_genes>', conversation, re.DOTALL)
    
    # Improve gene extraction to handle special cases
    all_genes = []
    for gene_list in gene_lists:
        # Clean and normalize the gene list
        # Replace common separators with commas
        cleaned_list = re.sub(r'[\]\[\)\(]', '', gene_list)
        cleaned_list = re.sub(r'\s+', ' ', cleaned_list)
        # Split by comma or space, depending on formatting
        genes = re.split(r',\s*|\s+', cleaned_list)
        all_genes.extend([g.strip() for g in genes if g.strip()])
    
    return sorted(set(all_genes))

def iterative_marker_analysis(
    major_cluster_info: str, 
    marker: Union[pd.DataFrame, Any], 
    comma_separated_genes: str, 
    annotation_history: str, 
    num_iterations: int = 2, 
    provider: str = "openrouter",
    model: Optional[str] = None,
    additional_task: Optional[str] = None,
    temperature: float = 0
) -> Tuple[str, List[Dict[str, str]]]:
    """
    Perform iterative marker analysis using the specified LLM provider.
    
    Args:
        major_cluster_info: Information about the cluster
        marker: DataFrame or other structure containing marker gene expression data
        comma_separated_genes: List of genes as comma-separated string
        annotation_history: Previous annotation history
        num_iterations: Maximum number of iterations
        provider: LLM provider to use ('openai', 'anthropic', or 'openrouter')
        model: Specific model from the provider to use
        additional_task: Optional additional task to perform during analysis
        temperature: Sampling temperature (0-1)
        
    Returns:
        tuple: (final_response_text, messages)
    """
    # Select the appropriate prompt based on whether there's an additional task
    if additional_task:
        prompt = prompt_hypothesis_generator_additional_task(
            major_cluster_info, 
            comma_separated_genes, 
            annotation_history, 
            additional_task
        )
        completion_marker = "FINAL ANALYSIS COMPLETED"
    else:
        prompt = prompt_hypothesis_generator(
            major_cluster_info, 
            comma_separated_genes, 
            annotation_history
        )
        completion_marker = "FINAL ANNOTATION COMPLETED"
    
    # Initialize the conversation history
    messages = [{"role": "user", "content": prompt}]
    
    # Iterative process
    for iteration in range(num_iterations):
        try:
            # Call the LLM
            conversation = call_llm(
                prompt=messages[-1]["content"],
                provider=provider,
                model=model,
                temperature=temperature,
                max_tokens=7000,
                # If not the first message, include conversation history
                additional_params={"messages": messages} if iteration > 0 else {}
            )
            
            # Check if the analysis is complete
            if completion_marker in conversation:
                print(f"Final annotation completed in iteration {iteration + 1}.")
                messages.append({"role": "assistant", "content": conversation})
                return conversation, messages
            
            # Extract gene lists and get marker info
            unique_genes = extract_genes_from_conversation(conversation)
            
            if unique_genes:
                # Get marker information for the requested genes
                retrieved_marker_info = get_marker_info(unique_genes, marker)
                
                # Append messages
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": retrieved_marker_info})
                
                print(f"Iteration {iteration + 1} completed.")
            else:
                # No genes to check, simply continue the conversation
                messages.append({"role": "assistant", "content": conversation})
                messages.append({"role": "user", "content": "Please continue your analysis and provide a final annotation."})
                
                print(f"Iteration {iteration + 1} completed (no genes to check).")
        
        except Exception as e:
            print(f"Error in iteration {iteration + 1}: {str(e)}")
            return f"Error occurred: {str(e)}", messages
    
    # Final response if max iterations reached
    try:
        final_response = call_llm(
            prompt="Please provide your final analysis based on all the information so far.",
            provider=provider,
            model=model,
            temperature=temperature,
            max_tokens=7000,
            additional_params={"messages": messages}
        )
        
        messages.append({"role": "assistant", "content": final_response})
        return final_response, messages
            
    except Exception as e:
        print(f"Error in final response: {str(e)}")
        return f"Error in final response: {str(e)}", messages

def prepare_analysis_data(full_result_path: str, marker_path: str, cluster_name: str, conversation_history_mode: str = "final") -> Tuple[pd.DataFrame, pd.DataFrame, str, str]:
    """
    Load and prepare data for marker analysis.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker_path: Path to the marker genes CSV file
        cluster_name: Name of the cluster to analyze
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
            - "full": Use the entire conversation history
            - "final": Extract only the part between "Step 6" and "FINAL ANNOTATION COMPLETED" (default)
            - "none": Don't include any conversation history
        
    Returns:
        tuple: (full_results, marker_data, top_markers_string, annotation_history)
    """
    # Load the full results - try first without index_col, then with index_col=0
    try:
        full_results = pd.read_csv(full_result_path)
    except:
        full_results = pd.read_csv(full_result_path, index_col=0)
    
    # Load marker data - handle both DataFrame and file path
    if isinstance(marker_path, pd.DataFrame):
        marker = marker_path
    else:
        try:
            marker = pd.read_csv(marker_path)
        except:
            marker = pd.read_csv(marker_path, index_col=0)
        
        # Clean up by removing the 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in marker.columns:
            marker.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Try to find the cluster column name - the old code assumed 'cluster'
    cluster_column = 'True Cell Type'  # Default based on the CSV file we checked
    if cluster_column not in full_results.columns and 'cluster' in full_results.columns:
        cluster_column = 'cluster'
    
    # Filter the results for the specified cluster
    cluster_data = full_results[full_results[cluster_column] == cluster_name]
    
    if cluster_data.empty:
        raise ValueError(f"No data found for cluster '{cluster_name}' in the full results file. Available clusters: {full_results[cluster_column].unique().tolist()}")
    
    # Get marker column from the CSV if it exists, otherwise use 'Marker List'
    marker_column = 'Marker List' if 'Marker List' in cluster_data.columns else None
    
    if marker_column and not cluster_data[marker_column].empty and not pd.isna(cluster_data[marker_column].iloc[0]):
        # Use the markers directly from the CSV
        top_markers_string = cluster_data[marker_column].iloc[0]
    else:
        # Use a default marker list or generate based on the marker data
        # This is a fallback in case markers aren't in the CSV
        top_markers_string = "CD14, CD11B, CD68, CSF1R, CX3CR1, CD163, MSR1, ITGAM, FCGR1A, CCR2"
        print(f"Warning: No marker list found for {cluster_name}, using default markers")
    
    # Extract conversation history if available
    annotation_history = ""
    if 'Conversation History' in cluster_data.columns and conversation_history_mode != "none":
        try:
            # Get the conversation history from the first row (should be the same for all rows of this cluster)
            if not cluster_data['Conversation History'].empty and not pd.isna(cluster_data['Conversation History'].iloc[0]):
                full_history = cluster_data['Conversation History'].iloc[0]
                
                # Process the conversation history based on the selected mode
                if conversation_history_mode == "full":
                    annotation_history = full_history
                elif conversation_history_mode == "final":
                    # First find the Step 6 line with various possible formats
                    step6_pattern = r'(?:\*\*Step 6:|Step 6:|STEP 6:|step 6:|Step6:).*?Concise Summary.*?\*\*'
                    # Then extract everything between that line and FINAL ANNOTATION COMPLETED
                    # but exclude the Step 6 line itself
                    match = re.search(f'{step6_pattern}(.*?)(?=FINAL ANNOTATION COMPLETED)', full_history, re.DOTALL)
                    
                    if match:
                        # Use group(1) to get only the content after the Step 6 line
                        annotation_history = match.group(1).strip()
                        print(f"Extracted final section of conversation history ({len(annotation_history)} characters)")
                    else:
                        # Fallback to original pattern if the specific format isn't found
                        step6_match = re.search(r'(?:Step 6|STEP 6|step 6|Step6).*?(?=FINAL ANNOTATION COMPLETED)', full_history, re.DOTALL)
                        if step6_match:
                            annotation_history = step6_match.group(0).strip()
                            print(f"Using original pattern - extracted final section ({len(annotation_history)} characters)")
                        else:
                            # If Step 6 pattern not found, use the full history
                            annotation_history = full_history
                            print(f"Step 6 pattern not found in conversation history, using full history")
                else:
                    # For any unrecognized mode, default to full history
                    annotation_history = full_history
                
                print(f"Using conversation history for cluster {cluster_name} (mode: {conversation_history_mode}, {len(annotation_history)} characters)")
            else:
                print(f"Conversation History column exists but is empty for cluster {cluster_name}")
        except Exception as e:
            print(f"Error extracting conversation history: {str(e)}")
    else:
        if conversation_history_mode == "none":
            print(f"Note: Conversation history extraction disabled (mode: none)")
        else:
            print(f"Note: 'Conversation History' column not found in the results file")
    
    return full_results, marker, top_markers_string, annotation_history

def generate_html_report(analysis_text: str) -> str:
    """
    Generate an HTML report from the analysis text.
    
    Args:
        analysis_text: The text output from the analysis
        
    Returns:
        str: HTML formatted report
    """
    # Split the text into sections
    sections = analysis_text.split(" | ")
    
    # Simple HTML template
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{ 
                font-family: Arial, sans-serif;
                max-width: 1000px;
                margin: 0 auto;
                padding: 20px;
                line-height: 1.6;
            }}
            .container {{ 
                background-color: #f8f9fa;
                padding: 30px;
                border-radius: 8px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.1);
            }}
            h1, h2, h3 {{ color: #333; }}
            pre {{ 
                background-color: #f5f5f5;
                padding: 15px;
                border-radius: 4px;
                overflow-x: auto;
            }}
            .highlight {{ 
                background-color: #e6f7ff;
                border-left: 4px solid #1890ff;
                padding: 10px 15px;
                margin: 20px 0;
            }}
        </style>
        <title>Cell Type Analysis Report</title>
    </head>
    <body>
        <div class="container">
            <h1>Cell Type Analysis Report</h1>
            <div class="content">
                {analysis_text.replace("\n", "<br>")}
            </div>
        </div>
    </body>
    </html>
    """
    
    return html

def save_html_report(report: str, filename: str) -> None:
    """
    Save the HTML report to a file.
    
    Args:
        report: HTML report content
        filename: Output filename
    """
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Report saved to {filename}")

def runCASSIA_annotationboost(
    full_result_path: str,
    marker: Union[str, pd.DataFrame],
    cluster_name: str,
    major_cluster_info: str,
    output_name: str,
    num_iterations: int = 5,
    model: Optional[str] = None,
    provider: str = "openrouter",
    temperature: float = 0,
    conversation_history_mode: str = "final"
) -> Tuple[str, List[Dict[str, str]]]:
    """
    Run annotation boost analysis for a given cluster.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker: Path to marker genes CSV file or DataFrame with marker data
        cluster_name: Name of the cluster to analyze
        major_cluster_info: General information about the dataset (e.g., "Human PBMC")
        output_name: Base name for the output HTML file
        num_iterations: Number of iterations for marker analysis (default=5)
        model: Model to use for analysis - if None, uses the provider's default
        provider: AI provider to use ('openai', 'anthropic', or 'openrouter')
        temperature: Sampling temperature (0-1)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    
    Returns:
        tuple: (analysis_result, messages_history)
    """
    try:
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter']:
            raise ValueError("Provider must be one of: 'openai', 'anthropic', or 'openrouter'")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode
        )
        
        # Run the iterative marker analysis
        analysis_text, messages = iterative_marker_analysis(
            major_cluster_info=major_cluster_info,
            marker=marker_data,
            comma_separated_genes=top_markers_string,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            provider=provider,
            model=model,
            temperature=temperature
        )
        
        # Generate and save the HTML report
        html_report = generate_html_report(analysis_text)
        
        # Ensure output_name has .html extension
        if not output_name.lower().endswith('.html'):
            output_name += '.html'
        
        save_html_report(html_report, output_name)
        
        return analysis_text, messages
    
    except Exception as e:
        error_msg = f"Error in runCASSIA_annotationboost: {str(e)}"
        print(error_msg)
        raise

def runCASSIA_annotationboost_additional_task(
    full_result_path: str,
    marker: Union[str, pd.DataFrame],
    cluster_name: str,
    major_cluster_info: str,
    output_name: str,
    num_iterations: int = 5,
    model: Optional[str] = None,
    provider: str = "openrouter",
    additional_task: str = "check if this is a cancer cluster",
    temperature: float = 0,
    conversation_history_mode: str = "final"
) -> Tuple[str, List[Dict[str, str]]]:
    """
    Run annotation boost analysis with an additional task for a given cluster.
    
    Args:
        full_result_path: Path to the full results CSV file
        marker: Path to marker genes CSV file or DataFrame with marker data
        cluster_name: Name of the cluster to analyze
        major_cluster_info: General information about the dataset (e.g., "Human PBMC")
        output_name: Base name for the output HTML file
        num_iterations: Number of iterations for marker analysis (default=5)
        model: Model to use for analysis - if None, uses the provider's default
        provider: AI provider to use ('openai', 'anthropic', or 'openrouter')
        additional_task: Additional task to perform during analysis
        temperature: Sampling temperature (0-1)
        conversation_history_mode: Mode for extracting conversation history ("full", "final", or "none")
    
    Returns:
        tuple: (analysis_result, messages_history)
    """
    try:
        # Validate provider input
        if provider.lower() not in ['openai', 'anthropic', 'openrouter']:
            raise ValueError("Provider must be one of: 'openai', 'anthropic', or 'openrouter'")
        
        # Prepare the data
        _, marker_data, top_markers_string, annotation_history = prepare_analysis_data(
            full_result_path, marker, cluster_name, conversation_history_mode
        )
        
        # Run the iterative marker analysis with additional task
        analysis_text, messages = iterative_marker_analysis(
            major_cluster_info=major_cluster_info,
            marker=marker_data,
            comma_separated_genes=top_markers_string,
            annotation_history=annotation_history,
            num_iterations=num_iterations,
            provider=provider,
            model=model,
            additional_task=additional_task,
            temperature=temperature
        )
        
        # Generate and save the HTML report
        html_report = generate_html_report(analysis_text)
        
        # Ensure output_name has .html extension
        if not output_name.lower().endswith('.html'):
            output_name += '.html'
        
        save_html_report(html_report, output_name)
        
        return analysis_text, messages
    
    except Exception as e:
        error_msg = f"Error in runCASSIA_annotationboost_additional_task: {str(e)}"
        print(error_msg)
        raise 