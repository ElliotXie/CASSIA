from typing import List, Dict, Any, Optional
import re
import csv
import os
import concurrent.futures
import json
from datetime import datetime
try:
    from CASSIA.core.llm_utils import call_llm
    from CASSIA.reports.generate_hypothesis_report import create_html_report
    from .summarize_hypothesis_runs import summarize_runs
except ImportError:
    try:
        from .llm_utils import call_llm
        from .generate_hypothesis_report import create_html_report
        from .summarize_hypothesis_runs import summarize_runs
    except ImportError:
        # Fallback for running the script directly for testing
        from llm_utils import call_llm
        from generate_hypothesis_report import create_html_report
        from summarize_hypothesis_runs import summarize_runs

def _create_candidate_generation_prompt(marker_genes: List[str], species: str, tissue: str) -> str:
    """Creates a prompt for the first LLM to generate candidate cell types."""
    marker_str = ", ".join(marker_genes[:50])  # Use top 50 markers
    prompt = (
        f"You are a professional biologist analyzing a sample from {species} {tissue}. "
        "Based on the following list of top 50 marker genes, "
        "what could this cell cluster be? Give me the top 3 most likely detailed subtypes."
        f"\n\nMarker genes: {marker_str}"
    )
    return prompt

def _create_summarization_prompt(candidate_text: str) -> str:
    """Creates a prompt for the second LLM to parse and format the candidate text."""
    prompt = (
        "You are an expert data formatter. Your task is to parse the following text, which contains cell type hypotheses, "
        "and structure it into a ranked list. Extract the top 3 cell types, their rank, and the reasoning provided.\n\n"
        f"Raw text to parse:\n---\n{candidate_text}\n---\n\n"
        "Please provide your answer in the following strict format, listing the most likely cell type first:\n\n"
        "1. **Cell Type:** [Detailed cell type name]\n"
        "   **Reasoning:** [Brief explanation based on the provided text]\n\n"
        "2. **Cell Type:** [Detailed cell type name]\n"
        "   **Reasoning:** [Brief explanation based on the provided text]\n\n"
        "3. **Cell Type:** [Detailed cell type name]\n"
        "   **Reasoning:** [Brief explanation based on the provided text]"
    )
    return prompt

def _parse_llm_response(response: str) -> List[Dict[str, Any]]:
    """Parses the structured LLM response to extract cell type hypotheses."""
    parsed_results = []
    # This pattern is designed to match the structured output from the formatter LLM.
    pattern = re.compile(
        r"^\s*(\d+)\.\s*\*\*Cell Type:\*\*\s*(.*?)\s*\*\*Reasoning:\*\*\s*(.*?)(?=\n\s*\d+\.|\Z)",
        re.DOTALL | re.MULTILINE
    )
    matches = pattern.findall(response)
    
    for match in matches:
        try:
            rank = int(match[0])
            cell_type = match[1].strip()
            reasoning = match[2].strip()
            parsed_results.append({
                "rank": rank,
                "cell_type": cell_type,
                "reasoning": reasoning
            })
        except (ValueError, IndexError):
            # Skip any malformed entries
            continue
        
    return parsed_results

def generate_hypothesis(
    marker_genes: List[str],
    species: str,
    tissue: str,
    generator_provider: str = "openrouter",
    formatter_provider: str = "openrouter",
    generator_api_key: Optional[str] = None,
    formatter_api_key: Optional[str] = None,
    generator_model: Optional[str] = "openai/gpt-4o-2024-11-20",    #google/gemini-2.5-flash, openai/gpt-4o
    formatter_model: Optional[str] = "google/gemini-2.5-flash",
    generator_temp: float = 0.5,
    formatter_temp: float = 0.0,
    max_tokens: int = 4096,
    **kwargs
) -> Dict[str, Any]:
    """
    Generates cell type hypotheses using a two-step LLM pipeline.
    
    1. Generator LLM: Brainstorms candidate cell types.
    2. Formatter LLM: Parses the brainstormed text into a structured format.

    Args:
        marker_genes: A list of marker gene symbols for a cluster.
        species: The species of the sample (e.g., "human").
        tissue: The tissue of the sample (e.g., "bone marrow").
        generator_provider: The LLM provider for the generator step (e.g., "openrouter", "openai").
        formatter_provider: The LLM provider for the formatter step.
        generator_api_key: API key for the generator provider.
        formatter_api_key: API key for the formatter provider.
        generator_model: The model for the creative generation step.
        formatter_model: The model for the structuring/formatting step.
        generator_temp: Temperature for the generator model (higher for more creativity).
        formatter_temp: Temperature for the formatter model (lower for more determinism).
        max_tokens: Maximum tokens for the LLM responses.
        **kwargs: Additional parameters for the llm_utils.call_llm function.

    Returns:
        A dictionary containing raw and parsed responses.
    """
    # Step 1: Generate candidate cell types with the first LLM
    candidate_prompt = _create_candidate_generation_prompt(marker_genes, species, tissue)
    
    candidate_response_text = call_llm(
        prompt=candidate_prompt,
        provider=generator_provider,
        model=generator_model,
        api_key=generator_api_key,
        temperature=generator_temp,
        max_tokens=max_tokens,
        system_prompt="You are a professional biologist assisting with single-cell data analysis.",
        additional_params=kwargs
    )
    
    print("\n--- Raw Response from Generator Agent ---")
    try:
        print(candidate_response_text)
    except UnicodeEncodeError:
        print("[Note: Raw response contains characters that cannot be displayed in this terminal. The full response is saved in the JSON and HTML report.]")
    print("---------------------------------------\n")
    
    # Step 2: Format the candidates with the second LLM
    summarization_prompt = _create_summarization_prompt(candidate_response_text)
    
    formatted_response_text = call_llm(
        prompt=summarization_prompt,
        provider=formatter_provider,
        model=formatter_model,
        api_key=formatter_api_key,
        temperature=formatter_temp,
        max_tokens=max_tokens,
        system_prompt="You are an expert data formatter. Your task is to extract and structure information precisely as requested.",
        additional_params=kwargs
    )
    
    # Step 3: Parse the final structured response
    parsed_response = _parse_llm_response(formatted_response_text)
    
    return {
        "raw_response": candidate_response_text,
        "parsed_response": parsed_response
    }

def process_marker_file(
    file_path: str,
    species: str,
    tissue: str,
    output_json_path: str,
    generator_model: str,
    formatter_model: str,
    **kwargs
) -> None:
    """
    Processes a CSV or TSV file, generates hypotheses, and saves results to a JSON file.

    Args:
        file_path: The path to the marker file.
        species: The species of the sample.
        tissue: The tissue of the sample.
        output_json_path: Path to save the output JSON file.
        generator_model: The model for the creative generation step.
        formatter_model: The model for the structuring/formatting step.
        **kwargs: Additional arguments to pass to the generate_hypothesis function.
    """
    all_hypotheses = {}
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        return

    tasks_to_process = []
    with open(file_path, mode='r', encoding='utf-8') as infile:
        try:
            dialect = csv.Sniffer().sniff(infile.read(2048))
            infile.seek(0)
        except csv.Error:
            dialect = 'excel'
            infile.seek(0)
            
        reader = csv.reader(infile, dialect)
        
        try:
            first_row = next(reader)
            if ',' not in first_row[1]:
                pass
            else:
                infile.seek(0)
                reader = csv.reader(infile, dialect)
        except StopIteration:
            return

        for row in reader:
            if not row or len(row) < 2:
                continue
            
            cluster_name = row[0].strip()
            marker_genes = [gene.strip() for gene in row[1].split(',') if gene.strip()]
            
            if not marker_genes:
                continue
            
            tasks_to_process.append({'cluster_name': cluster_name, 'marker_genes': marker_genes})

    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_cluster = {}
        for task in tasks_to_process:
            cluster_name = task['cluster_name']
            marker_genes = task['marker_genes']
            print(f"--- Submitting cluster for processing: {cluster_name} ---")
            future = executor.submit(generate_hypothesis, marker_genes, species, tissue, **kwargs)
            future_to_cluster[future] = cluster_name
            
        for future in concurrent.futures.as_completed(future_to_cluster):
            cluster_name = future_to_cluster[future]
            try:
                data = future.result()
                all_hypotheses[cluster_name] = data
                print(f"--- Successfully processed cluster: {cluster_name} ---")
            except Exception as exc:
                print(f'--- Cluster {cluster_name} generated an exception: {exc} ---')
                all_hypotheses[cluster_name] = {
                    "raw_response": None,
                    "parsed_response": [{"error": str(exc)}]
                }

    # Save the consolidated results to a JSON file
    with open(output_json_path, 'w', encoding='utf-8') as f:
        json.dump(all_hypotheses, f, indent=4)
    print(f"\n--- All results for this run saved to {output_json_path} ---")

def run_multi_analysis(
    marker_file_path: str,
    species: str,
    tissue: str,
    num_runs: int,
    generator_model_name: str,
    formatter_model_name: str
):
    """
    Runs the full analysis pipeline multiple times and generates individual reports.
    It concludes by creating a manifest file for the summarization step.
    """
    if not os.path.exists(marker_file_path):
        print(f"Error: Marker file not found at {marker_file_path}")
        return

    script_dir = os.path.dirname(os.path.abspath(__file__))
    generated_json_files = []
    run_timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    for i in range(num_runs):
        print(f"\n\n--- Starting Run {i+1} of {num_runs} ---")
        
        # --- Generate dynamic filename for this run ---
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        safe_species = species.replace(" ", "_")
        safe_tissue = tissue.replace(" ", "_")
        safe_gen_model = generator_model_name.replace("/", "_")
        
        base_filename = f"{timestamp}_run{i+1}_{safe_species}_{safe_tissue}_{safe_gen_model}"
        # Ensure all output is in the script's directory
        output_json_path = os.path.join(script_dir, f"{base_filename}.json")
        output_html_path = os.path.join(script_dir, f"{base_filename}.html")
        generated_json_files.append(output_json_path)
        # ---

        process_marker_file(
            file_path=marker_file_path,
            species=species,
            tissue=tissue,
            output_json_path=output_json_path,
            generator_model=generator_model_name,
            formatter_model=formatter_model_name
        )
        
        # Automatically generate the HTML report for this run
        print(f"\n--- Automatically generating HTML report for Run {i+1} ---")
        try:
            create_html_report(output_json_path, output_html_path)
        except Exception as e:
            print(f"Could not generate HTML report for Run {i+1}. Error: {e}")
    
    # --- Instructions for Final Consolidation ---
    if num_runs > 1 and generated_json_files:
        manifest_filename = os.path.join(script_dir, f"run_manifest_{run_timestamp}.csv")
        with open(manifest_filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['json_file_path'])
            for fname in generated_json_files:
                writer.writerow([fname])
        
        print("\n\n--- All runs complete. ---")
        print(f"A manifest file '{manifest_filename}' has been created with the paths to all generated JSON files.")
        print("\nTo generate a consolidated summary report, run the following command in your terminal:")
        print(f"\npython CASSIA/CASSIA_python/CASSIA/summarize_hypothesis_runs.py {manifest_filename}\n")

if __name__ == '__main__':
    """
    Main entry point to run the hypothesis generation pipeline.
    """
    # --- CONTEXT & PARAMETERS ---
    marker_file_path = r"C:\\Users\\ellio\\OneDrive - UW-Madison\\Revision_cassia\\3_csscore\\bone_top_50_genes (1).csv"
    species = "human"
    tissue = "bone marrow"
    num_runs = 3 # Set to 1 to disable multi-run and manifest generation.
    generator_model_name = "google/gemini-2.5-flash"
    formatter_model_name = "google/gemini-2.5-flash"
    # --------------------------

    run_multi_analysis(
        marker_file_path=marker_file_path,
        species=species,
        tissue=tissue,
        num_runs=num_runs,
        generator_model_name=generator_model_name,
        formatter_model_name=formatter_model_name
    )

    print("\n--- Main Analysis Script Finished ---")
