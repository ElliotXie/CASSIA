import json
import argparse
import os
import csv
import concurrent.futures
from datetime import datetime
from typing import List, Dict, Any

try:
    from CASSIA.core.llm_utils import call_llm
    from CASSIA.reports.generate_hypothesis_report import create_html_report, HTML_TEMPLATE, CLUSTER_TEMPLATE
except ImportError:
    try:
        from .llm_utils import call_llm
        from .generate_hypothesis_report import create_html_report, HTML_TEMPLATE, CLUSTER_TEMPLATE
    except ImportError:
        from llm_utils import call_llm
        from generate_hypothesis_report import create_html_report, HTML_TEMPLATE, CLUSTER_TEMPLATE

def _create_consolidation_prompt(cluster_name: str, all_run_results: List[Dict]) -> str:
    """Creates a prompt to consolidate results from multiple runs for a single cluster."""
    
    formatted_runs = ""
    for i, run_result in enumerate(all_run_results):
        formatted_runs += f"### Run {i+1} Hypotheses:\n"
        parsed_response = run_result.get('parsed_response', [])
        if parsed_response:
            for item in parsed_response:
                rank = item.get('rank')
                cell_type = item.get('cell_type')
                reasoning = item.get('reasoning')
                formatted_runs += f"- Rank {rank}: {cell_type} (Reason: {reasoning})\n"
        else:
            formatted_runs += "- No valid hypotheses generated.\n"
        formatted_runs += "\n"

    prompt = (
        "You are an expert biologist summarizing experimental replicates. Below are results from multiple analysis runs "
        f"for the same cell cluster: **{cluster_name}**. The predicted cell types and reasoning are provided for each run.\n\n"
        "Your task is to consolidate these results into a single, definitive list. Perform the following steps:\n"
        "1. Identify and group synonymous cell types (e.g., 'T-helper cell' is the same as 'Helper T Cell').\n"
        "2. Analyze the rankings and reasoning across all runs to determine a final, consolidated rank for each unique cell type.\n"
        "3. Provide a new, summarized reasoning for your final ranking, referencing the evidence from across the different runs.\n\n"
        "Here are the raw results:\n"
        f"-------------------\n{formatted_runs}-------------------\n\n"
        "Based on your analysis, provide a final list of the top 3 most likely cell types in the following strict format:\n\n"
        "1. **Cell Type:** [Consolidated cell type name]\n"
        "   **Reasoning:** [Synthesized reasoning based on evidence from all runs]\n\n"
        "2. **Cell Type:** [Consolidated cell type name]\n"
        "   **Reasoning:** [Synthesized reasoning based on evidence from all runs]\n\n"
        "3. **Cell Type:** [Consolidated cell type name]\n"
        "   **Reasoning:** [Synthesized reasoning based on evidence from all runs]"
    )
    return prompt

def _parse_summary_response(response: str) -> List[Dict[str, Any]]:
    """Parses the structured LLM response to extract final cell type hypotheses."""
    # This reuses the same robust parsing logic
    import re
    parsed_results = []
    pattern = re.compile(
        r"^\s*(\d+)\.\s*\*\*Cell Type:\*\*\s*(.*?)\s*\*\*Reasoning:\*\*\s*(.*?)(?=\n\s*\d+\.|\Z)",
        re.DOTALL | re.MULTILINE
    )
    matches = pattern.findall(response)
    for match in matches:
        parsed_results.append({
            "rank": int(match[0]),
            "cell_type": match[1].strip(),
            "reasoning": match[2].strip()
        })
    return parsed_results

def _consolidate_one_cluster(cluster_name: str, runs: List[Dict], **kwargs) -> Dict[str, Any]:
    """
    Helper function to run consolidation for a single cluster.
    This is designed to be called in parallel.
    """
    prompt = _create_consolidation_prompt(cluster_name, runs)
    
    summary_response = call_llm(
        prompt=prompt,
        provider="openrouter",
        model="google/gemini-2.5-flash",
        temperature=0.1
    )
    return {
        "parsed_response": _parse_summary_response(summary_response),
        "raw_response": summary_response
    }

def summarize_runs(manifest_csv_path: str, output_basename: str, **kwargs):
    """
    Summarizes hypothesis generation results from multiple JSON files into a final,
    parallelized report.

    Args:
        manifest_csv_path: Path to a CSV file containing paths to the JSON files to be summarized.
        output_basename: Base name for the output summary files.
    """
    json_files = []
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    with open(manifest_csv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader) # Skip header
        for row in reader:
            if row:
                json_files.append(row[0])

    if not json_files:
        print("No JSON files found in the manifest. Aborting.")
        return

    all_data = {}
    for f in json_files:
        if not os.path.exists(f):
            print(f"Warning: File not found: {f}. Skipping.")
            continue
        with open(f, 'r', encoding='utf-8') as infile:
            run_name = os.path.basename(f).replace('.json', '')
            data = json.load(infile)
            for cluster_name, results in data.items():
                if cluster_name not in all_data:
                    all_data[cluster_name] = []
                all_data[cluster_name].append(results)

    final_summary = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_cluster = {}
        for cluster_name, runs in all_data.items():
            print(f"--- Submitting consolidation task for cluster: {cluster_name} ---")
            future = executor.submit(_consolidate_one_cluster, cluster_name, runs, **kwargs)
            future_to_cluster[future] = cluster_name
        
        for future in concurrent.futures.as_completed(future_to_cluster):
            cluster_name = future_to_cluster[future]
            try:
                result = future.result()
                final_summary[cluster_name] = result
                print(f"--- Successfully consolidated cluster: {cluster_name} ---")
            except Exception as exc:
                print(f'--- Cluster {cluster_name} generated an exception during consolidation: {exc} ---')
                final_summary[cluster_name] = {"error": str(exc)}

    # Save the final summary data
    final_json_path = os.path.join(script_dir, f"{output_basename}_summary.json")
    with open(final_json_path, 'w', encoding='utf-8') as f:
        json.dump(final_summary, f, indent=4)
    print(f"\nFinal consolidated data saved to {final_json_path}")

    # Generate the final HTML report
    final_html_path = os.path.join(script_dir, f"{output_basename}_summary.html")
    create_html_report(final_json_path, final_html_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Consolidate multiple CASSIA hypothesis generation runs into a single summary report."
    )
    parser.add_argument(
        "manifest_file",
        help="Path to the manifest CSV file that lists the JSON files to process."
    )
    parser.add_argument(
        "-o", "--output",
        default=f"{datetime.now().strftime('%Y%m%d-%H%M%S')}_consolidation",
        help="Base name for the output summary files (e.g., 'final_report')."
    )
    args = parser.parse_args()

    summarize_runs(args.manifest_file, args.output) 