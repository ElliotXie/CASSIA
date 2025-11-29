import os
import glob
from datetime import datetime

# Import the core functions from our other scripts
from hypothesis_generation import run_multi_analysis
from summarize_hypothesis_runs import summarize_runs

def run_full_pipeline_test():
    """
    Runs the full, end-to-end hypothesis generation and summarization pipeline
    by directly calling the core functions.
    """
    print("--- Starting Full Pipeline Test ---")

    # --- Test Parameters ---
    # This is now the single source of truth for the test run.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    marker_file_path = r"C:\\Users\\ellio\\OneDrive - UW-Madison\\Revision_cassia\\3_csscore\\bone_top_50_genes (1).csv"
    species = "human"
    tissue = "bone marrow"
    num_runs = 3
    generator_model_name = "openai/gpt-4o-2024-11-20"
    formatter_model_name = "google/gemini-2.5-flash"
    # -----------------------

    # --- Step 1: Run the multi-analysis ---
    print("\n[Step 1/2] Executing multi-run analysis...")
    run_multi_analysis(
        marker_file_path=marker_file_path,
        species=species,
        tissue=tissue,
        num_runs=num_runs,
        generator_model_name=generator_model_name,
        formatter_model_name=formatter_model_name
    )
    print("[Step 1/2] Multi-run analysis complete.")

    # --- Step 2: Find the manifest and run the consolidation ---
    print("\n[Step 2/2] Searching for manifest and starting consolidation...")
    
    # Find the most recently created manifest file in the CWD
    # The CWD should be the script's directory.
    manifest_files = glob.glob(os.path.join(script_dir, "run_manifest_*.csv"))
    if not manifest_files:
        print("Error: Could not find any manifest files to process.")
        print("\n--- Full pipeline test failed. ---")
        return

    latest_manifest = max(manifest_files, key=os.path.getctime)
    print(f"Found manifest file: {os.path.basename(latest_manifest)}")

    consolidation_basename = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}_{species.replace(' ','_')}_{tissue.replace(' ','_')}_consolidated"
    
    summarize_runs(
        manifest_csv_path=latest_manifest,
        output_basename=consolidation_basename,
        generator_model=generator_model_name,
        formatter_model=formatter_model_name
    )
    print("[Step 2/2] Consolidation complete.")

    print("\n--- Full Pipeline Test Finished Successfully ---")

if __name__ == "__main__":
    # To run this test, ensure you are in the CASSIA/CASSIA_python/CASSIA directory
    # or that your PYTHONPATH is set up correctly.
    run_full_pipeline_test() 