import os
import pandas as pd
from CASSIA.cell_type_comparison import compareCelltypes

# Make sure to set your OPENROUTER_API_KEY as an environment variable
# For example, in your terminal: export OPENROUTER_API_KEY='your_key_here'

def main():
    """
    Example script to test the cell type comparison function with discussion mode.
    """
    print("Starting cell type comparison with discussion mode test...")

    # --- Test Case 1: Triggering Discussion Mode ---
    # Example data that is likely to cause disagreement among models
    tissue = "Peripheral Blood"
    celltypes_to_compare = ["Naive B cell", "Plasma cell", "Plasmablast"]
    
    # A more contradictory marker set to challenge the models, including a T-cell marker (CD3)
    # and markers that change dynamically during B-cell differentiation (PAX5, CD27, etc.)
    marker_set = "CD19, CD20, CD38, SDC1, PAX5, IRF4, CD27, CD3"
    
    print(f"\n--- Running Test Case: {', '.join(celltypes_to_compare)} in {tissue} ---")
    print(f"Using marker set: {marker_set}")

    # --- Output folder setup ---
    base_results_dir = os.path.join(os.path.dirname(__file__), 'CASSIA', 'test_results')
    subfolder_name = 'budget_discussion_mode_test'
    results_dir = os.path.join(base_results_dir, subfolder_name)
    os.makedirs(results_dir, exist_ok=True)
    output_file = os.path.join(results_dir, 'results.csv')
    html_file = os.path.join(results_dir, 'results_report.html')
    
    # Run the comparison with the "budget" preset and discussion mode enabled
    comparison_results = compareCelltypes(
        tissue=tissue,
        celltypes=celltypes_to_compare,
        marker_set=marker_set,
        species="human",
        model_preset="budget",  # Use the new budget model group
        discussion_mode=True,      # Enable discussion mode
        discussion_rounds=3,       # Perform up to three rounds of discussion
        generate_html_report=True,
        output_file=output_file
    )

    # Display results
    if comparison_results and "results_df" in comparison_results:
        print("\n--- Final Results ---")
        final_df = comparison_results["results_df"]
        display_columns = ['model', 'round', 'status'] + [col for col in final_df.columns if '_score' in col]
        print(final_df[display_columns])
        print(f"\nTest finished. All rounds (initial and discussion) are saved in: {results_dir}")
        print(f"CSV: {output_file}")
        print(f"HTML: {html_file}")
    else:
        print("\nTest failed or returned no results.")

if __name__ == "__main__":
    # Check for API key before running
    if not os.environ.get('OPENROUTER_API_KEY'):
        print("ERROR: OPENROUTER_API_KEY environment variable is not set.")
        print("Please set it before running the script.")
        print("e.g., export OPENROUTER_API_KEY='your_api_key'")
    else:
        main() 