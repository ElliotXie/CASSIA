"""
Test file for testing the iterative_marker_analysis function in tools_function.py
with the development version that includes the provider parameter.
"""

import os
import pandas as pd
import io
import csv
import tempfile
from CASSIA.tools_function import iterative_marker_analysis, get_marker_info
from pathlib import Path
import inspect

# Confirm we're using the right version
print(f"Using iterative_marker_analysis from: {inspect.getfile(iterative_marker_analysis)}")
print(f"Function signature: {inspect.signature(iterative_marker_analysis)}")

def create_sample_marker_data():
    """Create a sample marker dataframe with 10 genes."""
    marker_data = pd.DataFrame({
        'gene': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'GZMB', 'CD19', 'MS4A1', 'FOXP3', 'IL2RA'],
        'p_val': [0.001, 0.001, 0.001, 0.001, 0.01, 0.01, 0.05, 0.05, 0.02, 0.02],
        'avg_log2FC': [3.5, 3.2, 2.5, 2.0, -3.0, -2.5, -4.0, -3.8, 0.5, 0.6],
        'pct.1': [0.92, 0.90, 0.85, 0.75, 0.08, 0.06, 0.05, 0.04, 0.25, 0.22],
        'pct.2': [0.05, 0.06, 0.08, 0.10, 0.65, 0.60, 0.70, 0.65, 0.10, 0.08]
    }).set_index('gene')
    
    return marker_data

def create_sample_annotation_history():
    """Create a simulated annotation history."""
    annotation_history = """
    User: Please analyze this cluster with markers CD3D, CD3E, CD4, IL7R.
    
    Assistant: Based on the markers provided, this appears to be a T cell population, particularly CD4+ T cells.
    The high expression of CD3D and CD3E indicates T cell lineage, while CD4 expression suggests CD4+ T cells specifically.
    IL7R (CD127) is often expressed on naive and memory T cells.
    
    To confirm this annotation, I would suggest checking:
    1. Additional T cell markers to confirm the lineage
    2. Markers to distinguish between naive and memory subsets
    3. Regulatory T cell markers to rule out Tregs
    
    User: Additional analysis showed low expression of CD8A and GZMB, but also low CD19 and MS4A1.
    """
    
    return annotation_history

def create_sample_full_results_csv():
    """Create a simulated full results CSV file for testing."""
    csv_content = """True Cell Type,Predicted Main Cell Type,Predicted Sub Cell Types,Possible Mixed Cell Types,Marker Number,Marker List,Iterations,Model,Provider,Tissue,Species,Additional Info,Conversation History
Cluster1,T cells,CD4+ T cells,None,10,"CD3D, CD3E, CD4, IL7R, CD8A, GZMB, CD19, MS4A1, FOXP3, IL2RA",1,gpt-4o,openai,PBMC,human,None,"User: Please analyze these markers. | Assistant: Based on the expression pattern, this appears to be a T cell population."
Cluster2,B cells,Naive B cells,None,10,"MS4A1, CD19, CD79A, CD79B, TCL1A, CD24, CD3D, CD3E, CD4, CD8A",1,gpt-4o,openai,PBMC,human,None,"User: Please analyze these markers. | Assistant: This shows clear B cell markers with MS4A1 and CD19 highly expressed."
"""
    
    # Create a temporary file to simulate the CSV
    temp_csv = tempfile.NamedTemporaryFile(delete=False, suffix='.csv', mode='w')
    temp_csv.write(csv_content)
    temp_csv.close()
    
    return temp_csv.name

def test_iterative_marker_analysis():
    """Test the iterative_marker_analysis function with sample data."""
    # Create sample data
    marker_data = create_sample_marker_data()
    annotation_history = create_sample_annotation_history()
    full_results_path = create_sample_full_results_csv()
    
    # Parameters for the function
    major_cluster_info = "Human PBMC"
    comma_separated_genes = "CD3D, CD3E, CD4, IL7R, CD8A, GZMB, CD19, MS4A1, FOXP3, IL2RA"
    
    # Check if API keys are set (choose an appropriate provider based on available keys)
    provider = None
    model = None
    
    if "OPENAI_API_KEY" in os.environ and os.environ["OPENAI_API_KEY"]:
        provider = "openai"
        model = "gpt-4.1-mini-2025-04-14" # Using a cheaper model for testing
    elif "ANTHROPIC_API_KEY" in os.environ and os.environ["ANTHROPIC_API_KEY"]:
        provider = "anthropic"
        model = "claude-3-haiku-20240307" # Using a cheaper model for testing
    elif "OPENROUTER_API_KEY" in os.environ and os.environ["OPENROUTER_API_KEY"]:
        provider = "openrouter"
        model = "deepseek/deepseek-chat-v3-0324" # Using a cheaper model for testing
    
    if provider is None:
        print("No API keys set, skipping the actual API call test")
        return
    
    print(f"Testing iterative_marker_analysis with provider: {provider}, model: {model}")
    
    try:
        # Run the function with only 1 iteration to minimize API usage
        result, messages = iterative_marker_analysis(
            major_cluster_info=major_cluster_info,
            marker=marker_data,
            comma_separated_genes=comma_separated_genes,
            annotation_history=annotation_history,
            num_iterations=1,
            model=model,
            provider=provider
        )
        
        # Print the result
        print("\nAPI Result:")
        print(result)
        
        # Basic verification
        assert isinstance(result, str)
        assert len(result) > 0
        assert len(messages) >= 2
        
        print("\nTest completed successfully!")
        
    except Exception as e:
        print(f"Error during test: {str(e)}")
    
    finally:
        # Clean up
        try:
            os.unlink(full_results_path)
        except:
            pass

def example_usage():
    """Example of how to use the iterative_marker_analysis function."""
    # 1. Create marker data
    marker_data = create_sample_marker_data()
    
    # 2. Define annotation history
    annotation_history = create_sample_annotation_history()
    
    # 3. Set up other required parameters
    major_cluster_info = "Human PBMC"
    comma_separated_genes = "CD3D, CD3E, CD4, IL7R, CD8A, GZMB, CD19, MS4A1, FOXP3, IL2RA"
    
    # 4. Choose provider and model (OpenAI example)
    provider = "openai"
    model = "gpt-4o"  # or any other model
    
    # 5. Call the function (commented out to avoid actual API calls during examples)
    """
    result, messages = iterative_marker_analysis(
        major_cluster_info=major_cluster_info,
        marker=marker_data,
        comma_separated_genes=comma_separated_genes,
        annotation_history=annotation_history,
        num_iterations=2,  # Number of iterations
        model=model,
        provider=provider
    )
    
    # 6. Process the result
    print(result)
    """
    
    print("Example data generated successfully")
    print(f"Sample marker data shape: {marker_data.shape}")
    print(f"Sample annotation history length: {len(annotation_history)}")
    
    return {
        "marker_data": marker_data,
        "annotation_history": annotation_history,
        "major_cluster_info": major_cluster_info,
        "comma_separated_genes": comma_separated_genes
    }

if __name__ == "__main__":
    # When run directly, show example usage without making API calls
    example_data = example_usage()
    
    # Uncomment below to run the actual test with API calls
    # test_iterative_marker_analysis() 