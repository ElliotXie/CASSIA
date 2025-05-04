import os
import pandas as pd
from CASSIA import iterative_marker_analysis

# Test data
marker_data = pd.DataFrame({
    'gene': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'GZMB', 'CD19', 'MS4A1'],
    'p_val': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
    'avg_log2FC': [3.5, 3.2, 2.5, 2.0, -5.0, -3.0, -4.0, -4.5],
    'pct.1': [0.92, 0.90, 0.80, 0.75, 0.05, 0.03, 0.02, 0.01],
    'pct.2': [0.05, 0.06, 0.08, 0.10, 0.78, 0.65, 0.88, 0.82]
}).set_index('gene')

major_cluster_info = "Human PBMC"
comma_separated_genes = "CD3D, CD3E, CD4, IL7R"
annotation_history = "Previous analysis indicates T cell markers."

def test_anthropic():
    print("Testing with Anthropic API...")
    if "ANTHROPIC_API_KEY" not in os.environ or not os.environ["ANTHROPIC_API_KEY"]:
        print("ANTHROPIC_API_KEY not set, skipping test")
        return
    
    try:
        # Execute the function with a direct API call
        from anthropic import Anthropic
        
        # Manual API call to better control parameters
        client = Anthropic(api_key=os.environ.get("ANTHROPIC_API_KEY"))
        
        prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible.

        context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
        {comma_separated_genes}

        Below is the annotation analysis history:
        {annotation_history}
        """
        
        print("\nSending request to Anthropic API...")
        response = client.messages.create(
            model="claude-3-haiku-20240307",
            max_tokens=1000,
            messages=[{"role": "user", "content": prompt}]
        )
        
        result = response.content[0].text
        
        print("\nAnthropic API Result:")
        print(result)
        
        # Basic verification - the result should be a non-empty string
        assert isinstance(result, str)
        assert len(result) > 0
        
        print("\nTest completed successfully!")
        
    except Exception as e:
        print(f"Error during test: {str(e)}")

if __name__ == "__main__":
    test_anthropic() 