import os
import pandas as pd
import json
from unittest.mock import patch
from subclustering import runCASSIA_subclusters

def test_subclustering_with_mock():
    """
    Test the subclustering functionality with a mocked API response.
    This simulates a custom API provider like DeepSeek without making actual API calls.
    """
    print("Setting up test environment...")
    # Set a dummy API key for testing
    os.environ['CUSTERMIZED_API_KEY'] = 'dummy_key'
    
    # Create a simple test dataframe
    df = pd.DataFrame({
        'cluster': ['cluster1', 'cluster2', 'cluster3', 'cluster4'],
        'markers': [
            'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3',
            'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13'
        ]
    })
    print("Test dataframe created:")
    print(df.head())
    
    # Sample structured response that mimics what DeepSeek or other custom APIs might return
    mock_response = [
        {
            'cluster': 1,
            'key_markers': 'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'explanation': 'The presence of IL7R, CD8A, and CD8B suggests a CD8+ T cell identity.',
            'most_likely_top2_cell_types': ['CD8+ memory T cells', 'Tissue-resident memory CD8+ T cells (TRM)']
        },
        {
            'cluster': 2,
            'key_markers': 'LAYN, HAVCR2 (TIM-3), TIGIT, IKZF2, KLRC2, KLRC3',
            'explanation': 'LAYN (Lag-3) and HAVCR2 (TIM-3) are markers of exhausted CD8+ T cells.',
            'most_likely_top2_cell_types': ['Exhausted CD8+ T cells', 'NK-like CD8+ T cells']
        },
        {
            'cluster': 3,
            'key_markers': 'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'explanation': 'GZMK, GZMH, PRF1, and NKG7 are markers of cytotoxic activity.',
            'most_likely_top2_cell_types': ['Effector CD8+ T cells', 'Central memory CD8+ T cells']
        },
        {
            'cluster': 4,
            'key_markers': 'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'explanation': 'WFDC2, CEACAM7, and CLDN8 are markers associated with epithelial cells.',
            'most_likely_top2_cell_types': ['Regulatory CD8+ T cells', 'Epithelial-like CD8+ T cells (rare subset)']
        }
    ]
    
    print("Mocking API call with structured response...")
    # Patch the call_llm function to return our mock response directly
    # DeepSeek and similar APIs often return structured data directly, not a JSON string
    with patch('subclustering.call_llm') as mock_call:
        # Return the structured data directly, not as a JSON string
        mock_call.return_value = mock_response
        
        # Run the subclustering function with our mock
        print("Running subclustering with mocked API...")
        runCASSIA_subclusters(
            marker=df,
            major_cluster_info='cd8 t cell',
            output_name='mock_subclustering_result',
            model='dummy-model',
            provider='https://dummy.api'
        )
    
    # Check if the output file was created
    output_file = 'mock_subclustering_result.csv'
    if os.path.exists(output_file):
        print(f"\nSuccess! Output file created: {output_file}")
        result_df = pd.read_csv(output_file)
        print("Output file contents:")
        print(result_df)
        return True
    else:
        print(f"\nError: Output file {output_file} was not created")
        return False

if __name__ == "__main__":
    test_subclustering_with_mock() 