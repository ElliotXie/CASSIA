import pytest
import pandas as pd
import json
import os
from unittest.mock import patch, MagicMock, call
from CASSIA.tools_function import (
    runCASSIA_batch,
    runCASSIA_batch_n_times
)

# Simulated data path
SIMULATED_MARKER_PATH = "CASSIA/test/simulated_markers.csv"

@pytest.fixture
def mock_marker_df():
    return pd.read_csv(SIMULATED_MARKER_PATH)

@pytest.fixture
def mock_runCASSIA_output():
    return ({
        "main_cell_type": "Test Cell Type",
        "sub_cell_types": "Test Subtype",
        "explanation": "This is a test explanation",
        "confidence_score": 90
    }, [("user", "test prompt"), ("assistant", "test response")])

def test_runCASSIA_batch(mock_marker_df, mock_runCASSIA_output):
    # Mock read_csv to return our test dataframe
    with patch('pandas.read_csv', return_value=mock_marker_df):
        # Mock runCASSIA to return our test output
        with patch('CASSIA.tools_function.runCASSIA', return_value=mock_runCASSIA_output):
            # Mock os.makedirs to avoid creating directories
            with patch('os.makedirs'):
                # Mock open and json.dump to avoid writing files
                mock_open = MagicMock()
                with patch('builtins.open', mock_open):
                    with patch('json.dump'):
                        # Call the function with test parameters
                        result = runCASSIA_batch(
                            marker=SIMULATED_MARKER_PATH,
                            output_name="test_output.json",
                            tissue="lung",
                            species="human",
                            max_workers=1  # Use 1 worker to make testing simpler
                        )
                        
                        # Verify the result
                        assert isinstance(result, dict)
                        assert len(result) == 3  # Should have 3 clusters
                        
                        # Check that each cluster has the expected structure based on actual output
                        for cluster in ['Cluster_0', 'Cluster_1', 'Cluster_2']:
                            assert cluster in result
                            assert 'analysis_result' in result[cluster]
                            assert 'conversation_history' in result[cluster]
                            assert 'iterations' in result[cluster]
                            
                            # Check analysis_result content
                            analysis = result[cluster]['analysis_result']
                            assert 'main_cell_type' in analysis
                            assert 'sub_cell_types' in analysis
                            assert 'explanation' in analysis
                            assert 'confidence_score' in analysis
                        
                        # Verify that open was called to write the output file
                        mock_open.assert_called()

def test_runCASSIA_batch_n_times(mock_marker_df, mock_runCASSIA_output):
    # Mock runCASSIA_batch to return a dict of results
    mock_batch_result = {
        'Cluster_0': {
            "main_cell_type": "Neuron",
            "sub_cell_types": "Excitatory Neuron",
            "explanation": "Test explanation",
            "confidence_score": 90
        },
        'Cluster_1': {
            "main_cell_type": "T Cell",
            "sub_cell_types": "CD4+ T Cell",
            "explanation": "Test explanation",
            "confidence_score": 85
        }
    }
    
    # Set up the mocks
    with patch('CASSIA.tools_function.runCASSIA_batch', return_value=mock_batch_result) as mock_runCASSIA_batch:
        # Mock os.makedirs to avoid creating directories
        with patch('os.makedirs'):
            # Mock open and json.dump to avoid writing files
            with patch('builtins.open', MagicMock()):
                with patch('json.dump'):
                    # Call the function with test parameters
                    result = runCASSIA_batch_n_times(
                        n=3,
                        marker=SIMULATED_MARKER_PATH,
                        output_name="test_output",
                        tissue="lung",
                        species="human",
                        batch_max_workers=1,  # Use 1 worker to make testing simpler
                        max_workers=1
                    )
                    
                    # Verify runCASSIA_batch was called 3 times
                    assert mock_runCASSIA_batch.call_count == 3
                    
                    # Verify the expected calls to runCASSIA_batch
                    expected_calls = [
                        call(
                            marker=SIMULATED_MARKER_PATH,
                            output_name=f"test_output_{i}.json",
                            model="google/gemini-2.5-flash-preview",
                            temperature=0,
                            tissue="lung",
                            species="human",
                            additional_info=None,
                            celltype_column=None,
                            gene_column_name=None,
                            max_workers=1,
                            provider="openrouter",
                            max_retries=1
                        )
                        for i in range(3)
                    ]
                    mock_runCASSIA_batch.assert_has_calls(expected_calls)
                    
                    # Verify the function returns None (based on actual implementation)
                    assert result is None 