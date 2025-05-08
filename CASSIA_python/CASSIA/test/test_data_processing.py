import pytest
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock
import os

from ..tools_function import (
    split_markers,
    get_top_markers,
    loadmarker,
    list_available_markers, # Assuming this needs testing
    prepare_analysis_data, # Might need mocking for file reads
    get_marker_info,
    get_cell_type_info, # Might need mocking for API calls
    get_cell_type_info_single # Might need mocking for API calls
)

# Simulated data path (adjust as needed)
SIMULATED_MARKER_PATH = "CASSIA_python/CASSIA/test/simulated_markers.csv"

@pytest.fixture
def simulated_marker_df():
    return pd.read_csv(SIMULATED_MARKER_PATH)

def test_split_markers():
    # Test comma+space separated list
    result1 = split_markers("CD3D, CD4, IL7R, CD8A")
    assert result1 == ["CD3D", "CD4", "IL7R", "CD8A"]
    
    # Test comma-only separated list
    result2 = split_markers("CD3D,CD4,IL7R,CD8A")
    assert result2 == ["CD3D", "CD4", "IL7R", "CD8A"]
    
    # Test space-only separated list
    result3 = split_markers("CD3D CD4 IL7R CD8A")
    assert result3 == ["CD3D", "CD4", "IL7R", "CD8A"]
    
    # Test with mixed delimiters and extra spaces
    result4 = split_markers("CD3D,  CD4 IL7R,   CD8A")
    assert result4 == ["CD3D", "CD4", "IL7R", "CD8A"]
    
    # Test with empty strings
    result5 = split_markers("CD3D, , CD4, , IL7R")
    assert result5 == ["CD3D", "CD4", "IL7R"]
    
    # Test with empty input
    result6 = split_markers("")
    assert result6 == []

def test_get_top_markers(simulated_marker_df):
    # Test Seurat format
    result = get_top_markers(simulated_marker_df, n_genes=2)
    assert isinstance(result, pd.DataFrame)
    assert 'cluster' in result.columns
    assert 'markers' in result.columns
    
    # Check correct clusters are present
    clusters = set(result['cluster'])
    assert 'Cluster_0' in clusters
    assert 'Cluster_1' in clusters
    assert 'Cluster_2' in clusters
    
    # Check marker content for a specific cluster
    cluster0_markers = result[result['cluster'] == 'Cluster_0']['markers'].iloc[0]
    assert 'SNAP25' in cluster0_markers
    assert 'SYT1' in cluster0_markers
    
    # Test handling of inf values
    # Create a test dataframe with an inf value
    test_df = simulated_marker_df.copy()
    # The test data already has an 'inf' value in the CSF1R row
    
    result_inf = get_top_markers(test_df, n_genes=3)
    # Check that Cluster_2 contains CSF1R which had an inf value
    cluster2_markers = result_inf[result_inf['cluster'] == 'Cluster_2']['markers'].iloc[0]
    assert 'CSF1R' in cluster2_markers
    
    # Test empty dataframe
    empty_df = pd.DataFrame(columns=simulated_marker_df.columns)
    empty_result = get_top_markers(empty_df)
    assert empty_result.empty

def test_loadmarker():
    # This test requires mocking since it depends on package data
    with patch('importlib.resources.files', return_value=MagicMock()):
        with patch('pandas.read_csv', return_value=pd.DataFrame({'a': [1, 2, 3]})):
            # Test with default parameter
            result = loadmarker()
            assert isinstance(result, pd.DataFrame)
            
            # Test with specific parameter
            result = loadmarker("processed")
            assert isinstance(result, pd.DataFrame)

def test_list_available_markers():
    # Mock the resources.files and iterdir functions
    mock_dir = MagicMock()
    mock_file1 = MagicMock()
    mock_file1.name = "marker1.csv"
    mock_file2 = MagicMock()
    mock_file2.name = "marker2.csv"
    mock_dir.iterdir.return_value = [mock_file1, mock_file2]
    
    with patch('importlib.resources.files', return_value=mock_dir):
        markers = list_available_markers()
        assert "marker1" in markers
        assert "marker2" in markers

def test_prepare_analysis_data():
    # Create mock dataframes
    mock_result_df = pd.DataFrame({
        'cluster': ['Cluster_0', 'Cluster_1'],
        'main_cell_type': ['Neuron', 'T cell'],
        'sub_cell_type': ['Excitatory neuron', 'CD4+ T cell'],
        'confidence_score': [90, 85]
    })
    
    mock_marker_df = pd.DataFrame({
        'cluster': ['Cluster_0', 'Cluster_0', 'Cluster_1', 'Cluster_1'],
        'gene': ['SNAP25', 'SYT1', 'CD3D', 'CD4'],
        'avg_log2FC': [5.0, 4.5, 6.0, 5.5]
    })
    
    # Mock the file reading operations
    with patch('pandas.read_csv', side_effect=[mock_result_df, mock_marker_df]):
        # Call the function with a specific cluster name
        result = prepare_analysis_data("mock_result_path.csv", "mock_marker_path.csv", "Cluster_0")
        
        # Verify the results
        assert result['cluster_name'] == 'Cluster_0'
        assert result['main_cell_type'] == 'Neuron'
        assert result['sub_cell_type'] == 'Excitatory neuron'
        assert result['confidence_score'] == 90
        assert len(result['markers'].split(',')) > 0
        
        # Test with a non-existent cluster
        with pytest.raises(KeyError):
            prepare_analysis_data("mock_result_path.csv", "mock_marker_path.csv", "Cluster_99")

def test_get_marker_info(simulated_marker_df):
    # Test with a list of genes
    gene_list = ["SNAP25", "CD3D", "CD68"]
    result = get_marker_info(gene_list, simulated_marker_df)
    
    # Check that the result includes information about each gene
    assert "SNAP25" in result
    assert "CD3D" in result
    assert "CD68" in result
    
    # Test with empty gene list
    empty_result = get_marker_info([], simulated_marker_df)
    assert empty_result == "No gene information found."
    
    # Test with genes not in the marker dataframe
    not_found_result = get_marker_info(["GENE1", "GENE2"], simulated_marker_df)
    assert "No information found for the provided genes" in not_found_result

@pytest.mark.parametrize("cell_type,expected_contains", [
    ("T cell", "T cell"),
    ("NonexistentCellType", "No information found")
])
def test_get_cell_type_info(cell_type, expected_contains):
    # Mock the requests.get function
    mock_response = MagicMock()
    mock_response.json.return_value = {"_embedded": {"terms": []}} if "Nonexistent" in cell_type else {"_embedded": {"terms": [{"label": cell_type, "description": f"Description of {cell_type}"}]}}
    
    with patch('requests.get', return_value=mock_response):
        result = get_cell_type_info(cell_type)
        assert expected_contains in result

def test_get_cell_type_info_single():
    # Similar to test_get_cell_type_info but for the single version
    mock_response = MagicMock()
    mock_response.json.return_value = {"_embedded": {"terms": [{"label": "T cell", "description": "Description of T cell"}]}}
    
    with patch('requests.get', return_value=mock_response):
        result = get_cell_type_info_single("T cell")
        assert "T cell" in result 