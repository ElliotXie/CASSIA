import pytest
import pandas as pd
import json
from unittest.mock import patch, MagicMock, call
from CASSIA.tools_function import (
    runCASSIA,
    run_single_analysis, # Helper, might test via runCASSIA_n_times
    runCASSIA_n_times,
    iterative_marker_analysis_openai,
    iterative_marker_analysis,
    iterative_marker_analysis_openrouter,
    iterative_marker_analysis_openrouter_additional_task
)

# Simulated data path
SIMULATED_MARKER_PATH = "CASSIA/test/simulated_markers.csv"

@pytest.fixture
def mock_llm_response():
    return """
    {
      "main_cell_type": "T Cell",
      "sub_cell_types": "CD4+ Helper T Cell",
      "explanation": "The markers CD3D, CD4, and IL7R strongly indicate T cells, specifically CD4+ helper T cells.",
      "confidence_score": 90
    }
    """

@pytest.fixture
def mock_conversation_history():
    return [
        ("user", "Analyze these markers: CD3D, CD4, IL7R"),
        ("assistant", """
        {
          "main_cell_type": "T Cell",
          "sub_cell_types": "CD4+ Helper T Cell",
          "explanation": "The markers CD3D, CD4, and IL7R strongly indicate T cells, specifically CD4+ helper T cells.",
          "confidence_score": 90
        }
        """)
    ]

@pytest.fixture
def mock_marker_df():
    return pd.read_csv(SIMULATED_MARKER_PATH)

def test_runCASSIA(mock_llm_response, mock_conversation_history):
    parsed_result = json.loads(mock_llm_response)
    
    # Mock openai_agent to return our test response
    with patch('CASSIA.tools_function.openai_agent', return_value=mock_llm_response):
        # Mock anthropic_agent to return our test response
        with patch('CASSIA.tools_function.claude_agent', return_value=mock_llm_response):
            # Mock openrouter_agent to return our test response
            with patch('CASSIA.tools_function.openrouter_agent', return_value=mock_llm_response):
                # Mock run_cell_type_analysis to return our test data
                with patch('CASSIA.tools_function.run_cell_type_analysis', return_value=(parsed_result, mock_conversation_history)):
                    # Mock run_cell_type_analysis_claude to return our test data
                    with patch('CASSIA.tools_function.run_cell_type_analysis_claude', return_value=(parsed_result, mock_conversation_history)):
                        # Test OpenAI provider
                        result, history = runCASSIA(
                            marker_list=["CD3D", "CD4", "IL7R"],
                            tissue="blood",
                            species="human",
                            provider="openai"
                        )
                        
                        assert result["main_cell_type"].lower() == "t cell"
                        if isinstance(result["sub_cell_types"], list):
                            # Check if any item in the list matches
                            assert any(s.lower() == "cd4+ helper t cell" for s in result["sub_cell_types"])
                        else:
                            assert result["sub_cell_types"].lower() == "cd4+ helper t cell"
                        assert result["confidence_score"] == 90
                        assert isinstance(history, list)
                        
                        # Test Anthropic provider
                        result, history = runCASSIA(
                            marker_list=["CD3D", "CD4", "IL7R"],
                            tissue="blood",
                            species="human",
                            provider="anthropic"
                        )
                        
                        assert result["main_cell_type"].lower() == "t cell"
                        if isinstance(result["sub_cell_types"], list):
                            # Check if any item in the list matches
                            assert any(s.lower() == "cd4+ helper t cell" for s in result["sub_cell_types"])
                        else:
                            assert result["sub_cell_types"].lower() == "cd4+ helper t cell"
                        assert result["confidence_score"] == 90
                        assert isinstance(history, list)
                        
                        # Test OpenRouter provider (default)
                        result, history = runCASSIA(
                            marker_list=["CD3D", "CD4", "IL7R"],
                            tissue="blood",
                            species="human",
                        )
                        
                        assert result["main_cell_type"].lower() == "t cell"
                        if isinstance(result["sub_cell_types"], list):
                            # Check if any item in the list matches
                            assert any(s.lower() == "cd4+ helper t cell" for s in result["sub_cell_types"])
                        else:
                            assert result["sub_cell_types"].lower() == "cd4+ helper t cell"
                        assert result["confidence_score"] == 90
                        assert isinstance(history, list)

def test_runCASSIA_n_times(mock_llm_response, mock_conversation_history):
    parsed_result = json.loads(mock_llm_response)
    
    # Mock runCASSIA to return our test data
    with patch('CASSIA.tools_function.runCASSIA', return_value=(parsed_result, mock_conversation_history)) as mock_runCASSIA:
        # Mock ThreadPoolExecutor to avoid actual parallel execution
        with patch('concurrent.futures.ThreadPoolExecutor') as mock_executor:
            # Set up mock for ThreadPoolExecutor.map
            mock_executor.return_value.__enter__.return_value.map.side_effect = lambda fn, args: [fn(arg) for arg in args]
            
            # Call the function
            results = runCASSIA_n_times(
                n=3,
                tissue="blood",
                species="human",
                additional_info=None,
                temperature=0,
                marker_list=["CD3D", "CD4", "IL7R"],
                model="gpt-4o",
                max_workers=3,
                provider="openai"
            )
            
            # Check that runCASSIA was called 3 times
            assert mock_runCASSIA.call_count == 3
            
            # Check the result structure
            assert isinstance(results, dict)
            assert len(results) == 3
            
            # Check each result's content
            for idx, result in results.items():
                assert isinstance(result, tuple)
                assert len(result) == 2  # (parsed_json, conversation_history)
                assert result[0]["main_cell_type"].lower() == "t cell"
                assert result[0]["sub_cell_types"].lower() == "cd4+ helper t cell"
                assert result[0]["confidence_score"] == 90

@pytest.fixture
def mock_major_cluster_info():
    return {
        'cluster_name': 'Cluster_1',  # Use Cluster_1 which has CD3D and CD4 in simulated_markers.csv
        'main_cell_type': 'T Cell',
        'sub_cell_type': 'CD4+ Helper T Cell',
        'confidence_score': 90,
        'markers': 'CD3D,CD4,IL7R'
    }

"""
def test_iterative_marker_analysis_openai(mock_major_cluster_info, mock_marker_df, mock_llm_response):
    # Extract genes from Cluster_1 in the actual simulated data
    cluster1_genes = mock_marker_df[mock_marker_df['cluster'] == 'Cluster_1']['gene'].tolist()
    comma_separated_genes = ','.join(cluster1_genes)
    
    # Mock get_marker_info to return marker information
    with patch('CASSIA.tools_function.get_marker_info', return_value="CD3D: T cell marker\nCD4: Helper T cell marker\nIL7R: Cytokine receptor"):
        # Mock openai_agent to return our test response
        with patch('CASSIA.tools_function.openai_agent', return_value=mock_llm_response):
            # Call the function
            result = iterative_marker_analysis_openai(
                major_cluster_info=mock_major_cluster_info,
                marker=mock_marker_df,
                comma_separated_genes=comma_separated_genes,
                annotation_history=[],
                num_iterations=2
            )
            
            # Check result structure
            assert isinstance(result, list)
            assert len(result) == 2  # 2 iterations
            
            # Check that each iteration has the expected format
            for entry in result:
                assert isinstance(entry, tuple)
                assert len(entry) == 2  # (role, content)
                assert entry[0] in ["user", "assistant"]

def test_iterative_marker_analysis(mock_major_cluster_info, mock_marker_df, mock_llm_response):
    # Extract genes from Cluster_1 in the actual simulated data
    cluster1_genes = mock_marker_df[mock_marker_df['cluster'] == 'Cluster_1']['gene'].tolist()
    comma_separated_genes = ','.join(cluster1_genes)
    
    # Mock get_marker_info to return marker information
    with patch('CASSIA.tools_function.get_marker_info', return_value="CD3D: T cell marker\nCD4: Helper T cell marker\nIL7R: Cytokine receptor"):
        # Mock claude_agent to return our test response
        with patch('CASSIA.tools_function.claude_agent', return_value=mock_llm_response):
            # Call the function
            result = iterative_marker_analysis(
                major_cluster_info=mock_major_cluster_info,
                marker=mock_marker_df,
                comma_separated_genes=comma_separated_genes,
                annotation_history=[],
                num_iterations=2
            )
            
            # Check result structure
            assert isinstance(result, list)
            assert len(result) == 2  # 2 iterations
            
            # Check that each iteration has the expected format
            for entry in result:
                assert isinstance(entry, tuple)
                assert len(entry) == 2  # (role, content)
                assert entry[0] in ["user", "assistant"]

def test_iterative_marker_analysis_openrouter(mock_major_cluster_info, mock_marker_df, mock_llm_response):
    # Extract genes from Cluster_1 in the actual simulated data
    cluster1_genes = mock_marker_df[mock_marker_df['cluster'] == 'Cluster_1']['gene'].tolist()
    comma_separated_genes = ','.join(cluster1_genes)
    
    # Mock get_marker_info to return marker information
    with patch('CASSIA.tools_function.get_marker_info', return_value="CD3D: T cell marker\nCD4: Helper T cell marker\nIL7R: Cytokine receptor"):
        # Mock openrouter_agent to return our test response
        with patch('CASSIA.tools_function.openrouter_agent', return_value=mock_llm_response):
            # Call the function
            result = iterative_marker_analysis_openrouter(
                major_cluster_info=mock_major_cluster_info,
                marker=mock_marker_df,
                comma_separated_genes=comma_separated_genes,
                annotation_history=[],
                num_iterations=2
            )
            
            # Check result structure
            assert isinstance(result, list)
            assert len(result) == 2  # 2 iterations
            
            # Check that each iteration has the expected format
            for entry in result:
                assert isinstance(entry, tuple)
                assert len(entry) == 2  # (role, content)
                assert entry[0] in ["user", "assistant"]

def test_iterative_marker_analysis_openrouter_additional_task(mock_major_cluster_info, mock_marker_df, mock_llm_response):
    # Extract genes from Cluster_1 in the actual simulated data
    cluster1_genes = mock_marker_df[mock_marker_df['cluster'] == 'Cluster_1']['gene'].tolist()
    comma_separated_genes = ','.join(cluster1_genes)
    
    # Mock get_marker_info to return marker information
    with patch('CASSIA.tools_function.get_marker_info', return_value="CD3D: T cell marker\nCD4: Helper T cell marker\nIL7R: Cytokine receptor"):
        # Mock openrouter_agent to return our test response
        with patch('CASSIA.tools_function.openrouter_agent', return_value=mock_llm_response):
            # Call the function
            result = iterative_marker_analysis_openrouter_additional_task(
                major_cluster_info=mock_major_cluster_info,
                marker=mock_marker_df,
                comma_separated_genes=comma_separated_genes,
                annotation_history=[],
                num_iterations=2,
                additional_task="check if this is a cancer cluster"
            )
            
            # Check result structure
            assert isinstance(result, list)
            assert len(result) == 2  # 2 iterations
            
            # Check that each iteration has the expected format
            for entry in result:
                assert isinstance(entry, tuple)
                assert len(entry) == 2  # (role, content)
                assert entry[0] in ["user", "assistant"]
                
                # Check that the additional task is included in user queries
                if entry[0] == "user":
                    assert "cancer" in entry[1].lower()
""" 