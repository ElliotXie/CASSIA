"""
Integration tests for CASSIA's core cell type analysis functionality.
Using mock API responses to simulate model calls.
"""

import unittest
import json
import os
import re
from unittest.mock import patch, MagicMock
import pandas as pd
from CASSIA import runCASSIA, set_custom_provider
from io import StringIO

class MockResponse:
    """Mock response for API calls."""
    
    def __init__(self, content, status_code=200):
        self.content = content
        self.status_code = status_code

    def json(self):
        return self.content

class MockAPIClient:
    """Mock OpenAI and similar API clients."""
    
    def __init__(self, *args, **kwargs):
        self.chat = MagicMock()
        self.chat.completions = MagicMock()
        self.chat.completions.create = MagicMock(return_value=self._mock_completion())
        self.messages = MagicMock()
        self.messages.create = MagicMock(return_value=self._mock_completion())
    
    def _mock_completion(self):
        """Create a mock completion response."""
        response = MagicMock()
        message = MagicMock()
        message.content = """
I'll analyze these markers step by step.

1. Key Functional Markers:
   - CD3D, CD3E, CD2: T cell receptor complex and T cell activation
   - IL7R: Cytokine receptor important for T cell development
   - CCR7: Chemokine receptor for T cell migration
   - LTB: Lymphotoxin-beta, involved in lymphoid development

2. Key Cell Type Markers:
   - CD3D, CD3E: Pan-T cell markers
   - CD4: Helper T cell marker
   - IL7R: Naive T cells
   - CCR7: Central memory T cells

3. Database Cross-reference:
   These markers match T cell profiles in single-cell databases.

4. Most Probable General Cell Type:
   T lymphocytes

5. Top 3 Most Probable Sub-Cell Types:
   1. CD4+ Naive T cells - highest expression of CD4, IL7R, CCR7
   2. CD4+ Central Memory T cells - CD4+, CCR7+, but less IL7R
   3. CD4+ Effector T cells - CD4+ but lower CCR7

The expression pattern strongly indicates CD4+ Naive T cells.

6. Summary:
   This cluster represents CD4+ Naive T cells based on the expression of T cell markers (CD3D, CD3E), the CD4 co-receptor, and naive T cell markers IL7R and CCR7.

FINAL ANNOTATION COMPLETED
"""
        
        response.content = [message]
        response.choices = [MagicMock()]
        response.choices[0].message = MagicMock()
        response.choices[0].message.content = message.content
        return response

class MockRequestsPost:
    """Mock for requests.post."""
    
    def __init__(self, url, headers=None, json=None):
        self.url = url
        self.headers = headers
        self.json_data = json
        self.status_code = 200
    
    def json(self):
        """Return mock response JSON."""
        return {
            "choices": [
                {
                    "message": {
                        "content": """
I'll analyze these markers step by step.

1. Key Functional Markers:
   - CD3D, CD3E, CD2: T cell receptor complex and T cell activation
   - IL7R: Cytokine receptor important for T cell development
   - CCR7: Chemokine receptor for T cell migration
   - LTB: Lymphotoxin-beta, involved in lymphoid development

2. Key Cell Type Markers:
   - CD3D, CD3E: Pan-T cell markers
   - CD4: Helper T cell marker
   - IL7R: Naive T cells
   - CCR7: Central memory T cells

3. Database Cross-reference:
   These markers match T cell profiles in single-cell databases.

4. Most Probable General Cell Type:
   T lymphocytes

5. Top 3 Most Probable Sub-Cell Types:
   1. CD4+ Naive T cells - highest expression of CD4, IL7R, CCR7
   2. CD4+ Central Memory T cells - CD4+, CCR7+, but less IL7R
   3. CD4+ Effector T cells - CD4+ but lower CCR7

The expression pattern strongly indicates CD4+ Naive T cells.

6. Summary:
   This cluster represents CD4+ Naive T cells based on the expression of T cell markers (CD3D, CD3E), the CD4 co-receptor, and naive T cell markers IL7R and CCR7.

FINAL ANNOTATION COMPLETED
"""
                    }
                }
            ]
        }


class TestCellTypeAnalysis(unittest.TestCase):
    """Integration tests for CASSIA's cell type analysis."""
    
    def setUp(self):
        """Set up test environment."""
        # Set some environment variables for testing
        os.environ["OPENAI_API_KEY"] = "test-openai-api-key"
        
        # Set up a custom provider for testing
        set_custom_provider(
            api_key="test-deepseek-api-key",
            provider_name="deepseek",
            base_url="https://api.deepseek.com"
        )
        
        # Create sample marker data
        self.markers = ["CD3D", "CD3E", "CD2", "CD4", "IL7R", "CCR7", "LTB"]
    
    @patch('openai.OpenAI', MockAPIClient)
    @patch('anthropic.Anthropic', MockAPIClient)
    @patch('requests.post', MockRequestsPost)
    def test_runCASSIA_openai(self):
        """Test runCASSIA with OpenAI provider."""
        result, _ = runCASSIA(
            model="gpt-4o",
            temperature=0, 
            marker_list=self.markers, 
            tissue="blood", 
            species="human", 
            provider="openai"
        )
        
        # Check that the result contains expected fields
        self.assertIn("main_cell_type", result)
        self.assertIn("sub_cell_types", result)
        self.assertEqual(result["num_markers"], len(self.markers))
        
        # The mock response is set up for CD4+ T cells
        self.assertIn("T", result["main_cell_type"])
    
    @patch('openai.OpenAI', MockAPIClient)
    @patch('anthropic.Anthropic', MockAPIClient)
    @patch('requests.post', MockRequestsPost)
    def test_runCASSIA_custom_provider(self):
        """Test runCASSIA with a custom provider."""
        result, _ = runCASSIA(
            model="deepseek-llm",
            temperature=0, 
            marker_list=self.markers, 
            tissue="blood", 
            species="human", 
            provider="deepseek"
        )
        
        # Check that the result contains expected fields
        self.assertIn("main_cell_type", result)
        self.assertIn("sub_cell_types", result)
        self.assertEqual(result["num_markers"], len(self.markers))
        
        # The mock response is set up for CD4+ T cells
        self.assertIn("T", result["main_cell_type"])
    
    @patch('openai.OpenAI', MockAPIClient)
    @patch('anthropic.Anthropic', MockAPIClient)
    @patch('requests.post', MockRequestsPost)
    def test_runCASSIA_tissue_blind(self):
        """Test runCASSIA in tissue-blind mode."""
        result, _ = runCASSIA(
            model="gpt-4o",
            temperature=0, 
            marker_list=self.markers, 
            tissue="none",  # tissue blind
            species="human", 
            provider="openai"
        )
        
        # Check that the result contains expected fields
        self.assertIn("main_cell_type", result)
        self.assertIn("sub_cell_types", result)
        
        # In tissue-blind mode, "possible_tissues" should be present
        self.assertIn("possible_tissues", result)
    
    @patch('openai.OpenAI', MockAPIClient)
    @patch('anthropic.Anthropic', MockAPIClient)
    @patch('requests.post', MockRequestsPost)
    def test_runCASSIA_with_additional_info(self):
        """Test runCASSIA with additional context information."""
        result, _ = runCASSIA(
            model="gpt-4o",
            temperature=0, 
            marker_list=self.markers, 
            tissue="blood", 
            species="human", 
            additional_info="Sample from peripheral blood of healthy donor",
            provider="openai"
        )
        
        # Check that the result contains expected fields
        self.assertIn("main_cell_type", result)
        self.assertIn("sub_cell_types", result)
        self.assertEqual(result["num_markers"], len(self.markers))
        
        # The mock response is set up for CD4+ T cells
        self.assertIn("T", result["main_cell_type"])


class TestBatchProcessing(unittest.TestCase):
    """Test batch processing functionality with mocks."""
    
    def setUp(self):
        """Set up test environment."""
        # Create a small test CSV for batch processing
        self.test_csv = """cluster,markers
1,CD3D,CD3E,CD2,CD4,IL7R,CCR7,LTB
2,CD19,MS4A1,CD79A,CD79B,HLA-DRA,CD74
3,CD14,LYZ,CSF1R,FCGR3A,S100A8,S100A9"""
        
        # Save to a temp file
        with open("tests/test_data.csv", "w") as f:
            f.write(self.test_csv)
    
    def tearDown(self):
        """Clean up after tests."""
        try:
            os.remove("tests/test_data.csv")
        except:
            pass
        try:
            os.remove("test_output.json")
        except:
            pass
        try:
            os.remove("test_output_summary.csv")
        except:
            pass
        try:
            os.remove("test_output_full.csv")
        except:
            pass
    
    @patch('openai.OpenAI', MockAPIClient)
    @patch('anthropic.Anthropic', MockAPIClient)
    @patch('requests.post', MockRequestsPost)
    def test_batch_processing(self):
        """Test batch processing functionality with mocks."""
        # This test is skipped as it requires more complex setup
        # and would be better suited for a full integration test
        self.skipTest("Skipping batch processing test - implementation would require more complex mocking")


if __name__ == "__main__":
    unittest.main() 