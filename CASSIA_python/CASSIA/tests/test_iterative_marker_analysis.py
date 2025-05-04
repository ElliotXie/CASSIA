"""
Unit tests for the iterative_marker_analysis function in tools_function.py
"""

import unittest
import os
import pandas as pd
import json
from unittest.mock import patch, MagicMock, call

# Import the function to test
from CASSIA import iterative_marker_analysis, get_marker_info

class TestIterativeMarkerAnalysis(unittest.TestCase):
    """Tests for the iterative_marker_analysis function."""
    
    def setUp(self):
        """Set up test environment before each test."""
        # Clear any environment variables that might affect tests
        for var in ['OPENAI_API_KEY', 'ANTHROPIC_API_KEY', 'OPENROUTER_API_KEY']:
            if var in os.environ:
                del os.environ[var]
                
        # Set mock API keys for testing
        os.environ["OPENAI_API_KEY"] = "test-openai-api-key"
        os.environ["ANTHROPIC_API_KEY"] = "test-anthropic-api-key"
        os.environ["OPENROUTER_API_KEY"] = "test-openrouter-api-key"
        
        # Sample marker data for testing
        self.marker_data = pd.DataFrame({
            'gene': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'GZMB', 'CD19', 'MS4A1'],
            'p_val': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
            'avg_log2FC': [3.5, 3.2, 2.5, 2.0, -5.0, -3.0, -4.0, -4.5],
            'pct.1': [0.92, 0.90, 0.80, 0.75, 0.05, 0.03, 0.02, 0.01],
            'pct.2': [0.05, 0.06, 0.08, 0.10, 0.78, 0.65, 0.88, 0.82]
        }).set_index('gene')
        
        # Test input parameters
        self.major_cluster_info = "Human PBMC"
        self.comma_separated_genes = "CD3D, CD3E, CD4, IL7R"
        self.annotation_history = "Previous analysis indicates T cell markers."
        
    def tearDown(self):
        """Clean up after tests."""
        # Remove test environment variables
        for var in ['OPENAI_API_KEY', 'ANTHROPIC_API_KEY', 'OPENROUTER_API_KEY']:
            if var in os.environ:
                del os.environ[var]
                
    @patch('CASSIA.tools_function.OpenAI')
    def test_iterative_marker_analysis_openai(self, mock_openai):
        """Test iterative_marker_analysis with OpenAI provider."""
        # Set up mock responses
        mock_client = MagicMock()
        mock_openai.return_value = mock_client
        
        # First response asking to check more genes
        first_response = MagicMock()
        first_response.choices[0].message.content = """
        The T cell markers look promising. Let's check for CD8+ T cell markers to verify if this is a CD4+ or CD8+ population.
        
        <check_genes>
        CD8A, CD8B, GZMB, PRF1
        </check_genes>
        
        <reasoning>
        I need to check CD8 markers to distinguish between CD4+ and CD8+ T cells.
        </reasoning>
        """
        
        # Second response with final annotation
        second_response = MagicMock()
        second_response.choices[0].message.content = """
        Based on the gene expression patterns, this cluster shows high expression of CD3D, CD3E, CD4, and IL7R, while CD8A and GZMB are strongly downregulated.
        
        This pattern is consistent with CD4+ T cells, not CD8+ T cells.
        
        FINAL ANNOTATION COMPLETED
        
        This cluster represents CD4+ T cells, specifically likely naive or memory CD4+ T cells based on the IL7R expression.
        """
        
        # Configure the mock to return these responses in sequence
        mock_client.chat.completions.create.side_effect = [first_response, second_response]
        
        # Execute the function
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=2,
            model="gpt-4o",
            provider="openai"
        )
        
        # Verify function behavior
        self.assertIn("FINAL ANNOTATION COMPLETED", result)
        self.assertIn("CD4+ T cells", result)
        
        # Verify correct number of API calls (2 iterations)
        self.assertEqual(mock_client.chat.completions.create.call_count, 2)
        
        # Verify messages were built correctly (initial prompt, assistant, retrieved markers, final response)
        self.assertEqual(len(messages), 4)
        self.assertEqual(messages[0]["role"], "user")
        self.assertEqual(messages[1]["role"], "assistant")
        self.assertEqual(messages[2]["role"], "user")
        self.assertIn("CD8A", messages[2]["content"])  # The marker info should contain CD8A
    
    @patch('CASSIA.tools_function.anthropic')
    def test_iterative_marker_analysis_anthropic(self, mock_anthropic):
        """Test iterative_marker_analysis with Anthropic provider."""
        # Set up mock client and response for Anthropic
        mock_client = MagicMock()
        mock_anthropic.Anthropic.return_value = mock_client
        
        # First response
        first_response = MagicMock()
        first_response.content = [MagicMock(text="""
        The T cell markers are highly expressed. Let's check B cell markers to rule out any B cell contamination.
        
        <check_genes>
        CD19, MS4A1, CD79A
        </check_genes>
        
        <reasoning>
        I need to verify there's no B cell signature in this cluster.
        </reasoning>
        """)]
        
        # Second response with final annotation
        second_response = MagicMock()
        second_response.content = [MagicMock(text="""
        Based on the gene expression patterns, this cluster shows high expression of T cell markers (CD3D, CD3E, CD4, IL7R) 
        and downregulation of B cell markers (CD19, MS4A1).
        
        FINAL ANNOTATION COMPLETED
        
        This cluster represents CD4+ T cells, likely naive or central memory CD4+ T cells based on IL7R expression.
        """)]
        
        # Configure mock to return these responses
        mock_client.messages.create.side_effect = [first_response, second_response]
        
        # Execute the function
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=2,
            model="claude-3-5-sonnet-20241022",
            provider="anthropic"
        )
        
        # Verify function behavior
        self.assertIn("FINAL ANNOTATION COMPLETED", result)
        self.assertIn("CD4+ T cells", result)
        
        # Verify correct number of API calls
        self.assertEqual(mock_client.messages.create.call_count, 2)
        
        # Verify messages were built correctly
        self.assertEqual(len(messages), 4)
        self.assertEqual(messages[0]["role"], "user")
        self.assertEqual(messages[1]["role"], "assistant")
        self.assertEqual(messages[2]["role"], "user")
        self.assertIn("MS4A1", messages[2]["content"])  # The marker info should contain MS4A1
    
    @patch('CASSIA.tools_function.requests')
    def test_iterative_marker_analysis_openrouter(self, mock_requests):
        """Test iterative_marker_analysis with OpenRouter provider."""
        # Set up mock response for requests
        mock_response1 = MagicMock()
        mock_response1.status_code = 200
        mock_response1.json.return_value = {
            "choices": [
                {
                    "message": {
                        "content": """
                        The T cell markers look promising. Let's check for more specific T cell subset markers.
                        
                        <check_genes>
                        CCR7, LEF1, TCF7, FOXP3, CTLA4
                        </check_genes>
                        
                        <reasoning>
                        I need to check markers for naive, memory, and regulatory T cells.
                        </reasoning>
                        """
                    }
                }
            ]
        }
        
        mock_response2 = MagicMock()
        mock_response2.status_code = 200
        mock_response2.json.return_value = {
            "choices": [
                {
                    "message": {
                        "content": """
                        Based on the gene expression patterns, this is clearly a CD4+ T cell population.
                        
                        FINAL ANNOTATION COMPLETED
                        
                        This cluster represents CD4+ T cells, most likely a mixture of naive and central memory CD4+ T cells.
                        """
                    }
                }
            ]
        }
        
        # Configure mock to return these responses
        mock_requests.post.side_effect = [mock_response1, mock_response2]
        
        # Execute the function
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=2,
            model="anthropic/claude-3.5-sonnet",
            provider="openrouter"
        )
        
        # Verify function behavior
        self.assertIn("FINAL ANNOTATION COMPLETED", result)
        self.assertIn("CD4+ T cells", result)
        
        # Verify correct number of API calls
        self.assertEqual(mock_requests.post.call_count, 2)
        
        # Verify the URL used for OpenRouter
        for call_args in mock_requests.post.call_args_list:
            self.assertEqual(call_args[1]["url"], "https://openrouter.ai/api/v1/chat/completions")
        
        # Verify messages were built correctly
        self.assertEqual(len(messages), 4)
    
    @patch('CASSIA.tools_function.OpenAI')
    def test_iterative_marker_analysis_with_additional_task(self, mock_openai):
        """Test iterative_marker_analysis with an additional task parameter."""
        # Set up mock responses
        mock_client = MagicMock()
        mock_openai.return_value = mock_client
        
        # First response asking to check more genes
        first_response = MagicMock()
        first_response.choices[0].message.content = """
        Let's check some cancer-related markers to assess if this is a malignant T cell population.
        
        <check_genes>
        MYC, MKI67, PCNA, TOP2A
        </check_genes>
        
        <reasoning>
        These genes are associated with cell proliferation and cancer.
        </reasoning>
        """
        
        # Second response with final conclusion
        second_response = MagicMock()
        second_response.choices[0].message.content = """
        Based on the expression patterns, there's no clear evidence of malignancy.
        
        FINAL ANALYSIS COMPLETED
        
        This appears to be a normal, non-malignant CD4+ T cell population without cancer signatures.
        """
        
        # Configure the mock to return these responses
        mock_client.chat.completions.create.side_effect = [first_response, second_response]
        
        # Execute the function with additional task
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=2,
            model="gpt-4o",
            provider="openai",
            additional_task="check if this is a cancer cluster"
        )
        
        # Verify function behavior
        self.assertIn("FINAL ANALYSIS COMPLETED", result)
        self.assertIn("non-malignant", result)
        
        # Verify completion marker was correctly set to "FINAL ANALYSIS COMPLETED"
        self.assertEqual(mock_client.chat.completions.create.call_count, 2)
    
    @patch('CASSIA.tools_function.OpenAI')
    def test_iterative_marker_analysis_max_iterations(self, mock_openai):
        """Test that iterative_marker_analysis returns after max iterations even without completion marker."""
        # Set up mock responses that never include the completion marker
        mock_client = MagicMock()
        mock_openai.return_value = mock_client
        
        # First response
        first_response = MagicMock()
        first_response.choices[0].message.content = """
        Let's check some T cell subset markers.
        
        <check_genes>
        CCR7, SELL, FOXP3, RORC
        </check_genes>
        
        <reasoning>
        These will help identify specific T cell subsets.
        </reasoning>
        """
        
        # Second response (still no completion marker)
        second_response = MagicMock()
        second_response.choices[0].message.content = """
        Let's check another set of markers.
        
        <check_genes>
        GZMA, GZMK, PRF1, GNLY
        </check_genes>
        
        <reasoning>
        These will help identify effector functions.
        </reasoning>
        """
        
        # Final response when max iterations reached
        final_response = MagicMock()
        final_response.choices[0].message.content = """
        Based on the gene expression patterns, this is a CD4+ T cell population.
        """
        
        # Configure the mock to return these responses
        mock_client.chat.completions.create.side_effect = [first_response, second_response, final_response]
        
        # Execute the function with num_iterations=2
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=2,  # Should exit after 2 iterations
            model="gpt-4o",
            provider="openai"
        )
        
        # Verify function behavior - should have made a final call to get a conclusion
        self.assertEqual(mock_client.chat.completions.create.call_count, 3)
        
        # Result should be the final response
        self.assertEqual(result, final_response.choices[0].message.content)
        
        # Messages should include all exchanges
        self.assertEqual(len(messages), 5)  # Initial user, 2 iterations with assistant and user, final assistant
    
    def test_get_marker_info(self):
        """Test the get_marker_info helper function that's used by iterative_marker_analysis."""
        # Create a test DataFrame
        marker_df = pd.DataFrame({
            'p_val': [0.001, 0.005, 0.01],
            'avg_log2FC': [2.5, 1.8, -2.0],
            'pct.1': [0.8, 0.6, 0.1],
            'pct.2': [0.2, 0.3, 0.7]
        }, index=['CD3D', 'CD4', 'CD19'])
        
        # Test retrieving marker info
        gene_list = ['CD3D', 'CD19', 'CD8A']  # CD8A not in the DataFrame
        result = get_marker_info(gene_list, marker_df)
        
        # Verify result contains the expected genes
        self.assertIn('CD3D', result)
        self.assertIn('CD19', result)
        self.assertIn('CD8A', result)
        
        # Verify avg_log2FC values were formatted correctly
        self.assertIn('2.50e+00', result)  # CD3D
        self.assertIn('-2.00e+00', result)  # CD19
        
        # Verify 'NA' for missing gene
        self.assertIn('NA', result)  # CD8A should have NA values

if __name__ == '__main__':
    unittest.main() 
    
# Add a new test class for real API calls
class TestIterativeMarkerAnalysisReal(unittest.TestCase):
    """Tests for the iterative_marker_analysis function using real API calls.
    
    These tests use actual API keys from environment variables and make real API calls.
    Skip these tests if you don't have the required API keys set.
    """
    
    def setUp(self):
        """Set up test environment for real API calls."""
        # Check if API keys are set
        self.openai_key_present = "OPENAI_API_KEY" in os.environ and os.environ["OPENAI_API_KEY"]
        self.anthropic_key_present = "ANTHROPIC_API_KEY" in os.environ and os.environ["ANTHROPIC_API_KEY"]
        self.openrouter_key_present = "OPENROUTER_API_KEY" in os.environ and os.environ["OPENROUTER_API_KEY"]
        
        # Sample marker data for testing
        self.marker_data = pd.DataFrame({
            'gene': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'GZMB', 'CD19', 'MS4A1'],
            'p_val': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
            'avg_log2FC': [3.5, 3.2, 2.5, 2.0, -5.0, -3.0, -4.0, -4.5],
            'pct.1': [0.92, 0.90, 0.80, 0.75, 0.05, 0.03, 0.02, 0.01],
            'pct.2': [0.05, 0.06, 0.08, 0.10, 0.78, 0.65, 0.88, 0.82]
        }).set_index('gene')
        
        # Test input parameters
        self.major_cluster_info = "Human PBMC"
        self.comma_separated_genes = "CD3D, CD3E, CD4, IL7R"
        self.annotation_history = "Previous analysis indicates T cell markers."
        
    def test_iterative_marker_analysis_openai_real(self):
        """Test iterative_marker_analysis with real OpenAI API call."""
        if not self.openai_key_present:
            self.skipTest("OPENAI_API_KEY environment variable not set")
        
        # Execute the function with num_iterations=1 to minimize API usage
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=1,  # Just one iteration to minimize API usage
            model="gpt-3.5-turbo",  # Using a cheaper model for testing
            provider="openai"
        )
        
        print("\nOpenAI API Result:")
        print(result)
        
        # Basic verification - the result should be a non-empty string
        self.assertIsInstance(result, str)
        self.assertTrue(len(result) > 0)
        
        # Verify messages were built correctly
        self.assertGreaterEqual(len(messages), 2)
        self.assertEqual(messages[0]["role"], "user")
    
    def test_iterative_marker_analysis_anthropic_real(self):
        """Test iterative_marker_analysis with real Anthropic API call."""
        if not self.anthropic_key_present:
            self.skipTest("ANTHROPIC_API_KEY environment variable not set")
        
        # Execute the function with num_iterations=1 to minimize API usage
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=1,  # Just one iteration to minimize API usage
            model="claude-3-haiku-20240307",  # Using a cheaper model for testing
            provider="anthropic"
        )
        
        print("\nAnthropic API Result:")
        print(result)
        
        # Basic verification - the result should be a non-empty string
        self.assertIsInstance(result, str)
        self.assertTrue(len(result) > 0)
        
        # Verify messages were built correctly
        self.assertGreaterEqual(len(messages), 2)
        self.assertEqual(messages[0]["role"], "user")
    
    def test_iterative_marker_analysis_openrouter_real(self):
        """Test iterative_marker_analysis with real OpenRouter API call."""
        if not self.openrouter_key_present:
            self.skipTest("OPENROUTER_API_KEY environment variable not set")
        
        # Execute the function with num_iterations=1 to minimize API usage
        result, messages = iterative_marker_analysis(
            major_cluster_info=self.major_cluster_info,
            marker=self.marker_data,
            comma_separated_genes=self.comma_separated_genes,
            annotation_history=self.annotation_history,
            num_iterations=1,  # Just one iteration to minimize API usage
            model="openai/gpt-3.5-turbo",  # Using a cheaper model for testing
            provider="openrouter"
        )
        
        print("\nOpenRouter API Result:")
        print(result)
        
        # Basic verification - the result should be a non-empty string
        self.assertIsInstance(result, str)
        self.assertTrue(len(result) > 0)
        
        # Verify messages were built correctly
        self.assertGreaterEqual(len(messages), 2)
        self.assertEqual(messages[0]["role"], "user")

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == 'real':
        # Run only the real API tests
        suite = unittest.TestLoader().loadTestsFromTestCase(TestIterativeMarkerAnalysisReal)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        # Run all tests
        unittest.main() 