"""
Tests for CSV and data handling functions in CASSIA.
"""

import unittest
import os
import tempfile
import pandas as pd
import csv
from CASSIA import write_csv, standardize_cell_types

class TestCSVHandling(unittest.TestCase):
    """Test CSV and data handling functions."""
    
    def setUp(self):
        """Set up test environment before each test."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()
        
        # Sample data for tests
        self.headers = ["Cluster", "Main Cell Type", "Sub Cell Type", "Score"]
        self.row_data = [
            {"Cluster": "Cluster 1", "Main Cell Type": "T cell", "Sub Cell Type": "CD4+ T cell", "Score": 98},
            {"Cluster": "Cluster 2", "Main Cell Type": "B cell", "Sub Cell Type": "Naive B cell", "Score": 95},
            {"Cluster": "Cluster 3", "Main Cell Type": "Macrophage", "Sub Cell Type": "M1 Macrophage", "Score": 92}
        ]
        
        # Path for test CSV file
        self.test_csv_path = os.path.join(self.temp_dir, "test_output.csv")
    
    def tearDown(self):
        """Clean up after tests."""
        # Remove the temporary file if it exists
        if os.path.exists(self.test_csv_path):
            os.remove(self.test_csv_path)
        
        # Remove the temporary directory
        if os.path.exists(self.temp_dir):
            os.rmdir(self.temp_dir)
    
    def test_write_csv(self):
        """Test writing data to a CSV file."""
        # Use the write_csv function
        write_csv(self.test_csv_path, self.headers, self.row_data)
        
        # Verify the file was created
        self.assertTrue(os.path.exists(self.test_csv_path), "CSV file was not created")
        
        # Read the file back and verify its contents
        with open(self.test_csv_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            # Check the header row
            header_row = next(reader)
            self.assertEqual(header_row, self.headers, "CSV headers don't match")
            
            # Check the data rows
            data_rows = list(reader)
            self.assertEqual(len(data_rows), len(self.row_data), "Wrong number of data rows")
            
            # Convert rows back to dictionaries and compare
            for i, row in enumerate(data_rows):
                expected_values = list(str(x) for x in self.row_data[i].values())
                self.assertEqual(row, expected_values, f"Data row {i} doesn't match")
    
    def test_write_csv_empty(self):
        """Test writing empty data to a CSV file."""
        # Write an empty dataset
        write_csv(self.test_csv_path, self.headers, [])
        
        # Verify the file was created
        self.assertTrue(os.path.exists(self.test_csv_path), "CSV file was not created")
        
        # Read the file back and verify it has just the header
        with open(self.test_csv_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)
            self.assertEqual(len(rows), 1, "Empty CSV should only have the header row")
            self.assertEqual(rows[0], self.headers, "CSV header doesn't match")


class TestCellTypeStandardization(unittest.TestCase):
    """Test cell type standardization functions."""
    
    def test_standardize_cell_types(self):
        """Test standardizing cell type names."""
        test_cases = [
            # Regular cell types
            ("T cell", "T cell"),
            ("T cells", "T cell"),
            ("T-cell", "T cell"),
            ("T-cells", "T cell"),
            # Cell types with additional terms
            ("CD4+ T cell", "CD4+ T cell"),
            ("CD4+ T cells", "CD4+ T cell"),
            ("CD4+ T-cell", "CD4+ T cell"),
            ("CD4+ T-cells", "CD4+ T cell"),
            # Cell types with whitespace
            ("  T cell  ", "T cell"),
            (" CD4+ T cells ", "CD4+ T cell"),
            # Mixed case
            ("t cell", "T cell"),
            ("cd4+ t cells", "CD4+ T cell"),
            # Complex cases
            ("Naive CD4+ T Helper cells", "Naive CD4+ T Helper cell"),
            ("Memory B-cells", "Memory B cell"),
            ("CD8+ Cytotoxic T-lymphocytes", "CD8+ Cytotoxic T lymphocyte"),
            # Edge cases
            ("", ""),
            (None, None)
        ]
        
        for input_name, expected_output in test_cases:
            self.assertEqual(
                standardize_cell_types(input_name), 
                expected_output,
                f"Failed for input: '{input_name}'"
            )
    
    def test_standardize_cell_types_preserves_capitalization(self):
        """Test that standardization preserves existing capitalization."""
        # Cases where capitalization should be preserved
        test_cases = [
            ("NK cell", "NK cell"),  # Acronyms stay capitalized
            ("CD4+ T cell", "CD4+ T cell"),  # CD markers stay capitalized
            ("ILC2", "ILC2"),  # Acronyms without "cell" stay the same
            ("ILC2 cell", "ILC2 cell"),  # Acronyms with "cell" stay capitalized
        ]
        
        for input_name, expected_output in test_cases:
            self.assertEqual(
                standardize_cell_types(input_name), 
                expected_output,
                f"Failed to preserve capitalization for: '{input_name}'"
            )


if __name__ == "__main__":
    unittest.main() 