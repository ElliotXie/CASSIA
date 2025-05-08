import pytest
import os
import tempfile
from ..tools_function import (
    safe_get,
    write_csv # Tested elsewhere
)

def test_safe_get():
    # Test with a normal dictionary
    test_dict = {
        "level1": {
            "level2": {
                "level3": "value"
            }
        }
    }
    
    # Test successful nested access
    assert safe_get(test_dict, "level1", "level2", "level3") == "value"
    
    # Test missing key at level 1
    assert safe_get(test_dict, "missing_key") == "N/A"
    
    # Test missing key at level 2
    assert safe_get(test_dict, "level1", "missing_key") == "N/A"
    
    # Test missing key at level 3
    assert safe_get(test_dict, "level1", "level2", "missing_key") == "N/A"
    
    # Test with non-dict object
    assert safe_get("not_a_dict", "any_key") == "N/A"
    
    # Test with empty dict
    assert safe_get({}, "any_key") == "N/A"
    
    # Test with None
    assert safe_get(None, "any_key") == "N/A"

# write_csv is implicitly tested by other workflows, 
# but could have specific unit tests here if needed. 
def test_write_csv():
    # Create a temporary file for testing
    with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as temp_file:
        temp_filename = temp_file.name
    
    try:
        # Define test data
        headers = ["Name", "Age", "City"]
        row_data = [
            ["Alice", "30", "New York"],
            ["Bob", "25", "San Francisco"]
        ]
        
        # Write data to the CSV file
        write_csv(temp_filename, headers, row_data)
        
        # Read the CSV file and verify its contents
        import csv
        with open(temp_filename, 'r', newline='', encoding='utf-8') as csv_file:
            reader = csv.reader(csv_file)
            read_data = list(reader)
            
            # Check that the headers match
            assert read_data[0] == headers
            
            # Check that the data rows match
            for i, row in enumerate(row_data):
                assert read_data[i+1] == row
    
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_filename):
            os.remove(temp_filename) 