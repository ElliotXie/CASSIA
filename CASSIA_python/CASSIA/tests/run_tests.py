"""
Test runner for CASSIA tests.
"""

import unittest
import os
import sys

def run_all_tests():
    """Run all tests in the tests directory."""
    # Get the directory where this script is located
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Add parent directory to path so we can import CASSIA
    parent_dir = os.path.dirname(tests_dir)
    sys.path.insert(0, parent_dir)
    
    # Discover and load all tests
    loader = unittest.TestLoader()
    suite = loader.discover(tests_dir, pattern="test_*.py")
    
    # Run the tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return success if all tests passed, otherwise fail
    return result.wasSuccessful()

if __name__ == "__main__":
    # Run all tests and exit with appropriate status code
    success = run_all_tests()
    sys.exit(0 if success else 1) 