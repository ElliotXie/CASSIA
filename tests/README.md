# CASSIA Test Suite

This directory contains tests for the CASSIA (Cell type Annotation by Single-cell data analysis Improved by Artificial intelligence) package.

## Test Organization

- `test_custom_providers.py`: Tests for the custom provider functionality
- `test_marker_handling.py`: Tests for marker data parsing and manipulation
- `test_data_handling.py`: Tests for CSV and data handling functions
- `test_integration.py`: Integration tests with mock API responses
- `test_tools_function.py`: Tests for the main functions in tools_function.py, including runCASSIA, runCASSIA_batch, split_markers, and API provider handling

## Test Data

The `data/` directory contains sample data files for testing:

- `test_markers.csv`: Sample Seurat format marker data with 5 clusters
- `test_results.json`: Sample annotated results in JSON format

## Running Tests

### Run all tests

To run all tests at once:

```bash
python tests/run_tests.py
```

### Run individual test files

To run a specific test file:

```bash
python -m unittest tests/test_custom_providers.py
```

### Run specific test cases

To run a specific test case:

```bash
python -m unittest tests.test_custom_providers.TestCustomProviders.test_set_custom_provider
```

## Adding New Tests

When adding new tests:

1. Create a new test file with the prefix `test_` (e.g., `test_new_feature.py`)
2. Import required modules and the CASSIA functions you want to test
3. Create test classes that inherit from `unittest.TestCase`
4. Write test methods with names starting with `test_`
5. Use assertions to verify expected behavior

For testing functions that call external APIs, use the mock classes provided in `test_integration.py` to simulate API responses.

## Test Dependencies

- `unittest`: Python's built-in testing framework
- `pandas`: For data manipulation tests
- `unittest.mock`: For mocking external API calls

## Mocking Guidelines

When writing tests that involve API calls or other external dependencies:

1. Use `@patch` to mock the functions or classes that make API calls
2. Create `MagicMock` objects to simulate responses
3. Use `side_effect` for mock functions that need different responses based on input
4. For batch functions, create temporary files rather than relying on actual API calls 