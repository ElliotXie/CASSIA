import pytest
from unittest.mock import patch
from ..tools_function import (
    construct_prompt_from_csv_subcluster,
    annotate_subclusters,
    subcluster_agent_annotate_subcluster,
    extract_subcluster_results_with_llm_multiple_output,
    extract_subcluster_results_with_llm,
    write_results_to_csv, # Tested separately
    runCASSIA_subclusters,
    runCASSIA_n_subcluster
)

def test_construct_prompt_from_csv_subcluster():
    # TODO: Add test logic
    pass

def test_annotate_subclusters():
    # TODO: Add test logic (requires mocking LLM calls)
    pass

def test_subcluster_agent_annotate_subcluster():
     # TODO: Add test logic (requires mocking LLM calls)
    pass

def test_extract_subcluster_results_with_llm_multiple_output():
    # TODO: Add test logic (requires mocking LLM calls)
    pass

def test_extract_subcluster_results_with_llm():
    # TODO: Add test logic (requires mocking LLM calls)
    pass

def test_runCASSIA_subclusters():
    # TODO: Add test logic (requires mocking annotate and extract functions)
    pass

def test_runCASSIA_n_subcluster():
    # TODO: Add test logic (requires mocking run_single_analysis or lower level functions)
    pass 