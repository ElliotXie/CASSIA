"""
CASSIA Test Suite - Test Utilities
==================================
Common utilities for running tests.
"""

import os
import sys
import yaml
from pathlib import Path
from dotenv import load_dotenv


def get_test_root() -> Path:
    """Get the root directory of the test suite."""
    return Path(__file__).parent.parent.parent


def get_cassia_python_path() -> Path:
    """Get the path to CASSIA Python source."""
    return get_test_root().parent / "CASSIA_python" / "CASSIA"


def setup_cassia_imports():
    """Add CASSIA Python to sys.path for imports."""
    cassia_path = str(get_cassia_python_path())
    if cassia_path not in sys.path:
        sys.path.insert(0, cassia_path)


def load_config() -> dict:
    """
    Load test configuration from YAML file.

    Returns:
        dict: Configuration dictionary
    """
    config_path = get_test_root() / "config" / "test_config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def setup_api_keys():
    """
    Load API keys from config directory.

    Tries in order:
    1. set_api_keys.py (Python file with keys)
    2. api_keys.env (environment file)

    Returns:
        bool: True if keys were loaded, False otherwise
    """
    config_path = get_test_root() / "config"

    # Try loading from set_api_keys.py first
    set_api_keys_path = config_path / "set_api_keys.py"
    if set_api_keys_path.exists():
        # Add config directory to path temporarily
        config_str = str(config_path)
        if config_str not in sys.path:
            sys.path.insert(0, config_str)
        try:
            import set_api_keys
            return True
        except ImportError:
            pass

    # Fallback to api_keys.env
    env_path = config_path / "api_keys.env"
    if env_path.exists():
        load_dotenv(env_path)
        return True

    return False


def get_llm_config() -> dict:
    """
    Get LLM configuration from test config.

    Returns:
        dict: LLM configuration with provider, model, temperature, etc.
    """
    config = load_config()
    return config.get('llm', {})


def validate_annotation_result(result: dict) -> dict:
    """
    Validate that an annotation result has the expected structure.

    Args:
        result: The annotation result dictionary

    Returns:
        dict: Validation results with 'valid' bool and 'errors' list
    """
    errors = []
    required_fields = ['main_cell_type']
    optional_fields = ['sub_cell_types', 'possible_mixed_cell_types']

    for field in required_fields:
        if field not in result:
            errors.append(f"Missing required field: {field}")
        elif result[field] is None or result[field] == '':
            errors.append(f"Empty required field: {field}")

    return {
        'valid': len(errors) == 0,
        'errors': errors,
        'has_main_cell_type': 'main_cell_type' in result and result['main_cell_type'],
        'has_sub_cell_types': 'sub_cell_types' in result and result['sub_cell_types'],
    }


def print_test_header(test_name: str):
    """Print a formatted test header."""
    print("\n" + "=" * 60)
    print(f"  CASSIA TEST: {test_name}")
    print("=" * 60)


def print_test_result(success: bool, message: str = ""):
    """Print test result status."""
    status = "PASSED" if success else "FAILED"
    symbol = "[OK]" if success else "[X]"
    print(f"\n{symbol} Test {status}")
    if message:
        print(f"    {message}")


def print_config_summary(config: dict):
    """Print a summary of the test configuration."""
    llm = config.get('llm', {})
    data = config.get('data', {})
    print(f"\nConfiguration:")
    print(f"  Provider: {llm.get('provider', 'N/A')}")
    print(f"  Model: {llm.get('model', 'N/A')}")
    print(f"  Tissue: {data.get('tissue', 'N/A')}")
    print(f"  Species: {data.get('species', 'N/A')}")
