"""
Configuration loader and validator for CASSIA tests.

This module provides utilities for loading and validating test configurations.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, Optional


class TestConfig:
    """Test configuration manager."""

    def __init__(self, config_path: Path):
        """
        Initialize configuration from JSON file.

        Args:
            config_path: Path to config.json file
        """
        self.config_path = config_path
        self.config = self._load_config()
        self._validate_config()

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from JSON file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with open(self.config_path, 'r') as f:
            return json.load(f)

    def _validate_config(self):
        """Validate required fields exist in configuration."""
        required = ['test_name', 'model', 'provider']
        missing = [field for field in required if field not in self.config]

        if missing:
            raise ValueError(
                f"Missing required configuration fields: {', '.join(missing)}"
            )

    def get(self, key: str, default=None):
        """
        Get configuration value.

        Args:
            key: Configuration key
            default: Default value if key not found

        Returns:
            Configuration value or default
        """
        return self.config.get(key, default)

    def __getitem__(self, key):
        """Get configuration value using bracket notation."""
        return self.config[key]

    def __contains__(self, key):
        """Check if key exists in configuration."""
        return key in self.config

    def to_dict(self) -> Dict[str, Any]:
        """Return configuration as dictionary."""
        return self.config.copy()


def load_config(config_path: Path) -> TestConfig:
    """
    Load and validate test configuration.

    Args:
        config_path: Path to config.json file

    Returns:
        TestConfig object

    Example:
        >>> config = load_config(Path("config.json"))
        >>> model = config['model']
    """
    return TestConfig(config_path)


def get_data_path(filename: str) -> Path:
    """
    Get path to data file in CASSIA data directory.

    Args:
        filename: Name of data file (e.g., 'processed.csv')

    Returns:
        Path to data file

    Example:
        >>> data_path = get_data_path("processed.csv")
        >>> df = pd.read_csv(data_path)
    """
    # Navigate up from shared/ to CASSIA/ then to data/
    data_dir = Path(__file__).parent.parent.parent / "data"

    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    data_file = data_dir / filename

    if not data_file.exists():
        raise FileNotFoundError(f"Data file not found: {data_file}")

    return data_file


def get_results_path(test_dir: Path, filename: str) -> Path:
    """
    Get path for results file with automatic directory creation.

    Args:
        test_dir: Test directory (e.g., Path to 01_runCASSIA_batch/)
        filename: Name of results file

    Returns:
        Path to results file

    Example:
        >>> results_path = get_results_path(
        ...     Path("01_runCASSIA_batch"),
        ...     "results.csv"
        ... )
    """
    results_dir = test_dir / "results"
    results_dir.mkdir(exist_ok=True)
    return results_dir / filename


def get_reports_path(test_dir: Path, filename: str) -> Path:
    """
    Get path for report file with automatic directory creation.

    Args:
        test_dir: Test directory
        filename: Name of report file

    Returns:
        Path to report file
    """
    reports_dir = test_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    return reports_dir / filename


def validate_api_key(provider: str) -> Optional[str]:
    """
    Validate that API key for provider exists in environment.

    Args:
        provider: Provider name ('openrouter', 'openai', 'anthropic')

    Returns:
        API key if found, None otherwise

    Raises:
        EnvironmentError: If API key not found
    """
    key_map = {
        'openrouter': 'OPENROUTER_API_KEY',
        'openai': 'OPENAI_API_KEY',
        'anthropic': 'ANTHROPIC_API_KEY'
    }

    env_var = key_map.get(provider.lower())
    if not env_var:
        raise ValueError(f"Unknown provider: {provider}")

    api_key = os.getenv(env_var)
    if not api_key:
        raise EnvironmentError(
            f"{env_var} not found in environment. "
            f"Set it with: export {env_var}='your-api-key'"
        )

    return api_key
