"""
CASSIA Test 22: Docker Installation Test
==========================================
Tests CASSIA R package installation in various Docker environments
to simulate fresh user experiences.

This test requires Docker to be installed and running.

Usage:
    python test_docker_install.py                    # Run all scenarios
    python test_docker_install.py --scenario conda   # Run specific scenario
    python test_docker_install.py --api-key KEY      # Include annotation test
"""

import sys
import time
import subprocess
from pathlib import Path

# Add shared utilities to path
sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import (
    load_config,
    setup_api_keys,
    print_test_header,
    print_test_result,
    get_test_mode
)
from result_manager import (
    create_results_dir,
    save_test_metadata,
    save_test_results,
    create_test_metadata,
    setup_logging,
    cleanup_logging
)


def check_docker_available() -> bool:
    """Check if Docker is available and running."""
    try:
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            timeout=30
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False


def run_docker_test():
    """Run Docker installation tests."""
    print_test_header("22 - Docker Installation Test")

    # Check Docker availability
    if not check_docker_available():
        print("\nERROR: Docker is not available or not running.")
        print("Please start Docker Desktop and try again.")
        return False

    print("\nDocker is available.")

    # Load configuration
    config = load_config()

    # Setup API keys (will be passed to containers)
    setup_api_keys()

    # Get API key from environment
    import os
    api_key = os.environ.get("OPENROUTER_API_KEY", "")

    # Create results directory
    results = create_results_dir("22_docker_install_test", get_test_mode())
    print(f"Results will be saved to: {results['base']}")

    # Setup logging
    logging_ctx = setup_logging(results['logs'])

    # Run the orchestrator
    start_time = time.time()
    errors = []
    test_results = {}

    try:
        script_path = Path(__file__).parent / "scripts" / "run_scenarios.py"

        cmd = [
            sys.executable,
            str(script_path),
            "--output", str(results['base']),
        ]

        if api_key:
            cmd.extend(["--api-key", api_key])
            print("\nAPI key found - will include annotation test")
        else:
            print("\nNo API key - skipping annotation test")

        print("\nRunning Docker scenarios...")
        print("This may take 10-30 minutes on first run (building images)")
        print("-" * 50)

        result = subprocess.run(
            cmd,
            capture_output=False,  # Let output stream to console
            timeout=3600  # 1 hour timeout for all scenarios
        )

        test_results["exit_code"] = result.returncode
        test_results["success"] = result.returncode == 0

    except subprocess.TimeoutExpired:
        errors.append("Test timed out after 1 hour")
        test_results["success"] = False
    except Exception as e:
        errors.append(f"Test failed: {str(e)}")
        test_results["success"] = False

    duration = time.time() - start_time

    # Cleanup logging
    cleanup_logging(logging_ctx)

    # Save metadata
    metadata = create_test_metadata(
        test_name="22_docker_install_test",
        config=config,
        duration_seconds=duration,
        status="passed" if test_results.get("success") else "failed",
        clusters_tested=[],
        errors=errors
    )
    save_test_metadata(results['base'], metadata)

    # Save results
    save_test_results(results['base'], test_results, "docker_test_results.json")

    # Print final result
    print("\n" + "=" * 50)
    if test_results.get("success"):
        print_test_result(True, f"Docker installation tests completed in {duration:.1f}s")
        return True
    else:
        print_test_result(False, f"Docker installation tests failed after {duration:.1f}s")
        for error in errors:
            print(f"  Error: {error}")
        return False


if __name__ == "__main__":
    success = run_docker_test()
    sys.exit(0 if success else 1)
