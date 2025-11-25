"""
CASSIA Test Runner - Run Individual Test
=========================================
Run a specific test by number or name.

Usage:
    python run_test.py 01
    python run_test.py 02_batch_annotation
    python run_test.py single_annotation
"""

import sys
import os
import subprocess
from pathlib import Path


def find_test_folder(test_id: str) -> Path:
    """Find the test folder matching the given ID or name."""
    test_root = Path(__file__).parent

    # Get all test folders (those starting with digits)
    test_folders = sorted([
        d for d in test_root.iterdir()
        if d.is_dir() and d.name[0].isdigit()
    ])

    # Try exact match first
    for folder in test_folders:
        if folder.name == test_id:
            return folder

    # Try matching by number prefix
    for folder in test_folders:
        if folder.name.startswith(test_id) or folder.name.startswith(f"{test_id}_"):
            return folder

    # Try matching by name substring
    for folder in test_folders:
        if test_id.lower() in folder.name.lower():
            return folder

    return None


def run_python_test(test_folder: Path) -> bool:
    """Run the Python test in the given folder."""
    # Find the test script
    test_scripts = list(test_folder.glob("test_*.py"))
    if not test_scripts:
        print(f"No Python test script found in {test_folder}")
        return False

    test_script = test_scripts[0]
    print(f"\nRunning: {test_script.name}")
    print("-" * 50)

    # Run the test
    result = subprocess.run(
        [sys.executable, str(test_script)],
        cwd=str(test_folder)
    )

    return result.returncode == 0


def list_available_tests():
    """List all available tests."""
    test_root = Path(__file__).parent

    test_folders = sorted([
        d for d in test_root.iterdir()
        if d.is_dir() and d.name[0].isdigit()
    ])

    print("\nAvailable Tests:")
    print("-" * 40)
    for folder in test_folders:
        print(f"  {folder.name}")
    print()


def main():
    if len(sys.argv) < 2:
        print("CASSIA Test Runner - Individual Test")
        print("=" * 40)
        print("\nUsage:")
        print("  python run_test.py <test_id>")
        print("\nExamples:")
        print("  python run_test.py 01")
        print("  python run_test.py batch")
        print("  python run_test.py 03_validator")
        list_available_tests()
        sys.exit(1)

    test_id = sys.argv[1]

    # Handle special cases
    if test_id in ["--list", "-l", "list"]:
        list_available_tests()
        sys.exit(0)

    # Find and run the test
    test_folder = find_test_folder(test_id)

    if not test_folder:
        print(f"Error: Test '{test_id}' not found")
        list_available_tests()
        sys.exit(1)

    print(f"\n{'='*50}")
    print(f"Running Test: {test_folder.name}")
    print(f"{'='*50}")

    success = run_python_test(test_folder)

    print(f"\n{'='*50}")
    if success:
        print("Test PASSED")
    else:
        print("Test FAILED")
    print(f"{'='*50}")

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
