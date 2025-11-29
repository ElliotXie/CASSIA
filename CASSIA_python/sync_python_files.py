#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python File Sync Script for CASSIA

This script automatically syncs the CASSIA Python package from the CASSIA_python directory
to the CASSIA_R package directory.

The new approach copies the entire CASSIA/ folder to maintain the organized structure.
The R package imports the package via:
    py_cassia <- reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))

Usage:
    python sync_python_files.py

Author: CASSIA Development Team
"""

import os
import shutil
from pathlib import Path
from datetime import datetime


def sync_python_files():
    """
    Sync the CASSIA Python package to CASSIA_R/inst/python/

    This copies the entire CASSIA/ folder, preserving the organized structure:
    - core/
    - engine/
    - agents/
    - evaluation/
    - comparison/
    - hypothesis/
    - reports/
    - imaging/
    - pipeline/
    - config/
    - data/
    - reference_agent/
    - __init__.py (with all public exports)

    Excludes:
    - __pycache__ directories
    - .pyc files
    - Test files and notebooks
    """

    # Get the script directory and set up paths
    script_dir = Path(__file__).parent.resolve()

    # Source directory (the CASSIA package)
    source_dir = script_dir / "CASSIA"

    # Destination directory (R package python files)
    dest_dir = script_dir.parent / "CASSIA_R" / "inst" / "python"

    # Destination CASSIA folder
    dest_cassia = dest_dir / "CASSIA"

    print("CASSIA Python Package Sync")
    print("=" * 60)
    print(f"Source:      {source_dir}")
    print(f"Destination: {dest_cassia}")
    print(f"Time:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Check if directories exist
    if not source_dir.exists():
        print(f"Error: Source directory not found: {source_dir}")
        return False

    if not dest_dir.exists():
        print(f"Error: Destination directory not found: {dest_dir}")
        return False

    # Patterns to exclude during copy
    def ignore_patterns(directory, files):
        """Return files/folders to ignore during copy."""
        ignored = []
        for f in files:
            # Ignore __pycache__ directories
            if f == "__pycache__":
                ignored.append(f)
            # Ignore .pyc files
            elif f.endswith(".pyc"):
                ignored.append(f)
            # Ignore .ipynb files
            elif f.endswith(".ipynb"):
                ignored.append(f)
            # Ignore test files (but not test_data or test fixtures)
            elif f.startswith("test_") and f.endswith(".py"):
                ignored.append(f)
            elif f.endswith("_test.py"):
                ignored.append(f)
            # Ignore VS Code workspace files
            elif f.endswith(".code-workspace"):
                ignored.append(f)
        return ignored

    # Remove existing CASSIA folder in destination if it exists
    if dest_cassia.exists():
        print("Removing old CASSIA folder...")
        shutil.rmtree(dest_cassia)

    # Copy the entire CASSIA folder
    print("Copying CASSIA package...")
    shutil.copytree(source_dir, dest_cassia, ignore=ignore_patterns)

    # Count copied files
    file_count = sum(1 for _ in dest_cassia.rglob("*.py"))
    folder_count = sum(1 for p in dest_cassia.rglob("*") if p.is_dir())

    print()
    print("Sync Summary:")
    print("-" * 30)
    print(f"Python files copied: {file_count}")
    print(f"Folders created:     {folder_count}")
    print()

    # List the structure
    print("Package structure:")
    for item in sorted(dest_cassia.iterdir()):
        if item.is_dir():
            subfiles = list(item.rglob("*.py"))
            print(f"  {item.name}/ ({len(subfiles)} files)")
        elif item.suffix == ".py":
            print(f"  {item.name}")

    print()
    print("Sync completed successfully!")
    return True


def verify_sync():
    """
    Verify that the sync was successful by checking key files exist.
    """
    script_dir = Path(__file__).parent.resolve()
    dest_cassia = script_dir.parent / "CASSIA_R" / "inst" / "python" / "CASSIA"

    print()
    print("Verifying sync...")

    # Key files that should exist
    key_files = [
        "__init__.py",
        "core/__init__.py",
        "core/llm_utils.py",
        "core/validation.py",
        "engine/__init__.py",
        "engine/tools_function.py",
        "engine/main_function_code.py",
        "agents/annotation_boost/annotation_boost.py",
        "agents/merging/merging_annotation.py",
        "agents/uncertainty/Uncertainty_quantification.py",
        "agents/subclustering/subclustering.py",
        "evaluation/scoring.py",
        "pipeline/pipeline.py",
        "reports/generate_reports.py",
    ]

    all_good = True
    for rel_path in key_files:
        file_path = dest_cassia / rel_path
        if file_path.exists():
            print(f"  OK: {rel_path}")
        else:
            print(f"  MISSING: {rel_path}")
            all_good = False

    if all_good:
        print()
        print("All key files verified!")
    else:
        print()
        print("Warning: Some files may not have synced correctly")

    return all_good


def main():
    """Main function"""
    print()
    print("Starting CASSIA Python package sync...")
    print()

    try:
        # Perform the sync
        success = sync_python_files()

        if success:
            # Verify the sync
            verify_sync()
            print()
            print("All done! The CASSIA Python package has been synced to the R package.")
            print("The R package can now import it via:")
            print('  py_cassia <- reticulate::import_from_path("CASSIA", ...)')
        else:
            print()
            print("Sync completed with errors. Please check the messages above.")

    except KeyboardInterrupt:
        print()
        print("Sync interrupted by user")
    except Exception as e:
        print()
        print(f"Unexpected error: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
