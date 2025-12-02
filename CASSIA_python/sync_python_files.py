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
import filecmp
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
    def should_ignore(filename):
        """Check if a file/folder should be ignored."""
        if filename == "__pycache__":
            return True
        if filename.endswith(".pyc"):
            return True
        if filename.endswith(".ipynb"):
            return True
        if filename.startswith("test_") and filename.endswith(".py"):
            return True
        if filename.endswith("_test.py"):
            return True
        if filename.endswith(".code-workspace"):
            return True
        return False

    # Track changes
    unchanged_files = []
    changed_files = []
    new_files = []
    removed_files = []

    # Get all source files (excluding ignored patterns)
    def get_all_files(base_dir):
        """Get all files relative to base_dir, excluding ignored patterns."""
        files = []
        for root, dirs, filenames in os.walk(base_dir):
            # Filter out ignored directories
            dirs[:] = [d for d in dirs if not should_ignore(d)]
            for f in filenames:
                if not should_ignore(f):
                    full_path = Path(root) / f
                    rel_path = full_path.relative_to(base_dir)
                    files.append(rel_path)
        return set(files)

    source_files = get_all_files(source_dir)
    dest_files = get_all_files(dest_cassia) if dest_cassia.exists() else set()

    # Find new and removed files
    new_file_paths = source_files - dest_files
    removed_file_paths = dest_files - source_files
    common_files = source_files & dest_files

    # Compare common files
    for rel_path in common_files:
        src_file = source_dir / rel_path
        dst_file = dest_cassia / rel_path
        if filecmp.cmp(src_file, dst_file, shallow=False):
            unchanged_files.append(rel_path)
        else:
            changed_files.append(rel_path)

    new_files = list(new_file_paths)
    removed_files = list(removed_file_paths)

    # Now perform the actual sync
    print("Syncing files...")
    print()

    # Copy changed and new files
    for rel_path in changed_files + new_files:
        src_file = source_dir / rel_path
        dst_file = dest_cassia / rel_path
        dst_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src_file, dst_file)

    # Remove deleted files
    for rel_path in removed_files:
        dst_file = dest_cassia / rel_path
        if dst_file.exists():
            dst_file.unlink()

    # Print results with emojis
    print("=" * 60)
    print("SYNC RESULTS")
    print("=" * 60)
    print()

    if unchanged_files:
        print(f"✅ UNCHANGED ({len(unchanged_files)} files):")
        for f in sorted(unchanged_files):
            print(f"   ✅ {f}")
        print()

    if changed_files:
        print(f"🔄 CHANGED ({len(changed_files)} files):")
        for f in sorted(changed_files):
            print(f"   🔄 {f}")
        print()

    if new_files:
        print(f"✨ NEW ({len(new_files)} files):")
        for f in sorted(new_files):
            print(f"   ✨ {f}")
        print()

    if removed_files:
        print(f"🗑️  REMOVED ({len(removed_files)} files):")
        for f in sorted(removed_files):
            print(f"   🗑️  {f}")
        print()

    # Summary
    print("-" * 30)
    print("Summary:")
    print(f"  ✅ Unchanged: {len(unchanged_files)}")
    print(f"  🔄 Changed:   {len(changed_files)}")
    print(f"  ✨ New:       {len(new_files)}")
    print(f"  🗑️  Removed:  {len(removed_files)}")

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
