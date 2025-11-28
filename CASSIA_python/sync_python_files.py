#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python File Sync Script for CASSIA

This script automatically syncs Python files from the CASSIA_python directory
to the CASSIA_R package directory, excluding files that should not be overwritten.

Updated to handle the new organized folder structure - collects files from
subdirectories and flattens them for the R package.

Usage:
    python sync_python_files.py

Author: CASSIA Development Team
"""

import os
import shutil
import filecmp
from pathlib import Path
from datetime import datetime

# Mapping of subdirectories to search for Python files
# Format: subdirectory -> list of files to sync (or None for all .py files)
SUBDIRECTORY_MAP = {
    "core": None,  # All .py files
    "engine": None,
    "evaluation": None,
    "comparison": None,
    "hypothesis": None,
    "reports": None,
    "imaging": None,
    "pipeline": None,
    "config": None,
    "agents/annotation_boost": None,
    "agents/uncertainty": None,
    "agents/merging": None,
    "agents/subclustering": None,
    # reference_agent is handled separately (has its own folder structure)
}

def sync_python_files():
    """
    Sync Python files from CASSIA_python/CASSIA/ to CASSIA_R/inst/python/

    Now supports the new organized folder structure:
    - Collects files from subdirectories (core/, engine/, agents/, etc.)
    - Flattens them to the R package's flat structure

    Excludes:
    - merging_annotation_code.py (R-specific version, should not be overwritten)
    - __init__.py files (R package has its own)
    - __pycache__ directories
    - .pyc files
    - Test files and notebooks
    """

    # Get the script directory and set up paths
    script_dir = Path(__file__).parent.resolve()

    # Source directory (where you make changes)
    source_dir = script_dir / "CASSIA"

    # Destination directory (R package python files)
    dest_dir = script_dir.parent / "CASSIA_R" / "inst" / "python"

    # Files to exclude from syncing
    excluded_files = {
        "merging_annotation_code.py",  # R-specific version
        "__init__.py",                 # R package has its own
    }

    # File patterns to exclude
    excluded_patterns = {
        "__pycache__",
        ".pyc",
        ".ipynb",         # Jupyter notebooks
        ".code-workspace", # VS Code workspace files
        "_test.py",       # Test files
        "test_",          # Test files
    }

    print("🔄 CASSIA Python File Sync (with subdirectory support)")
    print("=" * 60)
    print(f"📁 Source:      {source_dir}")
    print(f"📁 Destination: {dest_dir}")
    print(f"⏰ Time:        {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    # Check if directories exist
    if not source_dir.exists():
        print(f"❌ Error: Source directory not found: {source_dir}")
        return False

    if not dest_dir.exists():
        print(f"❌ Error: Destination directory not found: {dest_dir}")
        return False

    # Collect all Python files from subdirectories
    python_files = []  # List of (source_path, dest_filename) tuples

    # 1. First, collect from subdirectories
    print("📂 Scanning subdirectories...")
    for subdir, _ in SUBDIRECTORY_MAP.items():
        subdir_path = source_dir / subdir
        if subdir_path.exists():
            for file_path in subdir_path.iterdir():
                if file_path.is_file() and file_path.suffix == '.py':
                    python_files.append((file_path, file_path.name))
            print(f"   ✓ {subdir}/")
        else:
            print(f"   ⚠ {subdir}/ (not found)")

    # 2. Also check reference_agent folder (special handling - preserve structure)
    ref_agent_src = source_dir / "agents" / "reference_agent"
    ref_agent_dest = dest_dir / "reference_agent"
    if ref_agent_src.exists():
        # Create reference_agent folder in dest if needed
        if not ref_agent_dest.exists():
            ref_agent_dest.mkdir(parents=True, exist_ok=True)
        for file_path in ref_agent_src.iterdir():
            if file_path.is_file() and file_path.suffix == '.py':
                # For reference_agent, preserve subfolder structure
                python_files.append((file_path, f"reference_agent/{file_path.name}"))
        print(f"   ✓ agents/reference_agent/")

    # 3. Also check the old reference_agent at root (if it exists)
    ref_agent_root = source_dir / "reference_agent"
    if ref_agent_root.exists() and ref_agent_root.is_dir():
        for file_path in ref_agent_root.iterdir():
            if file_path.is_file() and file_path.suffix == '.py':
                python_files.append((file_path, f"reference_agent/{file_path.name}"))
        print(f"   ✓ reference_agent/ (root)")

    print()

    if not python_files:
        print("⚠️  No Python files found in subdirectories")
        return True

    # Statistics
    copied_files = []
    skipped_files = []
    updated_files = []
    same_files = []
    error_files = []

    print("📋 Processing files:")
    print("-" * 50)

    for source_file, dest_filename in python_files:
        filename = source_file.name

        # Check if file should be excluded
        should_exclude = False
        reason = ""

        # Check exact filename exclusions
        if filename in excluded_files:
            should_exclude = True
            reason = "excluded file"

        # Check pattern exclusions
        for pattern in excluded_patterns:
            if pattern in filename:
                should_exclude = True
                reason = f"matches pattern '{pattern}'"
                break

        if should_exclude:
            print(f"⏭️  SKIP: {dest_filename:<35} ({reason})")
            skipped_files.append(dest_filename)
            continue

        # Destination file path
        if "/" in dest_filename:
            # Subfolder file (like reference_agent)
            dest_file = dest_dir / dest_filename
            dest_file.parent.mkdir(parents=True, exist_ok=True)
        else:
            dest_file = dest_dir / dest_filename

        try:
            # Check if file needs updating
            needs_update = True
            if dest_file.exists():
                # Compare files to see if they're different
                if filecmp.cmp(source_file, dest_file, shallow=False):
                    needs_update = False
                    same_files.append(dest_filename)
                else:
                    print(f"🔄 UPDATE: {dest_filename:<33} (modified)")
                    updated_files.append(dest_filename)
            else:
                print(f"➕ NEW: {dest_filename:<36} (new file)")
                copied_files.append(dest_filename)

            # Copy the file if it needs updating
            if needs_update:
                shutil.copy2(source_file, dest_file)

        except Exception as e:
            print(f"❌ ERROR: {dest_filename:<34} ({str(e)})")
            error_files.append(dest_filename)

    # Print summary
    print()
    print("📊 Sync Summary:")
    print("-" * 30)
    print(f"✅ Total files processed: {len(python_files)}")
    print(f"➕ New files copied:      {len(copied_files)}")
    print(f"🔄 Files updated:         {len(updated_files)}")
    print(f"✓  Files unchanged:       {len(same_files)}")
    print(f"⏭️  Files skipped:         {len(skipped_files)}")
    print(f"❌ Errors:                {len(error_files)}")

    if copied_files:
        print(f"\n📋 New files: {', '.join(copied_files)}")

    if updated_files:
        print(f"\n📋 Updated files: {', '.join(updated_files)}")

    if error_files:
        print(f"\n⚠️  Error files: {', '.join(error_files)}")
        return False

    print(f"\n✅ Sync completed successfully!")
    return True

def verify_sync():
    """
    Verify that the sync was successful by comparing key files
    """
    script_dir = Path(__file__).parent.resolve()
    source_dir = script_dir / "CASSIA"
    dest_dir = script_dir.parent / "CASSIA_R" / "inst" / "python"

    print("\n🔍 Verifying sync...")

    # Key files to verify with their new locations
    key_files = [
        ("engine/tools_function.py", "tools_function.py"),
        ("agents/annotation_boost/annotation_boost.py", "annotation_boost.py"),
        ("engine/main_function_code.py", "main_function_code.py"),
        ("core/llm_utils.py", "llm_utils.py"),
        ("core/validation.py", "validation.py"),
        ("evaluation/scoring.py", "scoring.py"),
    ]

    all_good = True
    for src_rel, dest_name in key_files:
        source_file = source_dir / src_rel
        dest_file = dest_dir / dest_name

        if source_file.exists() and dest_file.exists():
            if filecmp.cmp(source_file, dest_file, shallow=False):
                print(f"✅ {dest_name} - synced correctly")
            else:
                print(f"❌ {dest_name} - files differ!")
                all_good = False
        elif source_file.exists():
            print(f"⚠️  {dest_name} - missing in destination")
            all_good = False
        else:
            print(f"ℹ️  {dest_name} - source not found at {src_rel}")

    if all_good:
        print("✅ All key files verified successfully!")
    else:
        print("⚠️  Some files may not have synced correctly")

    return all_good

def main():
    """Main function"""
    print("🚀 Starting CASSIA Python file sync...\n")

    try:
        # Perform the sync
        success = sync_python_files()

        if success:
            # Verify the sync
            verify_sync()
            print(f"\n🎉 All done! Python files have been synced to the R package.")
            print("💡 You can now use the updated Python functions in your R package.")
        else:
            print(f"\n💥 Sync completed with errors. Please check the messages above.")

    except KeyboardInterrupt:
        print(f"\n🛑 Sync interrupted by user")
    except Exception as e:
        print(f"\n💥 Unexpected error: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
