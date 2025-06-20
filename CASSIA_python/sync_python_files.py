#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python File Sync Script for CASSIA

This script automatically syncs Python files from the CASSIA_python directory 
to the CASSIA_R package directory, excluding files that should not be overwritten.

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
    Sync Python files from CASSIA_python/CASSIA/ to CASSIA_R/inst/python/
    
    Excludes:
    - merging_annotation_code.py (R-specific version, should not be overwritten)
    - __pycache__ directories
    - .pyc files
    - Test files and notebooks (optional)
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
        "__init__.py",                 # May have R-specific content
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
    
    print("🔄 CASSIA Python File Sync")
    print("=" * 50)
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
    
    # Get list of Python files in source directory
    python_files = []
    for file_path in source_dir.iterdir():
        if file_path.is_file() and file_path.suffix == '.py':
            python_files.append(file_path)
    
    if not python_files:
        print("⚠️  No Python files found in source directory")
        return True
    
    # Statistics
    copied_files = []
    skipped_files = []
    updated_files = []
    error_files = []
    
    print("📋 Processing files:")
    print("-" * 30)
    
    for source_file in python_files:
        filename = source_file.name
        
        # Check if file should be excluded
        should_exclude = False
        
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
            print(f"⏭️  SKIP: {filename:<30} ({reason})")
            skipped_files.append(filename)
            continue
        
        # Destination file path
        dest_file = dest_dir / filename
        
        try:
            # Check if file needs updating
            needs_update = True
            if dest_file.exists():
                # Compare files to see if they're different
                if filecmp.cmp(source_file, dest_file, shallow=False):
                    needs_update = False
                    print(f"✅ SAME: {filename:<30} (no changes)")
                else:
                    print(f"🔄 UPDATE: {filename:<28} (modified)")
                    updated_files.append(filename)
            else:
                print(f"➕ NEW: {filename:<31} (new file)")
                copied_files.append(filename)
            
            # Copy the file if it needs updating
            if needs_update:
                shutil.copy2(source_file, dest_file)
                
        except Exception as e:
            print(f"❌ ERROR: {filename:<29} ({str(e)})")
            error_files.append(filename)
    
    # Print summary
    print()
    print("📊 Sync Summary:")
    print("-" * 20)
    print(f"✅ Total files processed: {len(python_files)}")
    print(f"➕ New files copied:      {len(copied_files)}")
    print(f"🔄 Files updated:         {len(updated_files)}")
    print(f"⏭️  Files skipped:         {len(skipped_files)}")
    print(f"❌ Errors:                {len(error_files)}")
    
    if copied_files:
        print(f"\n📋 New files: {', '.join(copied_files)}")
    
    if updated_files:
        print(f"\n📋 Updated files: {', '.join(updated_files)}")
    
    if skipped_files:
        print(f"\n📋 Skipped files: {', '.join(skipped_files)}")
    
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
    
    # Key files to verify (these should be identical)
    key_files = [
        "tools_function.py",
        "annotation_boost.py", 
        "main_function_code.py",
        "llm_utils.py"
    ]
    
    all_good = True
    for filename in key_files:
        source_file = source_dir / filename
        dest_file = dest_dir / filename
        
        if source_file.exists() and dest_file.exists():
            if filecmp.cmp(source_file, dest_file, shallow=False):
                print(f"✅ {filename} - synced correctly")
            else:
                print(f"❌ {filename} - files differ!")
                all_good = False
        elif source_file.exists():
            print(f"⚠️  {filename} - missing in destination")
            all_good = False
        else:
            print(f"ℹ️  {filename} - not found in source")
    
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

if __name__ == "__main__":
    main() 