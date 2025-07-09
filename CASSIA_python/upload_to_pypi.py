#!/usr/bin/env python3
"""
CASSIA Package Uploader for Official PyPI

This script automates the process of:
1. Building the CASSIA package with a specified version.
2. Uploading the package to the official PyPI repository.

Usage:
    python upload_to_pypi.py --version <new-version>

Example:
    python upload_to_pypi.py --version 0.2.5
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

def run_command(command, check=True):
    """Run a shell command."""
    print(f"Running: {command}")
    try:
        subprocess.run(command, shell=True, check=check)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function to build and upload the package to PyPI."""
    parser = argparse.ArgumentParser(
        description="Build and upload CASSIA package to the official PyPI."
    )
    parser.add_argument(
        '--version', 
        required=True, 
        help='The new version number to build and upload (e.g., 0.2.5).'
    )
    args = parser.parse_args()
    new_version = args.version

    # Ensure the script is run from the correct directory
    script_dir = Path(__file__).parent.resolve()
    os.chdir(script_dir)
    print(f"Working directory: {os.getcwd()}")

    # --- Confirmation Step ---
    print("\n" + "="*60)
    print("⚠️  WARNING: You are about to upload to the OFFICIAL PyPI repository! ⚠️")
    print("="*60)
    print(f"Package: CASSIA")
    print(f"Version: {new_version}")
    print("This action is permanent and cannot be undone.")
    
    confirm = input("Are you sure you want to continue? (yes/no): ")
    if confirm.lower() != 'yes':
        print("Upload cancelled by user.")
        sys.exit(0)

    # --- Step 1: Build the package ---
    print("\n--- Step 1: Building package ---")
    build_script_path = script_dir / "update_and_build.py"
    if not build_script_path.exists():
        print(f"Error: Build script not found at {build_script_path}", file=sys.stderr)
        sys.exit(1)

    build_command = (
        f'python "{build_script_path}" --version {new_version} '
        '--clean --skip-install --force-reinstall'
    )
    run_command(build_command)

    # --- Step 2: Verify build artifacts ---
    print("\n--- Step 2: Verifying build artifacts ---")
    dist_dir = script_dir / "dist"
    if not dist_dir.exists() or not any(dist_dir.iterdir()):
        print("Error: 'dist' directory is empty or does not exist.", file=sys.stderr)
        print("Build step may have failed.", file=sys.stderr)
        sys.exit(1)
    
    print("Build artifacts found in 'dist/' directory.")

    # --- Step 3: Check for twine ---
    print("\n--- Step 3: Checking for twine ---")
    try:
        run_command("python -m twine --version")
    except subprocess.CalledProcessError:
        print("Error: 'twine' is not installed.", file=sys.stderr)
        print("Please install it using: pip install twine", file=sys.stderr)
        sys.exit(1)

    # --- Step 4: Upload to PyPI ---
    print("\n--- Step 4: Uploading to PyPI ---")
    upload_command = "python -m twine upload dist/*"
    
    print("Twine will now ask for your PyPI username and an API token.")
    print("Username: __token__")
    print("Password: The API token value (starting with pypi-)")

    try:
        run_command(upload_command)
    except subprocess.CalledProcessError:
        print("\nUpload failed. Please check your credentials and network.", file=sys.stderr)
        sys.exit(1)

    # --- Success Message ---
    print("\n" + "="*50)
    print("🎉 Successfully uploaded to PyPI! 🎉")
    print("="*50)
    print("\nThe new version should be available shortly.")
    print(f"You can check it at: https://pypi.org/project/CASSIA/{new_version}/")

if __name__ == "__main__":
    main() 