#!/usr/bin/env python3
"""
CASSIA Package Uploader for TestPyPI

This script automates the process of:
1. Automatically incrementing the version number.
2. Building the CASSIA package with the new version.
3. Uploading the package to the TestPyPI repository for testing.

Usage:
    # Auto-increment patch version (e.g., 0.2.2 -> 0.2.3)
    python upload_to_testpypi.py

    # Increment minor version (e.g., 0.2.2 -> 0.3.0)
    python upload_to_testpypi.py --increment minor
    
    # Manually specify a version
    python upload_to_testpypi.py --version 0.3.0
"""

import argparse
import os
import re
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

def get_current_version():
    """Get the current version from setup.py."""
    setup_path = "setup.py"
    if not os.path.exists(setup_path):
        return None
    
    with open(setup_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    match = re.search(r'version=[\'"]([^\'"]+)[\'"]', content)
    if match:
        return match.group(1)
    return None

def increment_version(version, increment_type='patch'):
    """Increment version number based on semantic versioning, with support for .devN suffixes."""
    
    base_version = version
    dev_part = None
    
    # Handle .dev suffixes
    if '.dev' in version:
        parts = version.split('.dev')
        base_version = parts[0]
        dev_part = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

    try:
        major, minor, patch = map(int, base_version.split('.'))
        
        if increment_type == 'major':
            major += 1
            minor = 0
            patch = 0
        elif increment_type == 'minor':
            minor += 1
            patch = 0
        elif increment_type == 'patch':
            if dev_part is not None:
                # If it's a dev version, just increment the dev number
                dev_part += 1
            else:
                # If it's a release, increment the patch number
                patch += 1
        
        new_base = f"{major}.{minor}.{patch}"
        if dev_part is not None:
            return f"{new_base}.dev{dev_part}"
        else:
            return new_base
            
    except (ValueError, IndexError):
        print(f"Warning: Could not parse base version '{base_version}'. Appending '.1'.")
        return f"{version}.1"

def main():
    """Main function to build and upload the package to TestPyPI."""
    parser = argparse.ArgumentParser(
        description="Build and upload CASSIA package to TestPyPI. Automatically increments version."
    )
    parser.add_argument(
        '--increment', 
        choices=['major', 'minor', 'patch'],
        default='patch',
        help="Type of version increment (default: patch)."
    )
    parser.add_argument(
        '--version', 
        help='Manually specify a version number, overriding the increment logic.'
    )
    args = parser.parse_args()

    # Ensure the script is run from the correct directory
    script_dir = Path(__file__).parent.resolve()
    os.chdir(script_dir)
    print(f"Working directory: {os.getcwd()}")

    if args.version:
        new_version = args.version
        print(f"Using manually specified version: {new_version}")
    else:
        current_version = get_current_version()
        if not current_version:
            print("Error: Could not determine current version from setup.py", file=sys.stderr)
            sys.exit(1)
        
        new_version = increment_version(current_version, args.increment)
        print(f"Current version: {current_version}")
        print(f"New version will be: {new_version} (increment type: {args.increment})")

    # --- Step 1: Build the package using the existing script ---
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

    # --- Step 2: Verify build artifacts exist ---
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

    # --- Step 4: Upload to TestPyPI ---
    print("\n--- Step 4: Uploading to TestPyPI ---")
    upload_command = f'python -m twine upload --repository testpypi -u __token__ -p "{TEST_PYPI_TOKEN}" dist/*'
    
    print("Attempting to upload with embedded API token...")

    try:
        run_command(upload_command)
    except subprocess.CalledProcessError:
        print("\nUpload failed. Please check your credentials and network.", file=sys.stderr)
        sys.exit(1)

    # --- Success Message ---
    print("\n" + "="*50)
    print("🎉 Successfully uploaded to TestPyPI! 🎉")
    print("="*50)
    print("\nTo install and verify the package, run:")
    print(f"pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple CASSIA=={new_version}")

if __name__ == "__main__":
    main() 