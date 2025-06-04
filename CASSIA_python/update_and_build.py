#!/usr/bin/env python3
"""
CASSIA Package Update and Build Script

This script automates the process of:
1. Updating version numbers in setup.py and __init__.py
2. Building the package (both wheel and source distribution)
3. Uninstalling the old version
4. Installing the new version
5. Verifying the installation

Usage:
    python update_and_build.py --version 0.2.2
    python update_and_build.py --version 0.2.2 --skip-install
    python update_and_build.py --version 0.2.2 --force-reinstall
"""

import argparse
import os
import re
import subprocess
import sys
import shutil
from pathlib import Path


def run_command(command, check=True, capture_output=False):
    """Run a command and optionally capture its output."""
    print(f"Running: {command}")
    try:
        if capture_output:
            result = subprocess.run(command, shell=True, check=check, 
                                  capture_output=True, text=True)
            return result.stdout.strip()
        else:
            subprocess.run(command, shell=True, check=check)
            return None
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}")
        print(f"Error: {e}")
        if capture_output and e.stdout:
            print(f"Stdout: {e.stdout}")
        if capture_output and e.stderr:
            print(f"Stderr: {e.stderr}")
        raise


def update_version_in_file(filepath, old_version, new_version):
    """Update version string in a file."""
    print(f"Updating version in {filepath}")
    
    if not os.path.exists(filepath):
        print(f"Warning: File {filepath} does not exist")
        return False
    
    # Read the file
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Update version patterns
    patterns = [
        (rf'version="{re.escape(old_version)}"', f'version="{new_version}"'),
        (rf"version='{re.escape(old_version)}'", f"version='{new_version}'"),
        (rf'__version__ = "{re.escape(old_version)}"', f'__version__ = "{new_version}"'),
        (rf"__version__ = '{re.escape(old_version)}'", f"__version__ = '{new_version}'"),
    ]
    
    updated = False
    for pattern, replacement in patterns:
        if re.search(pattern, content):
            content = re.sub(pattern, replacement, content)
            updated = True
            print(f"  ✓ Updated pattern: {pattern}")
    
    if updated:
        # Write the updated content back
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"  ✓ Successfully updated {filepath}")
        return True
    else:
        print(f"  ⚠ No version patterns found in {filepath}")
        return False


def get_current_version():
    """Get the current version from setup.py."""
    setup_path = "setup.py"
    if not os.path.exists(setup_path):
        return None
    
    with open(setup_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Look for version patterns
    patterns = [
        r'version=[\'"]([\d\.]+)[\'"]',
        r'__version__\s*=\s*[\'"]([\d\.]+)[\'"]'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, content)
        if match:
            return match.group(1)
    
    return None


def clean_build_artifacts():
    """Clean up build artifacts."""
    print("Cleaning build artifacts...")
    
    directories_to_clean = ['build', 'dist', 'CASSIA.egg-info']
    for dir_name in directories_to_clean:
        if os.path.exists(dir_name):
            print(f"  Removing {dir_name}/")
            shutil.rmtree(dir_name)
    
    print("  ✓ Build artifacts cleaned")


def build_package():
    """Build the package."""
    print("Building package...")
    run_command("python -m build")
    print("  ✓ Package built successfully")


def get_wheel_file(version):
    """Get the path to the built wheel file."""
    wheel_pattern = f"dist/cassia-{version}-py3-none-any.whl"
    if os.path.exists(wheel_pattern):
        return wheel_pattern
    
    # Look for any wheel file in dist directory
    dist_dir = Path("dist")
    if dist_dir.exists():
        wheel_files = list(dist_dir.glob("*.whl"))
        if wheel_files:
            return str(wheel_files[0])
    
    return None


def uninstall_package():
    """Uninstall the current CASSIA package."""
    print("Uninstalling current CASSIA package...")
    try:
        run_command("pip uninstall CASSIA -y")
        print("  ✓ Package uninstalled successfully")
        return True
    except subprocess.CalledProcessError:
        print("  ⚠ No existing package to uninstall (or uninstall failed)")
        return False


def install_package(wheel_file):
    """Install the package from wheel file."""
    print(f"Installing package from {wheel_file}...")
    run_command(f"pip install {wheel_file}")
    print("  ✓ Package installed successfully")


def verify_installation(expected_version):
    """Verify the installation."""
    print("Verifying installation...")
    try:
        # Test import and version
        version_output = run_command(
            'python -c "import CASSIA; print(CASSIA.__version__)"', 
            capture_output=True
        )
        
        if version_output == expected_version:
            print(f"  ✓ Installation verified - Version: {version_output}")
            return True
        else:
            print(f"  ✗ Version mismatch - Expected: {expected_version}, Got: {version_output}")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Installation verification failed: {e}")
        return False


def increment_version(version, increment_type='patch'):
    """Increment version number."""
    parts = version.split('.')
    if len(parts) != 3:
        raise ValueError("Version must be in format x.y.z")
    
    major, minor, patch = map(int, parts)
    
    if increment_type == 'major':
        major += 1
        minor = 0
        patch = 0
    elif increment_type == 'minor':
        minor += 1
        patch = 0
    elif increment_type == 'patch':
        patch += 1
    else:
        raise ValueError("increment_type must be 'major', 'minor', or 'patch'")
    
    return f"{major}.{minor}.{patch}"


def main():
    parser = argparse.ArgumentParser(description="Update and build CASSIA package")
    parser.add_argument('--version', help='New version number (e.g., 0.2.2)')
    parser.add_argument('--increment', choices=['major', 'minor', 'patch'], 
                       help='Auto-increment version (alternative to --version)')
    parser.add_argument('--skip-install', action='store_true', 
                       help='Skip installation after building')
    parser.add_argument('--force-reinstall', action='store_true', 
                       help='Force reinstall even if same version')
    parser.add_argument('--clean', action='store_true', 
                       help='Clean build artifacts before building')
    
    args = parser.parse_args()
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    print(f"Working directory: {os.getcwd()}")
    
    # Get current version
    current_version = get_current_version()
    if not current_version:
        print("Error: Could not determine current version from setup.py")
        sys.exit(1)
    
    print(f"Current version: {current_version}")
    
    # Determine new version
    if args.version:
        new_version = args.version
    elif args.increment:
        new_version = increment_version(current_version, args.increment)
    else:
        print("Error: Must specify either --version or --increment")
        sys.exit(1)
    
    print(f"New version: {new_version}")
    
    # Check if we need to update version
    if new_version == current_version and not args.force_reinstall:
        print("Version unchanged. Use --force-reinstall to rebuild same version.")
        response = input("Continue anyway? (y/N): ")
        if response.lower() != 'y':
            print("Aborted.")
            sys.exit(0)
    
    try:
        # Clean build artifacts if requested
        if args.clean:
            clean_build_artifacts()
        
        # Update version numbers
        if new_version != current_version:
            print(f"\nUpdating version from {current_version} to {new_version}")
            
            # Update setup.py
            setup_updated = update_version_in_file("setup.py", current_version, new_version)
            
            # Update __init__.py
            init_updated = update_version_in_file("CASSIA/__init__.py", current_version, new_version)
            
            if not (setup_updated or init_updated):
                print("Warning: No version files were updated")
        
        # Build package
        print(f"\nBuilding package...")
        build_package()
        
        # Find the built wheel file
        wheel_file = get_wheel_file(new_version)
        if not wheel_file:
            print("Error: Could not find built wheel file")
            sys.exit(1)
        
        print(f"Built wheel: {wheel_file}")
        
        # Install package if not skipped
        if not args.skip_install:
            print(f"\nInstalling package...")
            
            # Uninstall current version
            uninstall_package()
            
            # Install new version
            install_package(wheel_file)
            
            # Verify installation
            if verify_installation(new_version):
                print(f"\n🎉 Successfully updated CASSIA to version {new_version}!")
            else:
                print(f"\n❌ Installation verification failed")
                sys.exit(1)
        else:
            print(f"\n✓ Package built successfully. Install manually with:")
            print(f"  pip install {wheel_file}")
        
        print(f"\nSummary:")
        print(f"  Old version: {current_version}")
        print(f"  New version: {new_version}")
        print(f"  Wheel file: {wheel_file}")
        
        if not args.skip_install:
            print(f"  Status: Installed and verified ✓")
        else:
            print(f"  Status: Built only (installation skipped)")
            
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 