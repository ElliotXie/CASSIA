# CASSIA Package Update and Build Script

This script (`update_and_build.py`) automates the entire process of updating, building, and installing the CASSIA package.

## Features

✅ **Automatic version management** - Updates both `setup.py` and `__init__.py`  
✅ **Flexible versioning** - Specify exact version or auto-increment  
✅ **Clean builds** - Optional cleanup of build artifacts  
✅ **Safe installation** - Uninstalls old version before installing new  
✅ **Verification** - Confirms installation and version after install  
✅ **Skip options** - Build-only mode for testing  

## Basic Usage

### 1. Auto-increment version (recommended)
```bash
# Increment patch version (0.2.1 → 0.2.2)
python update_and_build.py --increment patch

# Increment minor version (0.2.1 → 0.3.0)
python update_and_build.py --increment minor

# Increment major version (0.2.1 → 1.0.0)
python update_and_build.py --increment major
```

### 2. Specify exact version
```bash
# Set specific version
python update_and_build.py --version 0.2.5
python update_and_build.py --version 1.0.0
```

## Advanced Options

### Build without installing
```bash
# Build package but don't install (for testing)
python update_and_build.py --increment patch --skip-install
```

### Clean build artifacts first
```bash
# Remove build/, dist/, *.egg-info/ before building
python update_and_build.py --increment patch --clean
```

### Force reinstall same version
```bash
# Rebuild and reinstall current version
python update_and_build.py --version 0.2.1 --force-reinstall
```

## Complete Examples

### For regular development updates:
```bash
python update_and_build.py --increment patch
```
This will:
1. Increment version from 0.2.1 → 0.2.2
2. Update `setup.py` and `__init__.py`
3. Build the package
4. Uninstall old version
5. Install new version
6. Verify installation

### For testing builds:
```bash
python update_and_build.py --increment patch --skip-install --clean
```
This will:
1. Clean build artifacts
2. Increment version
3. Update version files
4. Build package
5. Skip installation (build only)

### For major releases:
```bash
python update_and_build.py --increment major --clean
```
This will:
1. Clean build artifacts
2. Increment version from 0.2.1 → 1.0.0
3. Update version files
4. Build and install

## Output

The script provides detailed progress information:

```
Working directory: D:\newgit\CASSIA\CASSIA_python
Current version: 0.2.1
New version: 0.2.2

Updating version from 0.2.1 to 0.2.2
Updating version in setup.py
  ✓ Updated pattern: version="0\.2\.1"
  ✓ Successfully updated setup.py
Updating version in CASSIA/__init__.py
  ✓ Updated pattern: __version__ = "0\.2\.1"
  ✓ Successfully updated CASSIA/__init__.py

Building package...
  ✓ Package built successfully

Installing package...
  ✓ Package uninstalled successfully
  ✓ Package installed successfully
  ✓ Installation verified - Version: 0.2.2

🎉 Successfully updated CASSIA to version 0.2.2!

Summary:
  Old version: 0.2.1
  New version: 0.2.2
  Wheel file: dist/cassia-0.2.2-py3-none-any.whl
  Status: Installed and verified ✓
```

## Error Handling

The script includes comprehensive error handling:
- Validates version format
- Checks for required files
- Handles build failures
- Verifies installation success
- Provides helpful error messages

## File Dependencies

The script expects these files to exist:
- `setup.py` - Contains package configuration and version
- `CASSIA/__init__.py` - Contains `__version__` variable
- Python build tools (`pip`, `build` package)

## Tips

1. **Always commit your changes before running the script** - The script modifies version files
2. **Use `--skip-install` for testing** - Build without installing to test changes
3. **Use `--clean` for clean builds** - Removes old build artifacts
4. **Check the Summary section** - Confirms what was done

## Troubleshooting

### "Could not determine current version"
- Check that `setup.py` exists and has a `version=` line
- Ensure the version format is `x.y.z` (e.g., "0.2.1")

### "Installation verification failed"
- Check that the package imported correctly
- Verify Python environment is correct
- Try running manually: `pip show CASSIA`

### Build failures
- Ensure `python -m build` package is installed: `pip install build`
- Check for syntax errors in Python files
- Review build output for specific error messages 