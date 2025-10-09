# CASSIA Development Documentation

This folder contains documentation and tools for CASSIA developers.

## Quick Development Workflow

### Python Development

#### Fast Install for Testing (Recommended)

**EASIEST METHOD - No cd Required!**

1. Open `dev_docs/quick_install_dev.py` in your code editor
2. Run it (VSCode: Right-click → "Run Python File" or F5)
3. Done! Script auto-finds CASSIA_python directory

**Or from command line:**

```bash
# Run from dev_docs folder - no need to navigate anywhere else
cd dev_docs
python quick_install_dev.py --editable

# Or with full path from project root
python dev_docs/quick_install_dev.py --editable
```

**Editable Mode Benefits:**
- Changes to `.py` files take effect immediately
- No need to reinstall after each change
- Just restart your Python session or reload modules
- Perfect for rapid development

**Common Usage Patterns:**

```bash
# First time setup - full install with test
python quick_install_dev.py --editable --test

# Quick reinstall without dependencies (very fast)
python quick_install_dev.py -e --no-deps

# Force reinstall everything
python quick_install_dev.py -e --force

# Skip uninstall step (fastest, use with editable mode)
python quick_install_dev.py -e --skip-uninstall
```

#### Testing Your Changes

After installation, test your changes:

```python
# Quick import test
python -c "import CASSIA; print(CASSIA.__version__)"

# Test specific functionality
python -c "import CASSIA; print(CASSIA.runCASSIA_pipeline.__doc__)"

# Run test suite (if available)
cd CASSIA_python/CASSIA/test_code
python run_full_test.py
```

#### Development Workflow Example

```bash
# 1. Make your changes to Python files
vim CASSIA_python/CASSIA/tools_function.py

# 2. If NOT in editable mode, reinstall quickly
python dev_docs/quick_install_dev.py -e --no-deps

# 3. Test in Python
python
>>> import CASSIA
>>> # Test your changes
>>> exit()

# 4. Repeat steps 1-3 as needed
```

### R Development

The R package wraps the Python code, so Python changes need to be synced:

```bash
# 1. Make changes to Python code
# 2. Sync Python files to R package
cd CASSIA_python
python sync_python_files.py

# 3. Rebuild R package
python update_and_build.py --version <current_version>

# 4. Or manually in R:
# devtools::document("CASSIA_R")
# devtools::install("CASSIA_R")
```

## Project Structure for Developers

```
CASSIA/
├── CASSIA_python/          # Python implementation (source of truth)
│   ├── CASSIA/             # Main package code
│   │   ├── *.py            # Core modules
│   │   └── data/           # Data files and configs
│   ├── setup.py            # Package configuration
│   ├── sync_python_files.py      # Sync to R package
│   ├── update_and_build.py       # Build and version management
│   └── upload_to_pypi.py         # PyPI publishing
│
├── CASSIA_R/               # R wrapper package
│   ├── R/                  # R wrapper functions
│   │   ├── annotator.R     # Main R interface
│   │   ├── utils.R         # R utilities
│   │   └── *.R             # Other R wrappers
│   ├── inst/python/        # Bundled Python source (auto-synced)
│   └── DESCRIPTION         # R package metadata
│
└── dev_docs/               # Developer documentation (you are here)
    ├── README.md           # This file
    └── quick_install_dev.py # Fast development install script
```

## Making Changes

### Adding a New Python Function

1. **Add function to appropriate module:**
   ```python
   # CASSIA_python/CASSIA/tools_function.py
   def my_new_function(param1, param2):
       """Your docstring here."""
       # Implementation
       return result
   ```

2. **Export in `__init__.py`:**
   ```python
   # CASSIA_python/CASSIA/__init__.py
   from .tools_function import my_new_function
   ```

3. **Test locally:**
   ```bash
   python dev_docs/quick_install_dev.py -e
   python -c "import CASSIA; CASSIA.my_new_function()"
   ```

4. **Add R wrapper (if needed):**
   ```R
   # CASSIA_R/R/utils.R
   #' @export
   my_new_function <- function(param1, param2) {
     py_tools$my_new_function(param1, param2)
   }
   ```

5. **Sync and rebuild:**
   ```bash
   cd CASSIA_python
   python sync_python_files.py  # Sync to R package
   ```

### Modifying Existing Functions

1. **Edit the Python source**
2. **If in editable mode:** Just restart Python - no reinstall needed!
3. **If not in editable mode:** Run `python dev_docs/quick_install_dev.py -e --no-deps`
4. **If R changes needed:** Sync and rebuild as above

### Adding Dependencies

1. **Update `setup.py`:**
   ```python
   install_requires=[
       "numpy>=1.21.0",
       "your-new-package>=1.0.0",  # Add here
   ]
   ```

2. **Reinstall with dependencies:**
   ```bash
   python dev_docs/quick_install_dev.py -e --force
   ```

3. **Update R DESCRIPTION if needed:**
   ```R
   Imports:
       reticulate,
       your-r-package
   ```

## Testing

### Python Testing

```bash
# Run existing tests
cd CASSIA_python/CASSIA/test_code
python run_full_test.py

# Quick import test
python -c "import CASSIA; print('Import successful')"

# Test specific function
python -c "from CASSIA import runCASSIA; print(runCASSIA.__doc__)"
```

### R Testing

```R
# In R console
library(CASSIA)

# Run example
markers <- loadExampleMarkers(processed = FALSE)
result <- runCASSIA_pipeline(
  output_file_name = "test",
  tissue = "Brain",
  species = "Human",
  marker = markers,
  max_workers = 2
)
```

## Version Management

### For Development
- Use `--editable` mode - version updates not needed for local testing
- Version only matters when releasing

### For Release

```bash
# Update version everywhere
cd CASSIA_python
python update_and_build.py --version 0.3.2

# This updates:
# - setup.py
# - __init__.py
# - Builds distribution
# - Installs new version
```

## Publishing Workflow

### Python (PyPI)

```bash
cd CASSIA_python

# Test on TestPyPI first
python upload_to_testpypi.py

# Then production
python upload_to_pypi.py
```

### R (GitHub)

```bash
# Update version in DESCRIPTION
# Commit and push to GitHub
# Users install via: devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

## Common Issues & Solutions

### "ModuleNotFoundError: No module named 'CASSIA'"
```bash
# Solution: Install in editable mode
python dev_docs/quick_install_dev.py -e
```

### "Changes not taking effect"
- **If in editable mode:** Restart Python session
- **If not editable:** Reinstall with `python dev_docs/quick_install_dev.py -e --no-deps`

### "ImportError: cannot import name 'new_function'"
```bash
# Check __init__.py has the export
# Reinstall
python dev_docs/quick_install_dev.py -e --force
```

### R package not finding Python modules
```R
# Rebuild R package environment
library(CASSIA)
setup_cassia_env()
```

### API key issues during testing
```python
# Set temporary key for testing
import os
os.environ["OPENROUTER_API_KEY"] = "your-test-key"

# Or use set_api_key
import CASSIA
CASSIA.set_api_key("your-test-key", provider="openrouter")
```

## Tips for Efficient Development

1. **Always use editable mode** (`-e` flag) - saves tons of time
2. **Use `--no-deps`** when dependencies haven't changed - much faster
3. **Keep a test script** handy for quick manual testing
4. **Use `--test` flag** to verify imports after install
5. **Sync to R only when needed** - Python development is faster
6. **Use version control** - commit before major changes

## Environment Setup

### Recommended Development Environment

```bash
# Create virtual environment for development
python -m venv cassia_dev
source cassia_dev/bin/activate  # On Windows: cassia_dev\Scripts\activate

# Install in editable mode
cd CASSIA/CASSIA_python
python ../dev_docs/quick_install_dev.py -e --test

# Install development tools
pip install pytest black flake8 mypy
```

### Multiple Python Versions

```bash
# Test with different Python versions
python3.9 dev_docs/quick_install_dev.py -e
python3.10 dev_docs/quick_install_dev.py -e
python3.11 dev_docs/quick_install_dev.py -e
```

## Getting Help

- **Project Maintainer:** Elliot Xie (xie227@wisc.edu)
- **GitHub Issues:** https://github.com/ElliotXie/CASSIA/issues
- **Main README:** See `../README.md` for user documentation
- **Project Spec:** See `../CLAUDE.md` for detailed architecture

---

*Last Updated: October 2025*
