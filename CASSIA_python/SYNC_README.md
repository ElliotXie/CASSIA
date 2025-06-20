# CASSIA Python File Sync Tool

This tool automatically syncs Python files from the `CASSIA_python` development directory to the `CASSIA_R` package directory.

## 🎯 Purpose

When you make changes to Python files in the `CASSIA_python/CASSIA/` directory, you need to copy them to the `CASSIA_R/inst/python/` directory so the R package can use the updated code. This script automates that process.

## 📁 Directory Structure

```
CASSIA/
├── CASSIA_python/                  # ← Development directory (where you make changes)
│   ├── CASSIA/                     # ← Source Python files
│   │   ├── annotation_boost.py
│   │   ├── tools_function.py
│   │   ├── main_function_code.py
│   │   └── ...
│   └── sync_python_files.py       # ← This sync script
│
└── CASSIA_R/                       # ← R package directory
    └── inst/
        └── python/                 # ← Destination for Python files
            ├── annotation_boost.py
            ├── tools_function.py
            ├── merging_annotation_code.py  # ← NOT synced (R-specific)
            └── ...
```

## 🚀 Usage

### Basic Usage

Simply run the script from the `CASSIA_python` directory:

```bash
cd CASSIA/CASSIA_python
python sync_python_files.py
```

### What It Does

1. **Compares files** between source and destination
2. **Copies new files** that don't exist in the destination
3. **Updates modified files** that have changed
4. **Skips unchanged files** to save time
5. **Excludes specific files** that shouldn't be synced
6. **Verifies the sync** was successful

## 📋 Files Excluded from Sync

The following files are **NOT** synced and will be preserved in the R directory:

- `merging_annotation_code.py` - R-specific version with different functionality
- `__init__.py` - May contain R-specific imports
- Test files (`*_test.py`, `test_*`)
- Jupyter notebooks (`*.ipynb`)
- VS Code workspace files (`*.code-workspace`)
- Python cache files (`__pycache__`, `*.pyc`)

## 📊 Output Example

```
🚀 Starting CASSIA Python file sync...

🔄 CASSIA Python File Sync
==================================================
📁 Source:      /path/to/CASSIA_python/CASSIA
📁 Destination: /path/to/CASSIA_R/inst/python
⏰ Time:        2025-06-18 21:38:55

📋 Processing files:
------------------------------
🔄 UPDATE: annotation_boost.py          (modified)
➕ NEW: debug_genes.py                  (new file)
✅ SAME: llm_utils.py                   (no changes)
⏭️  SKIP: merging_annotation_code.py     (excluded file)

📊 Sync Summary:
--------------------
✅ Total files processed: 15
➕ New files copied:      2
🔄 Files updated:         1
⏭️  Files skipped:         3
❌ Errors:                0

🔍 Verifying sync...
✅ tools_function.py - synced correctly
✅ annotation_boost.py - synced correctly
✅ All key files verified successfully!

🎉 All done! Python files have been synced to the R package.
```

## 🔧 Customization

### Adding Files to Exclude

Edit the `excluded_files` set in `sync_python_files.py`:

```python
excluded_files = {
    "merging_annotation_code.py",  # R-specific version
    "__init__.py",                 # May have R-specific content
    "your_custom_file.py",         # Add your exclusions here
}
```

### Adding Pattern Exclusions

Edit the `excluded_patterns` set to exclude files matching certain patterns:

```python
excluded_patterns = {
    "__pycache__",
    ".pyc",
    ".ipynb",
    "custom_pattern",  # Add patterns here
}
```

## 🛠️ Troubleshooting

### Common Issues

1. **"Source directory not found"**
   - Make sure you're running the script from the `CASSIA_python` directory
   - Check that the `CASSIA` subdirectory exists

2. **"Destination directory not found"**
   - Ensure the `CASSIA_R` directory exists at the same level as `CASSIA_python`
   - Check that `CASSIA_R/inst/python/` directory structure is intact

3. **Permission Errors**
   - Ensure you have write permissions to the destination directory
   - Close any files that might be open in editors

### Manual Verification

You can manually verify the sync by comparing key files:

```bash
# Compare a specific file
diff CASSIA_python/CASSIA/tools_function.py CASSIA_R/inst/python/tools_function.py

# Or use your preferred file comparison tool
```

## 💡 Tips

- **Run after every significant change** to keep the R package up to date
- **Check the output** to see what files were updated
- **The script is safe to run multiple times** - it only copies when needed
- **No backup is needed** - the script preserves excluded files automatically

## 🆘 Support

If you encounter issues:

1. Check the error messages in the script output
2. Verify directory structure matches the expected layout
3. Ensure file permissions allow copying
4. Make sure no files are locked by other applications

---

**Happy syncing!** 🎉 