---
title: Troubleshooting
---

Here we list some common errors users may encounter. If you have any other questions, feel free to leave a comment in the last section!

## Authentication Errors (Error 401)

### Common Causes
* API key not properly loaded
* Invalid or expired API key
* Insufficient API credits

### Solutions
```python
import os
# Verify API key is correctly set
key = os.environ.get("ANTHROPIC_API_KEY")
print(key)  # Should display your API key, not None or empty string

# Reset API key if needed
import CASSIA
CASSIA.set_api_key("your_anthropic_api_key", provider = "anthropic")
```

Make sure to:
* Verify you have sufficient credits on your platform account
* Add funds if necessary through the platform's billing portal

## File Not Found Errors

### Common Causes
* Incorrect file paths
* Missing input files
* CASSIA automatic file matching issues

### Solutions
```python
import os
# Verify file paths are correct and files exist
os.path.exists("your_input_file.csv")  # Returns True if file exists
```

When specifying paths:
* Use absolute paths when necessary:
  * Incorrect: "data/input.csv"
  * Correct: "C:/Users/YourName/Project/data/input.csv"
* Ensure input and output filenames match CASSIA's expectations:
  * Follow naming conventions exactly
  * Check file extensions match required format

## Permission Errors

### Common Causes
* File already exists and is locked/in use
* Insufficient write permissions
* Directory access restrictions

### Solutions
```python
import os
# If file exists, remove it before proceeding
if os.path.exists("output.csv"):
    os.remove("output.csv")
```

Additional steps:
* Check if file is already open in another program
* Close any programs using the file
* Wait a few seconds before retrying
* Use a different output filename if needed

## Best Practices
* Keep API keys secure and never share them
* Always backup important data before overwriting files
* Double-check file paths and permissions before running operations
* Maintain sufficient API credits for uninterrupted processing

