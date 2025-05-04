import os
import pandas as pd
import inspect

# Import the function from tools_function
from CASSIA.tools_function import iterative_marker_analysis

# Print the source file and signature of the function
print(f"Function source: {inspect.getfile(iterative_marker_analysis)}")
print(f"Function signature: {inspect.signature(iterative_marker_analysis)}")

# Check if the provider parameter is in the signature
sig = inspect.signature(iterative_marker_analysis)
if 'provider' in sig.parameters:
    print("SUCCESS: The 'provider' parameter exists in the function signature!")
else:
    print("ERROR: The 'provider' parameter is missing from the function signature.")
    print("Available parameters:", list(sig.parameters.keys())) 