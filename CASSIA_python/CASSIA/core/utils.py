"""
General utility functions for CASSIA.

This module provides common utility functions used across the CASSIA package
including data handling, formatting, and CSV operations.
"""

import re
import csv
import os

# Try to import extract_json_from_reply for rerun_formatting_agent
try:
    from CASSIA.engine.main_function_code import extract_json_from_reply
except ImportError:
    try:
        from ..engine.main_function_code import extract_json_from_reply
    except ImportError:
        try:
            from main_function_code import extract_json_from_reply
        except ImportError:
            extract_json_from_reply = None


def check_formatted_output(structured_output):
    """Check if the structured output has required fields."""
    return 'main_cell_type' in structured_output and 'sub_cell_types' in structured_output


def rerun_formatting_agent(agent, full_conversation_history):
    """Re-run the formatting agent on the full conversation history."""
    full_text = "\n\n".join([f"{role}: {message}" for role, message in full_conversation_history])
    formatted_result = agent(full_text, "user")
    if extract_json_from_reply is None:
        raise ImportError("extract_json_from_reply not available")
    return extract_json_from_reply(formatted_result)


def safe_get(dict_obj, *keys):
    """
    Safely get nested dictionary values.

    Args:
        dict_obj: Dictionary to traverse
        *keys: Keys to follow in sequence

    Returns:
        The value at the nested key path, or None if any key is missing
    """
    for key in keys:
        if isinstance(dict_obj, dict) and key in dict_obj:
            dict_obj = dict_obj[key]
        else:
            return None
    return dict_obj


def natural_sort_key(cell_type):
    """
    Create a sort key that handles numeric cluster names properly.

    Handles various formats:
    - Pure numbers: "0", "1", "10" → sorted as integers 0, 1, 10 (priority 0)
    - "cluster X": "cluster 0", "cluster 1", "cluster 10" → sorted by X numerically (priority 1)
    - "Cluster X": case-insensitive
    - Other text: sorted alphabetically (priority 2)

    Args:
        cell_type (str): The cell type or cluster name

    Returns:
        tuple: (sort_priority, numeric_value, string_value) for proper sorting
    """
    if not cell_type or not isinstance(cell_type, str):
        return (3, 0, str(cell_type))  # Non-string values go last

    cell_type_str = str(cell_type).strip()

    # Try to parse as pure integer
    try:
        return (0, int(cell_type_str), "")  # Pure numbers have priority 0
    except ValueError:
        pass

    # Try to extract number from "cluster X" or "Cluster X" format (case-insensitive)
    cluster_match = re.match(r'^cluster\s+(\d+)$', cell_type_str, re.IGNORECASE)
    if cluster_match:
        cluster_num = int(cluster_match.group(1))
        return (1, cluster_num, "")  # Cluster numbers have priority 1 (after pure numbers)

    # For any other text, sort alphabetically (priority 2)
    return (2, 0, cell_type_str.lower())


def clean_conversation_history(history_text):
    """
    Clean conversation history for safe CSV storage while preserving full content.

    Args:
        history_text (str): Raw conversation history text

    Returns:
        str: Cleaned text safe for CSV storage (no truncation)
    """
    if not history_text:
        return ""

    # Replace newlines with spaces (prevents row breaks in CSV/Excel)
    cleaned = history_text.replace('\n', ' ').replace('\r', ' ')

    # Collapse multiple spaces into single spaces
    cleaned = ' '.join(cleaned.split())

    # Double quotes will be handled by csv.writer's automatic escaping
    # No need to replace quotes - csv module handles this correctly

    # Return full content without truncation
    return cleaned


def write_csv(filename, headers, row_data):
    """
    Write data to a CSV file.

    Args:
        filename (str): Path to the output CSV file
        headers (list): List of column headers
        row_data (list): List of rows (each row is a list)
    """
    # Make sure the directory exists
    output_dir = os.path.dirname(filename)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(headers)
        writer.writerows(row_data)


def get_column_value(row, column_options, default=None):
    """
    Get a value from a row using multiple possible column names.

    This handles the common pattern where column names may vary (e.g.,
    'Conversation History' vs 'Conversation.History' vs 'conversation_history').

    Args:
        row: A pandas Series or dict-like object with column access
        column_options (list): List of possible column names to try
        default: Value to return if no column is found (default: None)

    Returns:
        The value from the first matching column, or default if none found

    Raises:
        KeyError: If no matching column found and default is not provided

    Example:
        >>> history = get_column_value(row, ['Conversation History', 'Conversation.History'])
    """
    for col in column_options:
        if col in row:
            return row[col]

    if default is not None:
        return default

    # Build helpful error message
    if hasattr(row, 'index'):
        available = list(row.index)
    elif hasattr(row, 'keys'):
        available = list(row.keys())
    else:
        available = "unknown"

    raise KeyError(f"Could not find any of {column_options}. Available columns: {available}")


# Common column name variants for convenience
MARKER_COLUMN_OPTIONS = ['Marker List', 'Marker.List', 'marker_list', 'Marker_List']
HISTORY_COLUMN_OPTIONS = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History']
REASONING_COLUMN_OPTIONS = ['Scoring_Reasoning', 'Scoring.Reasoning', 'scoring_reasoning', 'Scoring_reasoning']
CLUSTER_ID_COLUMN_OPTIONS = ['Cluster ID', 'Cluster.ID', 'cluster_id', 'Cluster_ID']
