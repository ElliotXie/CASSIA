"""
Utility functions for the Reference Agent module.

Provides file loading, path management, and common helper functions.
"""

import json
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Any
import yaml


def get_references_dir() -> Path:
    """Get the path to the references directory."""
    return Path(__file__).parent / "references_brain"


def get_reference_index_path() -> Path:
    """Get the path to the reference index.json file."""
    return get_references_dir() / "index.json"


def load_json_file(file_path: Path) -> Dict:
    """
    Load a JSON file and return its contents.

    Args:
        file_path: Path to the JSON file

    Returns:
        Dictionary with file contents, empty dict if file doesn't exist
    """
    if not file_path.exists():
        return {}

    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def load_reference_index() -> Dict:
    """
    Load the reference index.json file.

    Returns:
        Dictionary with reference index, empty dict if not found
    """
    return load_json_file(get_reference_index_path())


def load_markdown_file(file_path: Path) -> str:
    """
    Load a markdown file and return its contents.

    Args:
        file_path: Path to the markdown file

    Returns:
        String with file contents, empty string if file doesn't exist
    """
    if not file_path.exists():
        return ""

    with open(file_path, 'r', encoding='utf-8') as f:
        return f.read()


def parse_yaml_frontmatter(content: str) -> tuple[Dict, str]:
    """
    Parse YAML frontmatter from markdown content.

    Args:
        content: Full markdown content with optional YAML frontmatter

    Returns:
        Tuple of (frontmatter_dict, remaining_content)
    """
    frontmatter = {}
    remaining_content = content

    # Check for YAML frontmatter (starts with ---)
    if content.startswith('---'):
        # Find the closing ---
        parts = content.split('---', 2)
        if len(parts) >= 3:
            try:
                frontmatter = yaml.safe_load(parts[1]) or {}
                remaining_content = parts[2].strip()
            except yaml.YAMLError:
                # If YAML parsing fails, return empty frontmatter
                pass

    return frontmatter, remaining_content


def list_reference_files(category: Optional[str] = None) -> List[Path]:
    """
    List all reference markdown files, optionally filtered by category.

    Args:
        category: Optional category folder to filter by (e.g., 't_cell', 'b_cell')

    Returns:
        List of paths to markdown files
    """
    refs_dir = get_references_dir()

    if category:
        search_dir = refs_dir / category
    else:
        search_dir = refs_dir

    if not search_dir.exists():
        return []

    # Find all .md files recursively
    return list(search_dir.rglob("*.md"))


def get_reference_path(reference_id: str) -> Optional[Path]:
    """
    Get the path to a reference file by its ID.

    Args:
        reference_id: The reference ID from the index

    Returns:
        Path to the reference file, or None if not found
    """
    index = load_reference_index()

    if 'references' not in index:
        return None

    ref_info = index['references'].get(reference_id)
    if not ref_info or 'path' not in ref_info:
        return None

    ref_path = get_references_dir() / ref_info['path']
    if ref_path.exists():
        return ref_path

    return None


def normalize_marker_name(marker: str) -> str:
    """
    Normalize a marker gene name for comparison.

    Args:
        marker: Raw marker name

    Returns:
        Normalized marker name (uppercase, stripped)
    """
    return marker.strip().upper()


def normalize_marker_list(markers: List[str]) -> List[str]:
    """
    Normalize a list of marker gene names.

    Args:
        markers: List of raw marker names

    Returns:
        List of normalized marker names
    """
    return [normalize_marker_name(m) for m in markers]


def find_markers_in_text(text: str, markers: List[str]) -> List[str]:
    """
    Find which markers from a list appear in the given text.

    Args:
        text: Text to search in
        markers: List of marker names to look for

    Returns:
        List of markers found in the text
    """
    found = []
    text_upper = text.upper()

    for marker in markers:
        # Create a pattern that matches the marker as a whole word
        pattern = r'\b' + re.escape(normalize_marker_name(marker)) + r'\b'
        if re.search(pattern, text_upper):
            found.append(marker)

    return found


def calculate_marker_overlap(markers1: List[str], markers2: List[str]) -> float:
    """
    Calculate the overlap ratio between two marker lists.

    Args:
        markers1: First list of markers
        markers2: Second list of markers

    Returns:
        Overlap ratio (0.0 to 1.0)
    """
    if not markers1 or not markers2:
        return 0.0

    set1 = set(normalize_marker_list(markers1))
    set2 = set(normalize_marker_list(markers2))

    intersection = len(set1 & set2)
    union = len(set1 | set2)

    return intersection / union if union > 0 else 0.0


def format_reference_content(content: str, max_length: Optional[int] = None) -> str:
    """
    Format reference content for injection into prompts.

    Args:
        content: Raw reference content
        max_length: Optional maximum length to truncate to

    Returns:
        Formatted content string
    """
    # Clean up excessive whitespace
    content = re.sub(r'\n{3,}', '\n\n', content)
    content = content.strip()

    # Truncate if needed
    if max_length and len(content) > max_length:
        content = content[:max_length] + "\n\n[... content truncated ...]"

    return content
