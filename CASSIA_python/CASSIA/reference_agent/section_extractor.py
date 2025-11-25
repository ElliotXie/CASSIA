"""
Section Extractor for Reference Agent.

Parses markdown reference files and extracts relevant sections
based on cell type, markers, and context.
"""

import re
from typing import Dict, List, Optional, Tuple
from pathlib import Path

try:
    from .utils import (
        load_markdown_file,
        parse_yaml_frontmatter,
        find_markers_in_text,
        format_reference_content
    )
except ImportError:
    from utils import (
        load_markdown_file,
        parse_yaml_frontmatter,
        find_markers_in_text,
        format_reference_content
    )


class MarkdownSection:
    """Represents a section in a markdown document."""

    def __init__(self, title: str, level: int, content: str, start_line: int = 0):
        self.title = title
        self.level = level  # 1 for #, 2 for ##, etc.
        self.content = content
        self.start_line = start_line
        self.children: List['MarkdownSection'] = []

    def __repr__(self):
        return f"MarkdownSection(title='{self.title}', level={self.level})"

    def get_full_content(self, include_children: bool = True) -> str:
        """Get the full content including optional children."""
        result = f"{'#' * self.level} {self.title}\n\n{self.content}"

        if include_children:
            for child in self.children:
                result += f"\n\n{child.get_full_content(include_children=True)}"

        return result


def parse_markdown(content: str) -> Tuple[Dict, List[MarkdownSection]]:
    """
    Parse markdown content into frontmatter and sections.

    Args:
        content: Full markdown content

    Returns:
        Tuple of (frontmatter_dict, list_of_sections)
    """
    frontmatter, body = parse_yaml_frontmatter(content)

    sections = _parse_sections(body)

    return frontmatter, sections


def _parse_sections(content: str) -> List[MarkdownSection]:
    """
    Parse markdown content into a list of sections.

    Args:
        content: Markdown content without frontmatter

    Returns:
        List of MarkdownSection objects (flat list)
    """
    sections = []
    lines = content.split('\n')

    current_section = None
    current_content_lines = []
    current_start_line = 0

    # Pattern to match markdown headers
    header_pattern = re.compile(r'^(#{1,6})\s+(.+)$')

    for i, line in enumerate(lines):
        match = header_pattern.match(line)

        if match:
            # Save previous section if exists
            if current_section is not None:
                current_section.content = '\n'.join(current_content_lines).strip()
                sections.append(current_section)

            # Start new section
            level = len(match.group(1))
            title = match.group(2).strip()
            current_section = MarkdownSection(title, level, "", start_line=i)
            current_content_lines = []
            current_start_line = i
        else:
            current_content_lines.append(line)

    # Don't forget the last section
    if current_section is not None:
        current_section.content = '\n'.join(current_content_lines).strip()
        sections.append(current_section)

    return sections


def build_section_tree(sections: List[MarkdownSection]) -> List[MarkdownSection]:
    """
    Build a hierarchical tree from flat sections based on header levels.

    Args:
        sections: Flat list of sections

    Returns:
        List of root-level sections with children populated
    """
    if not sections:
        return []

    root_sections = []
    stack: List[MarkdownSection] = []

    for section in sections:
        # Pop sections from stack that are same level or deeper
        while stack and stack[-1].level >= section.level:
            stack.pop()

        if stack:
            # This section is a child of the last section in stack
            stack[-1].children.append(section)
        else:
            # This is a root-level section
            root_sections.append(section)

        stack.append(section)

    return root_sections


def extract_sections(
    file_path: Path,
    cell_type_guess: Optional[str] = None,
    markers: Optional[List[str]] = None,
    depth: str = 'summary',
    include_sections: Optional[List[str]] = None,
    exclude_sections: Optional[List[str]] = None
) -> Dict:
    """
    Extract relevant sections from a reference markdown file.

    Args:
        file_path: Path to the markdown file
        cell_type_guess: Preliminary cell type classification
        markers: List of marker genes to look for
        depth: 'summary' for key points, 'detailed' for full content
        include_sections: Optional list of section names to always include
        exclude_sections: Optional list of section names to exclude

    Returns:
        Dict with:
            - frontmatter: Dict of YAML frontmatter
            - sections: List of extracted section contents
            - full_content: Combined content string
    """
    content = load_markdown_file(file_path)
    if not content:
        return {'frontmatter': {}, 'sections': [], 'full_content': ''}

    frontmatter, sections = parse_markdown(content)

    # Default sections to include
    default_include = ['overview', 'pitfalls', 'common pitfalls', 'technical notes']
    include_patterns = include_sections or []
    include_patterns.extend(default_include)

    exclude_patterns = exclude_sections or []

    extracted = []

    for section in sections:
        section_name_lower = section.title.lower()

        # Check exclusions first
        if any(ex.lower() in section_name_lower for ex in exclude_patterns):
            continue

        should_include = False

        # Check if section name matches include patterns
        if any(inc.lower() in section_name_lower for inc in include_patterns):
            should_include = True

        # Check if section mentions the guessed cell type
        if cell_type_guess:
            cell_type_lower = cell_type_guess.lower()
            if cell_type_lower in section_name_lower:
                should_include = True
            # Also check content for cell type mention
            if cell_type_lower in section.content.lower():
                should_include = True

        # Check if section mentions markers
        if markers:
            marker_mentions = find_markers_in_text(
                section.title + " " + section.content,
                markers[:10]  # Check top 10 markers
            )
            if len(marker_mentions) >= 2:
                should_include = True

        # For detailed depth, be more inclusive
        if depth == 'detailed':
            # Include more sections in detailed mode
            if section.level <= 2:  # Include top-level sections
                should_include = True

        if should_include:
            extracted.append(section)

    # Combine content
    full_content = _combine_sections(extracted, depth)

    return {
        'frontmatter': frontmatter,
        'sections': [{'title': s.title, 'content': s.content} for s in extracted],
        'full_content': full_content
    }


def _combine_sections(sections: List[MarkdownSection], depth: str) -> str:
    """
    Combine extracted sections into a single content string.

    Args:
        sections: List of sections to combine
        depth: 'summary' or 'detailed'

    Returns:
        Combined content string
    """
    if not sections:
        return ""

    parts = []
    seen_titles = set()

    for section in sections:
        # Avoid duplicates
        if section.title in seen_titles:
            continue
        seen_titles.add(section.title)

        if depth == 'summary':
            # For summary, include just the section header and first paragraph
            content_lines = section.content.split('\n\n')
            summary_content = content_lines[0] if content_lines else ""
            parts.append(f"### {section.title}\n{summary_content}")
        else:
            # For detailed, include full section
            parts.append(section.get_full_content(include_children=False))

    return '\n\n---\n\n'.join(parts)


def extract_sections_by_keywords(
    file_path: Path,
    keywords: List[str],
    max_sections: int = 5
) -> List[Dict]:
    """
    Extract sections that contain specific keywords.

    Args:
        file_path: Path to the markdown file
        keywords: List of keywords to search for
        max_sections: Maximum number of sections to return

    Returns:
        List of dicts with 'title', 'content', and 'score'
    """
    content = load_markdown_file(file_path)
    if not content:
        return []

    _, sections = parse_markdown(content)

    scored_sections = []

    for section in sections:
        full_text = (section.title + " " + section.content).lower()
        score = sum(1 for kw in keywords if kw.lower() in full_text)

        if score > 0:
            scored_sections.append({
                'title': section.title,
                'content': section.content,
                'score': score
            })

    # Sort by score (descending) and return top N
    scored_sections.sort(key=lambda x: x['score'], reverse=True)
    return scored_sections[:max_sections]


def get_section_by_path(file_path: Path, section_path: str) -> Optional[str]:
    """
    Get a specific section by its path (e.g., "Subtype Differentiation/Tregs").

    Args:
        file_path: Path to the markdown file
        section_path: Path to the section using / as separator

    Returns:
        Section content or None if not found
    """
    content = load_markdown_file(file_path)
    if not content:
        return None

    _, sections = parse_markdown(content)

    path_parts = [p.strip().lower() for p in section_path.split('/')]

    for section in sections:
        if section.title.lower() == path_parts[0]:
            if len(path_parts) == 1:
                return section.get_full_content()

            # Build tree and search deeper
            tree = build_section_tree(sections)
            return _search_section_tree(tree, path_parts)

    return None


def _search_section_tree(sections: List[MarkdownSection], path_parts: List[str]) -> Optional[str]:
    """Recursively search for a section in the tree."""
    if not path_parts or not sections:
        return None

    target = path_parts[0]

    for section in sections:
        if section.title.lower() == target:
            if len(path_parts) == 1:
                return section.get_full_content()
            return _search_section_tree(section.children, path_parts[1:])

    return None
