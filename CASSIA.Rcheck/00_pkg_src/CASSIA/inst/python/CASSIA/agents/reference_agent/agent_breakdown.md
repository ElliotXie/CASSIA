# Reference Agent Breakdown

## Purpose

The Reference Agent provides intelligent reference document retrieval and context injection for cell type annotation tasks. It uses a two-step ReAct workflow to assess marker complexity and select relevant reference documents.

## Architecture

```
reference_agent/
├── __init__.py              # Module exports
├── reference_agent.py       # Main orchestrator (ReferenceAgent class)
├── complexity_scorer.py     # Two-step LLM complexity assessment
├── reference_selector.py    # Reference scoring and selection
├── section_extractor.py     # Markdown parsing and section extraction
├── utils.py                 # Utility functions
└── references/              # Reference document library
    ├── _router.md           # Library structure map (shown to LLM)
    ├── index.json           # Reference metadata index
    ├── t_cell/              # T cell references
    │   ├── _overview.md
    │   ├── cd4/_overview.md
    │   └── cd8/_overview.md
    ├── b_cell/              # B cell references
    │   └── _overview.md
    └── myeloid/             # Myeloid cell references
        └── _overview.md
```

## Two-Step ReAct Workflow

### Step 1: Complexity Assessment (No Router)
- LLM receives marker genes + tissue/species context
- Outputs:
  - `preliminary_cell_type`: Best guess cell type
  - `cell_type_range`: List of possible cell types
  - `complexity_score`: 0-100 (higher = more ambiguous)
  - `requires_reference`: Boolean decision

### Step 2: Reference Selection (With Router)
- Only called if Step 1 returns `requires_reference=True`
- LLM sees the router structure (`_router.md`)
- Selects 1-3 most relevant reference file paths
- Returns paths like `t_cell/cd4/treg.md`

## Key Components

### ReferenceAgent (reference_agent.py)
Main class that orchestrates the workflow.

```python
agent = ReferenceAgent(provider="openrouter", model="google/gemini-2.5-flash")
result = agent.get_reference_for_markers(
    markers=["CD3D", "CD4", "FOXP3", "IL2RA"],
    tissue="blood",
    species="human"
)
```

**Returns:**
- `should_use_reference`: bool
- `content`: Extracted reference content
- `references_used`: List of file paths
- `complexity_score`: 0-100
- `preliminary_cell_type`: str
- `reasoning`: str

### complexity_scorer.py
- `assess_complexity()`: Main entry point (runs both steps)
- `assess_complexity_step1()`: Step 1 only
- `select_references_step2()`: Step 2 only
- `quick_complexity_check()`: Rule-based check without LLM

### reference_selector.py
- `select_references()`: Score and rank references
- `find_references_by_category()`: Get references by category
- `ReferenceCandidate`: Scoring container class

### section_extractor.py
- `parse_markdown()`: Parse markdown into sections
- `extract_sections()`: Extract relevant sections by cell type/markers
- `MarkdownSection`: Section container class

### utils.py
- File loading (JSON, markdown)
- YAML frontmatter parsing
- Marker normalization and matching
- Content formatting

## Reference Document Format

Each reference markdown file uses YAML frontmatter:

```yaml
---
id: t_cell_overview
category: t_cell
cell_types:
  - T cell
  - T lymphocyte
trigger_markers:
  - CD3D
  - CD3E
exclusion_markers: []
---

# T Cell Annotation Guide

## Overview
...

## Common Pitfalls
...
```

## Completion Status

| Component | Status |
|-----------|--------|
| reference_agent.py | Complete |
| complexity_scorer.py | Complete |
| reference_selector.py | Complete |
| section_extractor.py | Complete |
| utils.py | Complete |
| _router.md | Complete (structure) |
| index.json | Complete |
| Reference files | ~15% complete |

### Missing Reference Files
The router describes ~40 files but only 6 exist:
- T cell subtypes: treg.md, th1.md, th2.md, th17.md, tfh.md, cytotoxic.md, exhausted.md
- B cell subtypes: naive/, memory/, germinal_center/, plasma/
- Myeloid subtypes: monocyte/, macrophage/, dendritic/, granulocyte/
- NK cells: entire directory
- ILC: entire directory
- Epithelial, endothelial, stromal, progenitor: all missing

## Usage Example

```python
from CASSIA.reference_agent import get_reference_content, format_reference_for_prompt

# Get reference content for markers
result = get_reference_content(
    markers=["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4"],
    tissue="blood",
    species="human",
    provider="openrouter"
)

if result['should_use_reference']:
    # Format for prompt injection
    prompt_addition = format_reference_for_prompt(result)
    print(f"Complexity: {result['complexity_score']}")
    print(f"Preliminary: {result['preliminary_cell_type']}")
```

## Configuration

Default models by provider:
- openrouter: `google/gemini-2.5-flash`
- openai: `gpt-4o-mini`
- anthropic: `claude-3-haiku-20240307`

Fast models are used by default for cost efficiency.
