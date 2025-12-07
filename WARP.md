# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

CASSIA is a multi-agent LLM framework for automated, interpretable single-cell RNA-seq cell type annotation. The project includes both Python (`CASSIA_python/`) and R (`CASSIA_R/`) implementations that work together via `reticulate`.

## Key Commands

### Python Package Development

Install the package in development mode:
```bash
cd CASSIA_python
pip install -e .
```

Run the main annotation pipeline:
```bash
cd CASSIA_python/CASSIA
python CASSIA_python_tutorial.py
```

Run evaluation/testing:
```bash
cd CASSIA_python/CASSIA
python LLM_evaluation_test.py
python run_full_test.py
```

### R Package Development

Install the R package from the local directory:
```bash
# From R console
devtools::install_local("CASSIA_R")

# Or using R CLI
Rscript -e "devtools::install_local('CASSIA_R')"
```

Run R tests:
```bash
# From R console
devtools::test("CASSIA_R")

# Or using R CLI
Rscript CASSIA_R/tests/testthat.R
```

### Test Suite Execution

Run individual tests:
```bash
cd Test
python run_test.py 01  # Run first test by number
python run_test.py batch  # Run by name
python run_test.py --list  # List all available tests
```

Run all tests:
```bash
cd Test
python run_all_tests.py
Rscript run_all_tests.R
```

## Architecture Overview

### Core Design Pattern: Multi-Agent LLM System

CASSIA uses a modular agent architecture where specialized agents handle different aspects of the annotation pipeline:

1. **Annotation Agent** (Core): Performs initial cell type prediction from marker genes
2. **Validator Agent**: Validates predictions against LLM reasoning
3. **Scoring Agent**: Quality-assesses annotations (0-100 scale)
4. **Annotation Boost Agent**: Refines low-scoring annotations iteratively
5. **Merging Agent**: Groups related cell types across clusters
6. **Reference Agent**: Intelligent retrieval-augmented generation for complex cases
7. **Uncertainty Quantification Agent**: Consensus scoring from multiple runs

### Python Package Structure

```
CASSIA_python/CASSIA/
├── core/                    # Shared utilities and validation
│   ├── llm_utils.py        # LLM API calls (OpenAI, Anthropic, OpenRouter)
│   ├── validation.py       # Input validation and error handling
│   ├── model_settings.py   # Model configuration and fuzzy alias resolution
│   ├── marker_utils.py     # Marker gene processing and ranking
│   └── logging_config.py   # CASSIA logging system
├── engine/
│   ├── tools_function.py   # Main entry points (runCASSIA, runCASSIA_batch, etc.)
│   └── main_function_code.py  # Core analysis logic per provider
├── agents/                  # Specialized agent implementations
│   ├── annotation_boost/    # Iterative refinement for low-scoring clusters
│   ├── merging/             # Cell type grouping and consolidation
│   ├── subclustering/       # Sub-cluster analysis
│   ├── uncertainty/         # Multi-run consensus and variance
│   └── reference_agent/     # RAG system for complex cases
├── evaluation/              # Scoring and analysis utilities
│   ├── scoring.py          # Batch scoring functionality
│   └── cell_type_comparison.py  # Cell type comparison logic
├── pipeline/
│   └── pipeline.py         # End-to-end workflow orchestration
└── reports/                 # Report generation
    ├── generate_reports.py  # HTML report creation
    └── generate_batch_report.py  # Batch report aggregation
```

### Data Flow

1. **Input**: Marker genes → DataFrame with columns [cluster/cell_type, genes]
2. **Processing**: 
   - Extract top N genes per cluster (configurable ranking method)
   - For each cluster, call annotation agent with LLM
   - Optional: Apply reference agent if markers indicate complex cell type
3. **Output**: Structured results with main_cell_type, sub_cell_types, confidence metrics
4. **Quality Control**: Score batch, identify low-confidence clusters, optionally boost

### Provider Architecture

The package supports multiple LLM providers through pluggable provider handlers:
- **OpenAI**: `run_cell_type_analysis()` - Uses gpt-4o and variants
- **Anthropic**: `run_cell_type_analysis_claude()` - Uses Claude models  
- **OpenRouter**: `run_cell_type_analysis_openrouter()` - Meta-routing for 100+ models
- **Custom OpenAI-compatible**: `run_cell_type_analysis_custom()` - Self-hosted LLMs

Model names support fuzzy aliasing (e.g., "gpt" resolves to full name) via `ModelSettings.resolve_model_name()`.

### R Package Integration

The R package (`CASSIA_R/`) wraps Python functions via `reticulate`:

```
CASSIA_R/
├── R/
│   ├── annotator.R      # Main annotation functions (runCASSIA, runCASSIA_batch, etc.)
│   ├── data.R           # Built-in marker datasets
│   ├── utils.R          # R utility functions
│   ├── model_settings.R # Model configuration
│   └── merge_annotations.R  # Annotation merging
├── inst/python/CASSIA/  # Python package (symlinked from CASSIA_python)
└── DESCRIPTION          # Package metadata
```

Environment setup is automatic on first load via `.onLoad()`, creating isolated Python environment (virtualenv or conda) as needed.

## Key Development Concepts

### Backward Compatibility

The `__init__.py` exports from multiple locations to maintain backward compatibility:
- Functions can be imported as `from CASSIA import runCASSIA`
- R package via `reticulate` accesses Python module directly
- Functions are re-exported even when moved to submodules

### Input Validation

All main entry points use centralized validators in `core/validation.py`:
- `validate_runCASSIA_inputs()` - Single annotation validation
- `validate_runCASSIA_batch_inputs()` - Batch processing validation
- Early fail-fast on invalid inputs with actionable error messages

### Marker Gene Processing

Flexible marker input handling:
- String lists, comma-separated, data frames with configurable column names
- Top-N ranking by: avg_log2FC, p_val_adj, pct_diff, Score
- Customizable ascending/descending sort per method

### Parallel Processing

Batch operations use `ThreadPoolExecutor` with configurable workers:
- Per-cluster annotation in parallel threads
- Built-in retry logic with exponential backoff
- Progress tracking via `BatchProgressTracker` class
- Failure isolation: one failed cluster doesn't block others

### Report Generation

Multi-format output:
- CSV files: Full detailed results and summary
- HTML reports: Interactive tables with modal dialogs for conversation history
- Structured data preservation: JSON serialization within CSV columns
- Multiple report styles: Per-iteration vs. total summary

### Extended Metadata

Reference information tracking:
- `reference_used` (bool): Whether RAG was applied
- `complexity_score` (0-100): Marker complexity assessment
- `preliminary_cell_type`: Initial prediction before reference lookup
- `references_used` (list): Which knowledge bases were retrieved

## Handling API Keys

Set API keys before analysis:
```python
CASSIA.set_api_key(key, provider)  # Python
CASSIA::setLLMApiKey(key, provider)  # R
```

Supported providers: "openai", "anthropic", "openrouter", or custom HTTP URL for self-hosted.

## Branch Information

The current development branch is `feature/shared-agent-class`. Reference agent updates and architectural refactoring are in progress.
