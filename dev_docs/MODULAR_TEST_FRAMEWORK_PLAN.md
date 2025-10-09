# CASSIA Modular Testing Framework - Comprehensive Plan

**Author**: Claude Code (with Elliot Yixuan Xie)
**Created**: 2025-10-07
**Status**: Implementation Ready
**Version**: 1.0

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Current State Analysis](#current-state-analysis)
3. [Proposed Architecture](#proposed-architecture)
4. [Directory Structure](#directory-structure)
5. [Test Module Specifications](#test-module-specifications)
6. [Shared Infrastructure](#shared-infrastructure)
7. [Standard Templates](#standard-templates)
8. [Implementation Roadmap](#implementation-roadmap)
9. [Success Criteria](#success-criteria)
10. [Migration Strategy](#migration-strategy)

---

## Executive Summary

### Problem Statement
The current CASSIA test suite in `test_code/` is unorganized with scattered test scripts, notebooks, and data files. There's no standardized approach, making it difficult to:
- Understand what each test does
- Reproduce test results
- Add new tests consistently
- Track test coverage
- Automate testing

### Proposed Solution
Create a **modular, standardized testing framework** where:
- Each CASSIA method has a dedicated test folder
- All tests use consistent structure and conventions
- Tests use real sample data from `data/` directory
- Results are organized and version-controlled
- Tests are executable independently or as a suite

### Key Benefits
1. **Clarity**: Each test's purpose and usage is immediately clear
2. **Reproducibility**: Standardized configs ensure consistent results
3. **Maintainability**: Easy to update, extend, or debug individual tests
4. **Automation**: Ready for CI/CD integration
5. **Documentation**: Self-documenting through consistent structure

---

## Current State Analysis

### Existing Test Files

#### In `CASSIA_python/CASSIA/test_code/`:

1. **test.py** (288 lines)
   - Tests `runCASSIA_batch` and `runCASSIA_pipeline`
   - Uses DeepSeek model with retry mechanism
   - Creates synthetic marker data
   - Has timeout monitoring and progress tracking

2. **test_batch_analysis.py** (142 lines)
   - Tests batch analysis with model settings
   - Demonstrates model name resolution
   - Uses sample marker data
   - Tests multiple model configurations

3. **run_batch_simple.py** (132 lines)
   - Simple batch run script
   - Uses real data from `data/unprocessed.csv`
   - Tests gemini and cheap models
   - Demonstrates local imports

4. **quick_test.py**
   - Quick validation test
   - Minimal test case

5. **test_import.py**
   - Import verification test
   - Package structure validation

6. **CASSIA_local_test.ipynb**
   - Interactive testing notebook
   - Exploratory analysis

7. **CASSIA_python_package_test.ipynb**
   - Package functionality testing
   - End-to-end workflows

8. **test_model_settings.ipynb**
   - Model configuration testing
   - Provider testing

9. **test_markers.csv**
   - Synthetic test data
   - 5 cell types with extended marker lists

### Sample Data Available

#### In `CASSIA_python/CASSIA/data/`:

1. **processed.csv**
   - 6 cell types (monocyte, plasma cell, cd8+ T cell, etc.)
   - Clean, well-formatted marker data
   - Suitable for standard testing

2. **unprocessed.csv**
   - Large dataset (3.6MB)
   - Multiple clusters
   - Suitable for stress testing

3. **subcluster_results.csv**
   - Subclustering analysis results
   - Good for subclustering tests

### CASSIA Core Methods (from analysis)

#### Annotation Methods:
1. `runCASSIA` - Single cluster annotation
2. `runCASSIA_batch` - Batch annotation
3. `runCASSIA_pipeline` - Full pipeline
4. `runCASSIA_annotationboost` - Deep-dive iterative analysis
5. `runCASSIA_annotationboost_additional_task` - Boost with custom task

#### Merging Methods:
6. `merge_annotations` - Merge annotations (broad/granular)
7. `merge_annotations_all` - Merge all levels

#### Uncertainty Quantification:
8. `runCASSIA_batch_n_times` - Multiple runs
9. `runCASSIA_similarity_score_batch` - Similarity scoring
10. `runCASSIA_n_times_similarity_score` - Combined UQ

#### Subclustering:
11. `runCASSIA_subclusters` - Subcluster annotation
12. `runCASSIA_n_subcluster` - Multiple subcluster runs
13. `annotate_subclusters` - Annotate existing subclusters

#### Comparison Methods:
14. `compareCelltypes` - Cell type comparison
15. `symphonyCompare` - Multi-agent discussion

#### Utilities:
16. `call_llm` - Direct LLM calls
17. `resolve_model_name` - Model name resolution
18. `get_recommended_model` - Model recommendations
19. `generate_subclustering_report` - Report generation
20. `generate_html_report` - HTML report generation

---

## Proposed Architecture

### Design Principles

1. **Modularity**: Each test is independent and self-contained
2. **Consistency**: All tests follow the same structure and conventions
3. **Configurability**: Tests are driven by JSON configurations
4. **Reusability**: Shared utilities eliminate code duplication
5. **Traceability**: Results are timestamped and versioned
6. **Documentation**: Each test is self-documenting

### Test Categories

#### Category 1: Core Annotation (Priority: High)
- `01_runCASSIA_batch` - Basic batch annotation
- `02_runCASSIA_pipeline` - Full pipeline
- `03_annotation_boost` - Iterative deep-dive

#### Category 2: Post-Processing (Priority: High)
- `04_merge_annotations` - Annotation merging

#### Category 3: Advanced Analysis (Priority: Medium)
- `05_uncertainty_quantification` - UQ analysis
- `06_subclustering` - Subcluster annotation

#### Category 4: Comparison & Debate (Priority: Medium)
- `07_celltype_comparison` - Cell type comparison
- `08_symphony_compare` - Multi-agent symphony

#### Category 5: Utilities (Priority: Low)
- `09_llm_utils` - LLM utilities
- `10_model_settings` - Model configuration
- `11_report_generation` - Report generation

---

## Directory Structure

```
CASSIA_python/CASSIA/
│
├── test/                                    # NEW - Main test directory
│   │
│   ├── README.md                           # Master test documentation
│   ├── run_all_tests.py                    # Master test runner
│   ├── test_results_summary.json           # Aggregated results
│   │
│   ├── shared/                             # Shared utilities
│   │   ├── __init__.py
│   │   ├── test_config.py                 # Config loader
│   │   ├── test_utils.py                  # Logging, timing, validation
│   │   └── sample_data.py                 # Data loading helpers
│   │
│   ├── 01_runCASSIA_batch/
│   │   ├── README.md                      # Test-specific docs
│   │   ├── test_batch.py                  # Test script
│   │   ├── config.json                    # Configuration
│   │   ├── requirements.txt               # Dependencies (if any)
│   │   ├── results/                       # Test outputs
│   │   │   ├── .gitkeep
│   │   │   ├── 20251007_120534_batch_full.csv
│   │   │   ├── 20251007_120534_batch_summary.csv
│   │   │   └── 20251007_120534_test_log.txt
│   │   └── reports/                       # Generated reports
│   │       ├── .gitkeep
│   │       └── 20251007_120534_report.html
│   │
│   ├── 02_runCASSIA_pipeline/
│   │   ├── README.md
│   │   ├── test_pipeline.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 03_annotation_boost/
│   │   ├── README.md
│   │   ├── test_annotation_boost.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 04_merge_annotations/
│   │   ├── README.md
│   │   ├── test_merge.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 05_uncertainty_quantification/
│   │   ├── README.md
│   │   ├── test_uq_batch.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 06_subclustering/
│   │   ├── README.md
│   │   ├── test_subclustering.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 07_celltype_comparison/
│   │   ├── README.md
│   │   ├── test_comparison.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 08_symphony_compare/
│   │   ├── README.md
│   │   ├── test_symphony.py
│   │   ├── config.json
│   │   ├── results/
│   │   │   └── .gitkeep
│   │   └── reports/
│   │       └── .gitkeep
│   │
│   ├── 09_llm_utils/
│   │   ├── README.md
│   │   ├── test_llm_utils.py
│   │   ├── config.json
│   │   └── results/
│   │       └── .gitkeep
│   │
│   ├── 10_model_settings/
│   │   ├── README.md
│   │   ├── test_model_settings.py
│   │   ├── config.json
│   │   └── results/
│   │       └── .gitkeep
│   │
│   └── 11_report_generation/
│       ├── README.md
│       ├── test_reports.py
│       ├── config.json
│       ├── results/
│       │   └── .gitkeep
│       └── reports/
│           └── .gitkeep
│
├── test_code_legacy/                       # ARCHIVED - Old test files
│   ├── README_LEGACY.md
│   ├── test.py
│   ├── test_batch_analysis.py
│   ├── run_batch_simple.py
│   ├── quick_test.py
│   ├── test_import.py
│   ├── test_markers.csv
│   └── *.ipynb
│
└── data/                                    # Sample data (unchanged)
    ├── processed.csv
    ├── unprocessed.csv
    └── subcluster_results.csv
```

---

## Test Module Specifications

### 01_runCASSIA_batch

**Purpose**: Test basic batch annotation functionality
**Method**: `CASSIA.runCASSIA_batch()`
**Expected Runtime**: 3-5 minutes
**Data**: `processed.csv` (6 clusters)

**What it tests**:
- Basic batch annotation
- Multi-cluster processing
- Output format (full + summary CSVs)
- Error handling
- Worker parallelization

**Expected Outputs**:
- `[timestamp]_batch_full.csv` - Full results
- `[timestamp]_batch_summary.csv` - Summary results
- `[timestamp]_test_log.txt` - Test log

**Key Parameters**:
```json
{
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "tissue": "large intestine",
  "species": "human",
  "max_workers": 4,
  "n_genes": 50
}
```

### 02_runCASSIA_pipeline

**Purpose**: Test full end-to-end pipeline
**Method**: `CASSIA.runCASSIA_pipeline()`
**Expected Runtime**: 10-15 minutes
**Data**: `processed.csv` (6 clusters)

**What it tests**:
- Full workflow integration
- Batch annotation
- Quality scoring
- Conditional annotation boost
- Annotation merging
- Report generation

**Expected Outputs**:
- `[timestamp]_pipeline_results.csv`
- `[timestamp]_pipeline_scored.csv`
- `[timestamp]_pipeline_merged.csv`
- `[timestamp]_pipeline_report.html`

**Key Parameters**:
```json
{
  "annotation_model": "google/gemini-2.5-flash-preview",
  "annotation_provider": "openrouter",
  "score_threshold": 97,
  "merge_annotations": true,
  "merge_detail_level": "broad"
}
```

### 03_annotation_boost

**Purpose**: Test iterative marker analysis for ambiguous clusters
**Method**: `CASSIA.runCASSIA_annotationboost()`
**Expected Runtime**: 5-10 minutes
**Data**: `processed.csv` (monocyte cluster)

**What it tests**:
- Iterative hypothesis generation
- Gene checking workflow
- Conversation history management
- Report generation
- Search strategies (breadth/depth)

**Expected Outputs**:
- `[timestamp]_boost_conversation.json`
- `[timestamp]_boost_summary.html`
- `[timestamp]_boost_raw_conversation.txt`

**Key Parameters**:
```json
{
  "cluster_name": "monocyte",
  "conversation_history_mode": "final",
  "search_strategy": "breadth",
  "report_style": "per_iteration",
  "max_iterations": 5
}
```

### 04_merge_annotations

**Purpose**: Test annotation merging at different granularity levels
**Method**: `CASSIA.merge_annotations()`
**Expected Runtime**: 2-3 minutes
**Data**: Results from 01_runCASSIA_batch

**What it tests**:
- Broad-level merging
- Granular-level merging
- Cell type consolidation
- Consensus building

**Expected Outputs**:
- `[timestamp]_merged_broad.csv`
- `[timestamp]_merged_granular.csv`
- `[timestamp]_merged_all.csv`

**Key Parameters**:
```json
{
  "detail_level": "broad",
  "merge_model": "google/gemini-2.5-flash-preview",
  "merge_provider": "openrouter"
}
```

### 05_uncertainty_quantification

**Purpose**: Test uncertainty quantification through multiple runs
**Method**: `CASSIA.runCASSIA_batch_n_times()`
**Expected Runtime**: 15-20 minutes
**Data**: `processed.csv` (subset: 3 clusters)

**What it tests**:
- Multiple stochastic runs (n=5)
- Annotation stability
- Similarity scoring
- Variance metrics
- Confidence intervals

**Expected Outputs**:
- `[timestamp]_uq_run_1.csv` to `run_5.csv`
- `[timestamp]_uq_similarity_matrix.csv`
- `[timestamp]_uq_stability_report.html`

**Key Parameters**:
```json
{
  "n_iterations": 5,
  "temperature": 0.7,
  "similarity_threshold": 0.8
}
```

### 06_subclustering

**Purpose**: Test hierarchical subclustering analysis
**Method**: `CASSIA.runCASSIA_subclusters()`
**Expected Runtime**: 5-8 minutes
**Data**: `subcluster_results.csv`

**What it tests**:
- Subcluster annotation
- Hierarchical analysis
- Report generation
- Multi-level annotations

**Expected Outputs**:
- `[timestamp]_subcluster_annotations.csv`
- `[timestamp]_subcluster_report.html`

**Key Parameters**:
```json
{
  "parent_cluster": "monocyte",
  "n_subclusters": 3,
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter"
}
```

### 07_celltype_comparison

**Purpose**: Test cell type comparison functionality
**Method**: `CASSIA.compareCelltypes()`
**Expected Runtime**: 3-5 minutes
**Data**: Custom marker set

**What it tests**:
- Multi-celltype comparison
- Marker-based scoring
- Comparison matrix generation
- HTML report

**Expected Outputs**:
- `[timestamp]_comparison_results.csv`
- `[timestamp]_comparison_report.html`

**Key Parameters**:
```json
{
  "tissue": "large intestine",
  "celltypes": ["monocyte", "macrophage", "dendritic cell"],
  "species": "human",
  "generate_html_report": true
}
```

### 08_symphony_compare

**Purpose**: Test multi-agent discussion/debate mode
**Method**: `CASSIA.symphonyCompare()`
**Expected Runtime**: 8-12 minutes
**Data**: Ambiguous marker set

**What it tests**:
- Multi-agent consensus
- Discussion rounds
- Debate synthesis
- Confidence evolution

**Expected Outputs**:
- `[timestamp]_symphony_discussion.json`
- `[timestamp]_symphony_consensus.csv`
- `[timestamp]_symphony_report.html`

**Key Parameters**:
```json
{
  "discussion_mode": true,
  "discussion_rounds": 3,
  "model_list": ["google/gemini-2.5-flash-preview"]
}
```

### 09_llm_utils

**Purpose**: Test LLM utility functions
**Method**: `CASSIA.call_llm()`
**Expected Runtime**: 1-2 minutes
**Data**: Simple test prompts

**What it tests**:
- Basic LLM calls
- Provider switching
- Error handling
- Response parsing

**Expected Outputs**:
- `[timestamp]_llm_responses.json`
- `[timestamp]_llm_test_log.txt`

### 10_model_settings

**Purpose**: Test model configuration and resolution
**Method**: `CASSIA.resolve_model_name()`, etc.
**Expected Runtime**: < 1 minute
**Data**: N/A

**What it tests**:
- Model name resolution
- Provider mapping
- Model recommendations
- Configuration validation

**Expected Outputs**:
- `[timestamp]_model_resolution.json`

### 11_report_generation

**Purpose**: Test report generation utilities
**Method**: `CASSIA.generate_html_report()`
**Expected Runtime**: 1-2 minutes
**Data**: Sample annotation results

**What it tests**:
- HTML report generation
- Metric calculation
- Visualization creation
- Report formatting

**Expected Outputs**:
- `[timestamp]_test_report.html`

---

## Shared Infrastructure

### shared/test_config.py

```python
"""
Configuration loader and validator for CASSIA tests.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any

class TestConfig:
    """Test configuration manager."""

    def __init__(self, config_path: Path):
        self.config_path = config_path
        self.config = self._load_config()
        self._validate_config()

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from JSON file."""
        with open(self.config_path, 'r') as f:
            return json.load(f)

    def _validate_config(self):
        """Validate required fields."""
        required = ['test_name', 'model', 'provider']
        for field in required:
            if field not in self.config:
                raise ValueError(f"Missing required field: {field}")

    def get(self, key: str, default=None):
        """Get config value."""
        return self.config.get(key, default)

    def __getitem__(self, key):
        return self.config[key]

def load_config(config_path: Path) -> TestConfig:
    """Load and validate test configuration."""
    return TestConfig(config_path)

def get_data_path(filename: str) -> Path:
    """Get path to data file."""
    data_dir = Path(__file__).parent.parent.parent / "data"
    return data_dir / filename

def get_results_path(test_name: str, filename: str) -> Path:
    """Get path for results file."""
    results_dir = Path(__file__).parent.parent / test_name / "results"
    results_dir.mkdir(exist_ok=True)
    return results_dir / filename
```

### shared/test_utils.py

```python
"""
Common testing utilities for CASSIA tests.
"""

import logging
import time
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from typing import Any, Dict

def setup_logging(test_name: str, level=logging.INFO) -> logging.Logger:
    """Setup logging for test."""
    logger = logging.getLogger(test_name)
    logger.setLevel(level)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger

def get_timestamp() -> str:
    """Get timestamp string for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def log_test_start(logger: logging.Logger, config: Dict) -> float:
    """Log test start and return start time."""
    logger.info("=" * 80)
    logger.info(f"Starting test: {config['test_name']}")
    logger.info(f"Model: {config['model']}")
    logger.info(f"Provider: {config['provider']}")
    logger.info(f"Data: {config.get('data_file', 'N/A')}")
    logger.info("=" * 80)
    return time.time()

def log_test_end(logger: logging.Logger, start_time: float, success: bool = True):
    """Log test completion."""
    elapsed = time.time() - start_time
    status = "SUCCESS" if success else "FAILED"
    logger.info("=" * 80)
    logger.info(f"Test {status}")
    logger.info(f"Elapsed time: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
    logger.info("=" * 80)

def save_results(
    data: Any,
    test_name: str,
    results_dir: Path,
    prefix: str = ""
) -> Path:
    """Save test results with timestamp."""
    timestamp = get_timestamp()

    if isinstance(data, pd.DataFrame):
        filename = f"{timestamp}_{prefix}results.csv"
        output_path = results_dir / filename
        data.to_csv(output_path, index=False)
    elif isinstance(data, dict):
        filename = f"{timestamp}_{prefix}results.json"
        output_path = results_dir / filename
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
    else:
        filename = f"{timestamp}_{prefix}results.txt"
        output_path = results_dir / filename
        with open(output_path, 'w') as f:
            f.write(str(data))

    return output_path

def validate_output(
    data: pd.DataFrame,
    expected_columns: list,
    min_rows: int = 1
) -> bool:
    """Validate output DataFrame."""
    # Check columns
    missing_cols = set(expected_columns) - set(data.columns)
    if missing_cols:
        raise ValueError(f"Missing columns: {missing_cols}")

    # Check row count
    if len(data) < min_rows:
        raise ValueError(f"Expected at least {min_rows} rows, got {len(data)}")

    return True

class Timer:
    """Context manager for timing code blocks."""

    def __init__(self, name: str, logger: logging.Logger = None):
        self.name = name
        self.logger = logger
        self.start_time = None
        self.elapsed = None

    def __enter__(self):
        self.start_time = time.time()
        if self.logger:
            self.logger.info(f"Starting: {self.name}")
        return self

    def __exit__(self, *args):
        self.elapsed = time.time() - self.start_time
        if self.logger:
            self.logger.info(
                f"Completed: {self.name} ({self.elapsed:.2f}s)"
            )
```

### shared/sample_data.py

```python
"""
Sample data loading utilities for CASSIA tests.
"""

import pandas as pd
from pathlib import Path
from typing import Optional, List

class SampleDataLoader:
    """Loader for sample CASSIA data."""

    def __init__(self):
        self.data_dir = Path(__file__).parent.parent.parent / "data"

    def load_processed(self) -> pd.DataFrame:
        """Load processed.csv (6 clean clusters)."""
        return pd.read_csv(self.data_dir / "processed.csv")

    def load_unprocessed(self) -> pd.DataFrame:
        """Load unprocessed.csv (large dataset)."""
        return pd.read_csv(self.data_dir / "unprocessed.csv")

    def load_subcluster_results(self) -> pd.DataFrame:
        """Load subcluster results."""
        return pd.read_csv(self.data_dir / "subcluster_results.csv")

    def load_subset(
        self,
        dataset: str = "processed",
        clusters: Optional[List[str]] = None,
        n_genes: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Load a subset of data.

        Args:
            dataset: 'processed' or 'unprocessed'
            clusters: List of cluster names to include
            n_genes: Max genes per cluster
        """
        if dataset == "processed":
            df = self.load_processed()
        else:
            df = self.load_unprocessed()

        if clusters:
            # Filter by cluster
            # Assumes column name is 'Broad.cell.type' or similar
            cluster_col = df.columns[1]  # Second column usually has cell type
            df = df[df[cluster_col].isin(clusters)]

        if n_genes:
            # Take top N genes per cluster
            marker_col = df.columns[2]  # Third column usually has markers
            # This is simplified - actual implementation depends on format
            pass

        return df

    def get_single_cluster(
        self,
        cluster_name: str,
        dataset: str = "processed"
    ) -> pd.DataFrame:
        """Get data for a single cluster."""
        return self.load_subset(dataset=dataset, clusters=[cluster_name])

# Convenience function
def load_sample_data(dataset: str = "processed") -> pd.DataFrame:
    """Quick load sample data."""
    loader = SampleDataLoader()
    if dataset == "processed":
        return loader.load_processed()
    elif dataset == "unprocessed":
        return loader.load_unprocessed()
    elif dataset == "subcluster":
        return loader.load_subcluster_results()
    else:
        raise ValueError(f"Unknown dataset: {dataset}")
```

---

## Standard Templates

### Template: config.json

```json
{
  "test_name": "test_method_name",
  "description": "Test description goes here",
  "model": "google/gemini-2.5-flash-preview",
  "provider": "openrouter",
  "temperature": 0.7,
  "data_file": "processed.csv",
  "method_params": {
    "tissue": "large intestine",
    "species": "human",
    "max_workers": 4,
    "n_genes": 50,
    "additional_info": "Sample test data"
  },
  "validation": {
    "check_output_format": true,
    "expected_columns": ["cluster", "annotation"],
    "min_rows": 1
  },
  "reporting": {
    "generate_html": true,
    "save_logs": true,
    "verbose": true
  }
}
```

### Template: test_[method].py

```python
#!/usr/bin/env python3
"""
Test: [METHOD_NAME]
Description: [DESCRIPTION]
Expected Runtime: [ESTIMATE]

Usage:
    python test_[method].py

Requirements:
    - OPENROUTER_API_KEY environment variable set
    - Sample data in ../../data/
"""

import os
import sys
from pathlib import Path
import pandas as pd

# Add parent directories to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Import CASSIA
import CASSIA

# Import shared utilities
shared_dir = Path(__file__).parent.parent / "shared"
sys.path.insert(0, str(shared_dir))
from test_config import load_config, get_data_path
from test_utils import (
    setup_logging,
    log_test_start,
    log_test_end,
    save_results,
    validate_output,
    Timer
)
from sample_data import load_sample_data

def main():
    """Main test function."""
    # Setup
    test_dir = Path(__file__).parent
    config = load_config(test_dir / "config.json")
    logger = setup_logging(config['test_name'])

    # Check API key
    api_key = os.getenv('OPENROUTER_API_KEY')
    if not api_key:
        logger.error("OPENROUTER_API_KEY not set in environment")
        logger.error("Set it with: export OPENROUTER_API_KEY='your-key'")
        return 1

    # Set API key
    CASSIA.set_api_key(api_key, provider=config['provider'])

    # Load data
    logger.info(f"Loading data: {config['data_file']}")
    marker_data = load_sample_data(config['data_file'].replace('.csv', ''))
    logger.info(f"Loaded {len(marker_data)} rows")

    # Start test
    start_time = log_test_start(logger, config)

    try:
        # Run CASSIA method
        with Timer("CASSIA method execution", logger):
            result = CASSIA.[METHOD_NAME](
                marker=marker_data,
                model=config['model'],
                provider=config['provider'],
                **config['method_params']
            )

        # Validate results
        if config['validation']['check_output_format']:
            logger.info("Validating output format...")
            validate_output(
                result,
                config['validation']['expected_columns'],
                config['validation'].get('min_rows', 1)
            )
            logger.info("✓ Validation passed")

        # Save results
        logger.info("Saving results...")
        results_dir = test_dir / "results"
        results_dir.mkdir(exist_ok=True)

        output_path = save_results(result, config['test_name'], results_dir)
        logger.info(f"✓ Results saved: {output_path}")

        # Optional: Generate report
        if config['reporting'].get('generate_html', False):
            logger.info("Generating HTML report...")
            # Implementation depends on method
            pass

        log_test_end(logger, start_time, success=True)
        return 0

    except Exception as e:
        logger.error(f"Test failed: {str(e)}", exc_info=True)
        log_test_end(logger, start_time, success=False)
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

### Template: README.md (per test)

```markdown
# Test: [METHOD_NAME]

## Overview

**Method**: `CASSIA.[METHOD_NAME]()`
**Purpose**: [Brief description]
**Category**: [Core Annotation | Post-Processing | Advanced Analysis | etc.]
**Priority**: [High | Medium | Low]

## Description

[Detailed description of what this test does]

## Prerequisites

- Python 3.8+
- CASSIA package installed
- OpenRouter API key set in environment

## Quick Start

```bash
# Set API key
export OPENROUTER_API_KEY='your-key-here'

# Run test
cd test/[XX_test_name]/
python test_[method].py
```

## Configuration

Edit `config.json` to customize:

- **model**: LLM model to use
- **provider**: API provider (openrouter, openai, anthropic)
- **data_file**: Input data file
- **method_params**: Method-specific parameters

## Expected Outputs

### Results Directory (`results/`)

- `[timestamp]_[method]_results.csv` - Main results
- `[timestamp]_test_log.txt` - Execution log

### Reports Directory (`reports/`)

- `[timestamp]_[method]_report.html` - HTML report (if enabled)

## Expected Runtime

**Typical**: [X-Y minutes]
**Maximum**: [Z minutes]

## Validation

The test validates:

- [ ] Output format (columns, dtypes)
- [ ] Minimum row count
- [ ] Required fields present
- [ ] No null values in critical columns
- [ ] [Method-specific checks]

## Sample Output

```csv
cluster,annotation,confidence,sub_types
monocyte,Classical Monocyte,0.95,"CD14+ Monocyte, CD16- Monocyte"
```

## Troubleshooting

### API Key Error
```
Error: OPENROUTER_API_KEY not set
Solution: export OPENROUTER_API_KEY='your-key'
```

### Data Not Found
```
Error: File not found: ../../data/processed.csv
Solution: Ensure you're running from the test directory
```

## Related Tests

- [Link to related tests]

## References

- [CASSIA Documentation](../../../data/CASSIA_Package_Documentation.md)
- [Method Source](../../../[method_file].py)
```

---

## Implementation Roadmap

### Phase 1: Infrastructure Setup (Days 1-2)

**Day 1**:
- [ ] Create `test/` directory structure
- [ ] Create all 11 method folders with subdirectories
- [ ] Add `.gitkeep` files to `results/` and `reports/` folders
- [ ] Create `shared/` directory
- [ ] Write `shared/__init__.py`
- [ ] Write `shared/test_config.py`
- [ ] Write `shared/test_utils.py`
- [ ] Write `shared/sample_data.py`
- [ ] Test shared utilities

**Day 2**:
- [ ] Create master `test/README.md`
- [ ] Create template files (config.json, test_script.py, README.md)
- [ ] Document naming conventions
- [ ] Create `.gitignore` for test outputs
- [ ] Initialize test results tracking

### Phase 2: Core Tests (Days 3-5)

**Day 3**:
- [ ] **01_runCASSIA_batch**
  - [ ] Write `config.json`
  - [ ] Write `test_batch.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

**Day 4**:
- [ ] **02_runCASSIA_pipeline**
  - [ ] Write `config.json`
  - [ ] Write `test_pipeline.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

**Day 5**:
- [ ] **03_annotation_boost**
  - [ ] Write `config.json`
  - [ ] Write `test_annotation_boost.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

### Phase 3: Post-Processing Tests (Day 6)

**Day 6**:
- [ ] **04_merge_annotations**
  - [ ] Write `config.json`
  - [ ] Write `test_merge.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

### Phase 4: Advanced Tests (Days 7-8)

**Day 7**:
- [ ] **05_uncertainty_quantification**
  - [ ] Write `config.json`
  - [ ] Write `test_uq_batch.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

**Day 8**:
- [ ] **06_subclustering**
  - [ ] Write `config.json`
  - [ ] Write `test_subclustering.py`
  - [ ] Write `README.md`
  - [ ] Run and validate
  - [ ] Document results

### Phase 5: Comparison Tests (Day 9)

**Day 9**:
- [ ] **07_celltype_comparison**
  - [ ] Write `config.json`
  - [ ] Write `test_comparison.py`
  - [ ] Write `README.md`
  - [ ] Run and validate

- [ ] **08_symphony_compare**
  - [ ] Write `config.json`
  - [ ] Write `test_symphony.py`
  - [ ] Write `README.md`
  - [ ] Run and validate

### Phase 6: Utility Tests (Day 10)

**Day 10**:
- [ ] **09_llm_utils**
  - [ ] Write `config.json`
  - [ ] Write `test_llm_utils.py`
  - [ ] Write `README.md`
  - [ ] Run and validate

- [ ] **10_model_settings**
  - [ ] Write `config.json`
  - [ ] Write `test_model_settings.py`
  - [ ] Write `README.md`
  - [ ] Run and validate

- [ ] **11_report_generation**
  - [ ] Write `config.json`
  - [ ] Write `test_reports.py`
  - [ ] Write `README.md`
  - [ ] Run and validate

### Phase 7: Integration & Automation (Day 11)

**Day 11**:
- [ ] Create `run_all_tests.py`
- [ ] Add progress tracking
- [ ] Add aggregated reporting
- [ ] Test full suite execution
- [ ] Optimize parallelization
- [ ] Add failure handling

### Phase 8: Migration & Documentation (Day 12)

**Day 12**:
- [ ] Move old files to `test_code_legacy/`
- [ ] Create `README_LEGACY.md`
- [ ] Update main CASSIA documentation
- [ ] Create `TESTING_GUIDE.md` in dev_docs
- [ ] Create `MIGRATION_NOTES.md`
- [ ] Update `CONTRIBUTING.md` if exists
- [ ] Create CI/CD integration guide

---

## Success Criteria

### For Each Individual Test:

1. **Functionality**:
   - [ ] Executes without errors
   - [ ] Produces expected outputs
   - [ ] Validates results correctly

2. **Configuration**:
   - [ ] Uses sample data from `data/`
   - [ ] Uses gemini-2.5-flash model
   - [ ] All parameters in config.json

3. **Documentation**:
   - [ ] Clear README with usage
   - [ ] Expected runtime documented
   - [ ] Sample outputs shown

4. **Organization**:
   - [ ] Results saved with timestamps
   - [ ] Logs captured properly
   - [ ] Follows naming conventions

5. **Reliability**:
   - [ ] Handles errors gracefully
   - [ ] API key validation
   - [ ] Data validation

### For the Overall Framework:

1. **Consistency**:
   - [ ] All tests follow same structure
   - [ ] Shared utilities used consistently
   - [ ] Naming conventions followed

2. **Usability**:
   - [ ] Easy to run individual tests
   - [ ] Easy to add new tests
   - [ ] Clear documentation

3. **Maintainability**:
   - [ ] DRY principle followed
   - [ ] Modular design
   - [ ] Easy to update

4. **Automation**:
   - [ ] Master test runner works
   - [ ] Can run subset of tests
   - [ ] Results aggregated properly

5. **Documentation**:
   - [ ] Master README complete
   - [ ] Each test documented
   - [ ] Migration guide created

---

## Migration Strategy

### Phase 1: Preparation
1. Create new `test/` structure
2. Build and test shared utilities
3. Document migration plan

### Phase 2: Parallel Development
1. Build new tests in `test/`
2. Keep old tests in `test_code/` (working)
3. Validate new tests match old behavior

### Phase 3: Validation
1. Run both old and new tests
2. Compare results
3. Fix discrepancies

### Phase 4: Migration
1. Move old files to `test_code_legacy/`
2. Update all documentation
3. Update CI/CD if applicable

### Phase 5: Cleanup
1. Archive notebooks separately
2. Remove deprecated code
3. Final documentation updates

---

## Appendix

### A. File Naming Conventions

**Test Scripts**:
- Format: `test_[method_name].py`
- Example: `test_batch.py`, `test_pipeline.py`

**Config Files**:
- Format: `config.json`
- One per test directory

**Results Files**:
- Format: `[YYYYMMDD]_[HHMMSS]_[description].[ext]`
- Example: `20251007_143022_batch_full.csv`

**Logs**:
- Format: `[YYYYMMDD]_[HHMMSS]_test_log.txt`

**Reports**:
- Format: `[YYYYMMDD]_[HHMMSS]_report.html`

### B. Environment Variables

Required:
- `OPENROUTER_API_KEY` - OpenRouter API key (primary)

Optional:
- `OPENAI_API_KEY` - For OpenAI tests
- `ANTHROPIC_API_KEY` - For Anthropic tests
- `CASSIA_TEST_DATA_DIR` - Override data directory
- `CASSIA_TEST_VERBOSE` - Enable verbose logging

### C. Git Configuration

Add to `.gitignore`:
```
# Test results
test/*/results/*
!test/*/results/.gitkeep

# Test reports
test/*/reports/*
!test/*/reports/.gitkeep

# Test logs
test/**/test_log*.txt

# Jupyter checkpoints
test/**/.ipynb_checkpoints/
```

### D. Dependencies

All tests use:
- `pandas>=1.3.0`
- `numpy>=1.21.0`
- `CASSIA` (local installation)

Individual tests may require:
- `matplotlib>=3.3.0` (for visualizations)
- `seaborn>=0.11.0` (for plotting)
- `tqdm>=4.60.0` (for progress bars)

---

**End of Plan**

Last Updated: 2025-10-07
Version: 1.0
Author: Claude Code & Elliot Yixuan Xie
