# CASSIA Python Codebase Reference

> **Purpose**: Reference document for AI agents and developers to quickly locate Python code, module definitions, and understand CASSIA's architecture.

---

## Quick Navigation

| Mode | File Location | Function |
|------|---------------|----------|
| **Single Annotation** | `CASSIA_python/CASSIA/engine/tools_function.py` | `runCASSIA()` |
| **Batch Mode (Default)** | `CASSIA_python/CASSIA/engine/tools_function.py` | `runCASSIA_batch()` |
| **Pipeline Mode** | `CASSIA_python/CASSIA/pipeline/pipeline.py` | `runCASSIA_pipeline()` |
| **Annotation Boost** | `CASSIA_python/CASSIA/agents/annotation_boost/annotation_boost.py` | `runCASSIA_annotationboost()` |
| **Merging/Grouping** | `CASSIA_python/CASSIA/agents/merging/merging_annotation.py` | `merge_annotations()` |

---

## Directory Structure

```
CASSIA_python/CASSIA/
├── __init__.py                          # Root module - backward-compatible exports
├── engine/                              # Core annotation engine
│   ├── main_function_code.py           # Multi-agent LLM system (3-agent workflow)
│   └── tools_function.py               # High-level API entry points
├── core/                                # Shared utilities
│   ├── llm_utils.py                    # Unified LLM interface (OpenAI, Anthropic, OpenRouter)
│   ├── model_settings.py               # Model configuration & fuzzy name resolution
│   ├── validation.py                   # Input validation for all APIs
│   ├── exceptions.py                   # Custom exception classes
│   ├── marker_utils.py                 # Marker gene processing & ranking
│   ├── anndata_utils.py                # AnnData/scanpy integration
│   ├── api_validation.py               # API key validation
│   ├── progress_tracker.py             # Batch progress tracking
│   ├── utils.py                        # General utilities
│   └── logging_config.py               # Logging configuration
├── pipeline/                            # Pipeline orchestration
│   └── pipeline.py                     # Full workflow: annotation → scoring → boosting
├── agents/                              # Specialized LLM agents
│   ├── annotation_boost/               # Deep analysis for low-scoring clusters
│   │   ├── annotation_boost.py         # Main annotation boost functions
│   │   └── super_annotation_boost.py   # Advanced scanpy-based tools
│   ├── merging/                         # Cluster grouping
│   │   └── merging_annotation.py       # Merge similar annotations
├── reports/                             # HTML report generation
│   ├── generate_reports.py             # Main report functions
│   ├── generate_batch_report.py        # Batch results HTML
├── config/                              # Configuration
│   └── set_api_keys.py                 # API key setup utilities
└── data/                                # Data files
    └── model_settings.json             # Model tier definitions & aliases
```

---

## Module Details

### 1. Default Batch Mode

**File**: `CASSIA_python/CASSIA/engine/tools_function.py`

**Function**: `runCASSIA_batch()`

```python
def runCASSIA_batch(
    marker,                    # DataFrame or CSV path with cluster markers
    tissue,                    # Tissue type (e.g., "brain", "lung")
    species,                   # Species (e.g., "human", "mouse")
    output_name,               # Base filename for outputs
    max_workers=10,            # Parallel threads
    ranking_method='avg_log2FC',  # Marker ranking method
    ascending=False,           # Sort direction
    top_gene=10,               # Top N markers per cluster
    model='balanced',          # LLM model tier
    provider='openai',         # LLM provider
    temperature=0,             # LLM temperature
    use_reference=False,       # Enable reference agent
    ...
)
```

**Outputs**:
- `{output_name}_full.csv` - All columns including conversation history
- `{output_name}_summary.csv` - Key columns only
- `{output_name}_report.html` - Interactive HTML report

---

### 2. Pipeline Mode

**File**: `CASSIA_python/CASSIA/pipeline/pipeline.py`

**Function**: `runCASSIA_pipeline()`

```python
def runCASSIA_pipeline(
    output_name,               # Base output name
    tissue,                    # Tissue type
    species,                   # Species
    marker,                    # Marker data
    score_threshold=50,        # Minimum acceptable score
    merge_annotations=False,   # Enable cluster grouping
    annotationboost_model='best',  # Model for deep analysis
    num_iterations=5,          # Annotation boost iterations
    conversation_history_mode='final',  # 'final', 'full', 'none'
    report_style='per_iteration',  # Report format
    ...
)
```

**Pipeline Steps**:
1. Initial batch annotation (`runCASSIA_batch`)
2. Optional annotation merging (`merge_annotations_all`)
3. Quality scoring (`runCASSIA_score_batch`)
4. Identify low-scoring clusters (< score_threshold)
5. Annotation boost for low-scoring clusters
6. HTML report generation

**Output Structure**:
```
CASSIA_Pipeline_[tissue]_[species]_[timestamp]/
├── 01_annotation_report/     # HTML reports
├── 02_annotation_boost/      # Boosted results
└── 03_csv_files/             # All CSV outputs
```

---

### 3. Annotation Boost

**File**: `CASSIA_python/CASSIA/agents/annotation_boost/annotation_boost.py`

**Function**: `runCASSIA_annotationboost()`

```python
def runCASSIA_annotationboost(
    full_result_path,          # CSV with cluster data
    cluster_name,              # Specific cluster to boost
    tissue,                    # Tissue type
    species,                   # Species
    num_iterations=5,          # Number of refinement iterations
    model='best',              # LLM model
    provider='openai',         # LLM provider
    conversation_history_mode='full',  # History mode
    report_style='per_iteration',  # Report format
    output_dir=None,           # Output directory
    ...
)
```

**Features**:
- Iterative hypothesis generation and refinement
- Dynamic marker checking via LLM
- Conversation history tracking with badge visualization
- Temperature annealing (optional)
- HTML report with step-by-step reasoning

**Super Annotation Boost** (`super_annotation_boost.py`):
- Enhanced with 7 scanpy-based analysis tools
- Automated tool selection
- Gene enrichment analysis (GSEA)
- PCA visualization

---

### 4. Core Annotation Engine (3-Agent System)

**File**: `CASSIA_python/CASSIA/engine/main_function_code.py`

**Function**: `run_cell_type_analysis()`

**Three-Agent Workflow**:
```
1. ANNOTATION AGENT
   - Takes: marker list, tissue, species, additional info
   - Task: Identify cell type step-by-step
   - System prompt: Expert computational biologist persona
   - Output: Detailed cell type identification + reasoning

2. VALIDATION AGENT
   - Takes: Annotation result + marker list + context
   - Task: Validate against marker profile
   - Output: Validation assessment

3. FORMATTING AGENT
   - Takes: Full conversation history
   - Task: Structure output in JSON format
   - Output: JSON with predicted cell types, confidence, reasoning
```

---

### 5. LLM Interface

**File**: `CASSIA_python/CASSIA/core/llm_utils.py`

**Function**: `call_llm()`

```python
def call_llm(
    prompt,                    # User prompt
    provider='openai',         # Provider name
    model='gpt-4o',            # Model name
    temperature=0,             # 0 (deterministic) to 1 (creative)
    system_prompt=None,        # Custom system instructions
    max_tokens=4096,           # Output length limit
    ...
)
```

**Supported Providers**:
- OpenAI (gpt-4o, gpt-5-mini, gpt-5.1)
- Anthropic (claude-sonnet-4-5, claude-opus-4.5)
- OpenRouter (gemini, llama, deepseek)
- Custom OpenAI-compatible endpoints

---

### 6. Model Settings

**File**: `CASSIA_python/CASSIA/core/model_settings.py`

**Class**: `ModelSettings`

**Tier Shortcuts**:
| Tier | Purpose | Example Models |
|------|---------|----------------|
| `"best"` | Highest quality | gpt-5.1, claude-opus-4.5, gemini-2.5-pro |
| `"balanced"` | Good quality/speed | gpt-4o, claude-sonnet-4-5 |
| `"fast"` | Quick responses | gpt-5-mini, gemini-flash |
| `"recommended"` | Internal defaults | Provider-specific |

**Method**: `resolve_model_name(name, provider) → (model, provider)`

---

### 7. Scoring & Evaluation

**File**: `CASSIA_python/CASSIA/evaluation/scoring.py`

**Function**: `runCASSIA_score_batch()`

```python
def runCASSIA_score_batch(
    input_file,                # CSV with annotations
    output_file,               # CSV with scores added
    model='balanced',          # Scoring model
    provider='openai',         # Provider
    max_workers=10,            # Parallel workers
    ...
)
```

**Score Criteria** (0-100):
- Correctness against marker list
- Marker balance (doesn't over-focus on one)
- Captures general picture of cell types
- Considers marker rank importance

---


### 12. Report Generation

**File**: `CASSIA_python/CASSIA/reports/generate_reports.py`

**Key Functions**:
- `generate_batch_html_report_from_data()` - Full batch results
- `generate_analysis_html_report()` - Individual cluster reports
- `runCASSIA_generate_score_report()` - Score summary reports

---

## API Entry Points Summary

| Function | Location | Use Case |
|----------|----------|----------|
| `runCASSIA()` | `engine/tools_function.py` | Single cluster annotation |
| `runCASSIA_batch()` | `engine/tools_function.py` | Standard batch analysis |
| `runCASSIA_pipeline()` | `pipeline/pipeline.py` | Complete workflow |
| `runCASSIA_annotationboost()` | `agents/annotation_boost/annotation_boost.py` | Deep refinement |
| `runCASSIA_n_times()` | `agents/uncertainty/Uncertainty_quantification.py` | Robustness check |
| `runCASSIA_batch_n_times()` | `agents/uncertainty/Uncertainty_quantification.py` | Batch robustness |
| `merge_annotations()` | `agents/merging/merging_annotation.py` | Group annotations |
| `generate_hypothesis()` | `hypothesis/hypothesis_generation.py` | Alternative hypotheses |
| `runCASSIA_score_batch()` | `evaluation/scoring.py` | Quality scoring |
| `call_llm()` | `core/llm_utils.py` | Raw LLM access |
