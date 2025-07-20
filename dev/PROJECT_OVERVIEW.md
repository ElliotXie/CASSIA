# CASSIA Project Overview

## Project Structure

```
CASSIA/
├── CASSIA_python/CASSIA/          # Core Python Package
├── CASSIA_R/                      # R Package Wrapper
└── dev/                           # Development Documentation
```

## Core Architecture

### Python Package (`CASSIA_python/CASSIA/`)
- **Main Module**: Single-cell annotation using LLMs
- **Key Files**:
  - `tools_function.py` - Core functions (runCASSIA_batch, runCASSIA_pipeline)
  - `model_settings.py` - Model name resolution system
  - `data/model_settings.json` - Model configuration
  - `llm_utils.py` - LLM API calls
  - `annotation_boost.py` - Advanced annotation
  - `merging_annotation.py` - Annotation merging
  - `subclustering.py` - Subcluster analysis

### R Package (`CASSIA_R/`)
- **Wrapper**: R interface to Python functions via reticulate
- **Key Files**:
  - `R/annotator.R` - Main R functions wrapping Python
  - `R/model_settings.R` - R interface to model settings
  - `inst/python/` - Bundled Python modules
  - `man/` - R documentation

## Model Settings System

### Purpose
Simplifies model selection with memorable names and provider control.

### Key Features
- **Simple Names**: `"gemini"`, `"claude"`, `"gpt4"` → resolve to full model names
- **Quality Shortcuts**: `"best"`, `"cheap"`, `"fast"` per provider
- **Provider Control**: Must specify provider for API key safety
- **Cost Optimization**: Different tiers per provider

### Implementation
```
User Input: model="gemini", provider="openrouter"
    ↓
model_settings.resolve_model_name()
    ↓
Output: "google/gemini-2.5-flash"
```

## API Integration

### Supported Providers
- **OpenAI**: `gpt-4o`, `gpt-4o-mini`, `gpt-3.5-turbo`
- **Anthropic**: `claude-3-5-sonnet-latest`, `claude-3-5-haiku-latest`
- **OpenRouter**: `google/gemini-2.5-flash`, `deepseek/deepseek-chat-v3-0324`, etc.

### Authentication
- Environment variables: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `OPENROUTER_API_KEY`
- R functions: `setOpenaiApiKey()`, `setAnthropicApiKey()`, `setOpenRouterApiKey()`

## Core Workflows

### 1. Batch Analysis
```python
runCASSIA_batch(marker=data, model="gemini", provider="openrouter")
```

### 2. Pipeline Analysis
```python
runCASSIA_pipeline(
    annotation_model="gemini", annotation_provider="openrouter",
    score_model="claude", score_provider="anthropic"
)
```

### 3. Annotation Boost
```python
runCASSIA_annotationboost(marker=data, model="best", provider="anthropic")
```

## Data Formats

### Input
- **Marker Data**: CSV/DataFrame with `gene`, `cluster`, `avg_log2FC`, `p_val_adj`
- **Seurat Objects**: Direct integration via R package

### Output
- **Full Results**: Detailed conversation history and analysis
- **Summary**: Concise predictions and confidence scores
- **Reports**: HTML formatted results

## Development Notes

### Python Package Development
- Local development: Add parent dir to `sys.path`, import directly
- Testing: Use `/test_code/` directory
- Model settings: Managed via `data/model_settings.json`

### R Package Development
- Uses reticulate for Python integration
- Python modules bundled in `inst/python/`
- Environment setup via `setup_cassia_env()`

### Key Dependencies
- **Python**: pandas, numpy, openai, anthropic, requests
- **R**: reticulate, Seurat, dplyr

## Migration Path

### For Existing Users
- **Backward Compatible**: All existing code continues to work
- **Gradual Migration**: Can mix old and new model names
- **Helper Functions**: `resolve_model_name()` shows actual models used

### Example Migration
```r
# Old way (still works)
runCASSIA_batch(model="google/gemini-2.5-flash", provider="openrouter")

# New way
runCASSIA_batch(model="gemini", provider="openrouter")
```

## Testing

### Python Testing
- `/CASSIA_python/CASSIA/test_code/`
- Real data: `/data/unprocessed.csv`
- Model settings tests in Jupyter notebooks

### R Testing
- `/CASSIA_R/test_r/`
- Integration tests with Python backend
- Seurat object compatibility

## Deployment

### Python Package
- Standard pip installation
- Environment setup via conda/virtualenv

### R Package
- CRAN-compatible structure
- Automatic Python environment setup
- Bundled Python dependencies in `inst/python/`

## Future Considerations

- **New LLM Providers**: Add to `model_settings.json`
- **Model Updates**: Update configuration as models evolve
- **API Changes**: Centralized in `llm_utils.py`
- **Performance**: Caching and parallel processing optimizations