# CASSIA

Collaborative Agent System for Single-cell Interpretable Annotation

## Overview

CASSIA provides R wrappers for a Python-based single-cell annotator. It uses LLMs to perform cell type annotation and analysis on single-cell RNA-seq data.

## Installation

```r
# Install from GitHub
devtools::install_github("ElliotXie/CASSIA", subdir = "CASSIA_R")
```

## Setup

CASSIA requires a Python environment. Set up automatically:

```r
library(CASSIA)
setup_cassia_env()
```

Set your LLM API key:

```r
setLLMApiKey("your-api-key", provider = "anthropic")
```

## Quick Start

```r
library(CASSIA)

# Load example marker data
markers <- loadExampleMarkers()

# Run annotation
result <- runCASSIA(
  marker_list = markers,
  tissue = "PBMC",
  species = "human"
)
```

## Features

- Single and batch cell type annotation
- Multiple LLM provider support (Anthropic, OpenAI, OpenRouter)
- Annotation boosting for iterative refinement
- Subclustering analysis
- Score reports and quality metrics
- Seurat integration

## Documentation

For full documentation, visit: https://github.com/ElliotXie/CASSIA

## License

MIT
