# CASSIA

CASSIA (Cell type Annotation using Specialized System with Integrated AI) is a Python package for automated cell type annotation in single-cell RNA sequencing data using large language models.

## Features

- Automated cell type annotation using multiple LLM providers (OpenAI, Anthropic, OpenRouter)
- Support for batch processing of multiple clusters
- Variance analysis for annotation reliability
- Detailed HTML report generation
- Score-based annotation quality assessment
- Support for marker gene validation
- Interactive analysis with step-by-step reasoning

## Installation



```bash

#Prepare dependency
pip install pandas openai requests numpy anthropic

#Install CASSIA
pip install CASSIA
pip install CASSIA_rag  # optional for the RAG agent
```

## Quick Start

```python
import CASSIA

# Set your API keys
CASSIA.set_api_key("your-openai-key", provider="openai")
CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")

# Run single cluster analysis
result = CASSIA.run_cell_type_analysis_wrapper(
    model="gpt-4",
    temperature=0,
    marker_list=["CD3D", "CD4", "IL7R"],
    tissue="PBMC",
    species="human",
    provider="openai"
)

# Run batch analysis
results = CASSIA.run_cell_type_analysis_batchrun(
    marker="path/to/marker_file.csv",
    output_name="results.json",
    model="gpt-4",
    tissue="PBMC",
    species="human"
)
```

## Input Formats

CASSIA accepts marker genes in several formats:
- CSV files with cluster and marker columns
- Pandas DataFrames
- Direct lists of marker genes
- Seurat or Scanpy differential expression results

## Key Functions

### Basic Analysis
- `run_cell_type_analysis_wrapper()`: Single cluster analysis with support for multiple LLM providers
- `run_cell_type_analysis_batchrun()`: Batch analysis of multiple clusters
- `run_analysis_n_times()`: Multiple runs for consensus analysis

### Advanced Features
- `process_cell_type_variance_analysis_batch()`: Analyze annotation variance
- `score_annotation_batch()`: Score annotation quality
- `generate_cell_type_analysis_report_wrapper()`: Generate detailed HTML reports

### Utility Functions
- `set_api_key()`: Set API keys for different providers
- `get_top_markers()`: Extract top markers from differential expression results
- `split_markers()`: Process marker lists into standardized format

## Detailed Usage Examples

### Single Cluster Analysis
```python
# Basic single cluster analysis
result, conversation = CASSIA.run_cell_type_analysis_wrapper(
    model="gpt-4",
    temperature=0,
    marker_list=["CD3D", "CD4", "IL7R", "FOXP3"],
    tissue="PBMC",
    species="human",
    provider="openai"
)

# Print the main cell type identified
print(result["main_cell_type"])
```

### Batch Analysis with Multiple Iterations
```python
# Run batch analysis multiple times for consensus
CASSIA.run_batch_analysis_n_times(
    n=3,  # Number of iterations
    marker="markers.csv",
    output_name="batch_results",
    model="gpt-4",
    tissue="lung",
    species="human",
    max_workers=10
)
```

### Generate Analysis Report
```python
CASSIA.generate_cell_type_analysis_report_wrapper(
    full_result_path="results.csv",
    marker="markers.csv",
    cluster_name="cluster_1",
    major_cluster_info="Human PBMC",
    output_name="analysis_report.html",
    model="gpt-4",
    provider="openai"
)
```

## Contributing

We welcome contributions! Please feel free to submit pull requests or open issues on our GitHub repository.

## License

This project is licensed under the MIT License.

## Support

For support, please open an issue on our [GitHub repository](https://github.com/elliotxe/CASSIA/issues).