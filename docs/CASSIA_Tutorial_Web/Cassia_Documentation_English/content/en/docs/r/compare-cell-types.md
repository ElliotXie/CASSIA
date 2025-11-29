---
title: Symphony Compare (Advanced)
---


`symphonyCompare` is an advanced module that orchestrates multiple AI models to compare cell types with automatic consensus building. It conducts a comprehensive cell type comparison using multiple AI models in parallel, automatically triggering discussion rounds when models disagree on the best matching cell type. Think of it as a virtual panel of expert biologists debating and reaching consensus.

### Usage Example

```r
# Basic usage - let Symphony Compare handle everything
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

# Access the results
cat("Consensus:", results$consensus, "\n")
cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
```

### Function Parameters

```r
symphonyCompare(
    tissue,             # Tissue type
    celltypes,          # Vector of 2-4 cell types to compare
    marker_set,         # Vector or string of marker genes
    species = "human",  # Species
    model_preset = "symphony", # Preset model configuration
    enable_discussion = TRUE,  # Enable automatic discussion rounds
    max_discussion_rounds = 2, # Maximum discussion rounds
    consensus_threshold = 0.8, # Threshold for consensus (0-1)
    generate_report = TRUE,    # Generate HTML report
    verbose = TRUE             # Print progress messages
)
```

### Parameter Details

- **`tissue`** (character string): The tissue type being analyzed (e.g., "blood", "brain", "liver").
- **`celltypes`** (character vector): A vector of 2-4 cell types to compare.
- **`marker_set`** (character vector or string): A list or string of marker genes.
- **`species`** (character string): The species of the sample (default: "human").
- **`model_preset`** (character string): Pre-configured model ensemble to use.
  - `"symphony"`: High-performance ensemble (Claude, GPT-4, Gemini Pro)
  - `"quartet"`: Balanced 4-model ensemble
  - `"budget"`: Cost-effective models
  - `"custom"`: Use custom_models list
- **`enable_discussion`** (logical): If `TRUE`, models will "discuss" and reconsider their votes if initial consensus is not reached.
- **`max_discussion_rounds`** (integer): Maximum number of discussion rounds allowed (default: 2).
- **`consensus_threshold`** (numeric): The confidence threshold required to declare consensus (0-1, default: 0.8).
- **`generate_report`** (logical): Whether to generate an HTML report of the analysis (default: `TRUE`).
- **`verbose`** (logical): Whether to print progress messages to the console (default: `TRUE`).

### Output Format

The function returns a list containing:
- **`results`**: List of all model responses and scores.
- **`consensus`**: The consensus cell type (if reached).
- **`confidence`**: Confidence level of the consensus (0-1).
- **`csv_file`**: Path to the generated CSV file with detailed results.
- **`html_file`**: Path to the generated interactive HTML report (if enabled).
- **`summary`**: Summary statistics of the comparison.
