---
title: "Basic Annotation with Marker Files"
---

This vignette covers the fundamental steps to use CASSIA for cell type annotation when you already have marker gene lists prepared. This is ideal when you've already performed clustering analysis and want to annotate your clusters.

## 1. Installation and Setup

Before starting, make sure you have CASSIA installed and configured. For detailed instructions, see the [**Setting Up CASSIA**](/docs/r/setting-up-cassia) documentation.

```r
library(CASSIA)

# Set up API key (OpenRouter recommended)
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

## 2. Working with Marker Files

CASSIA works with marker gene data from differential expression analysis. It accepts Seurat `FindAllMarkers` output, Scanpy `rank_genes_groups` output, or a simplified format. For detailed format specifications, see the [**Batch Processing**](/docs/r/batch-processing#marker-data-format) documentation.

### 2.1 Sample Data

For this vignette, we'll use the example data included with CASSIA, which contains clusters from a large intestine dataset with six distinct cell populations:
1. Monocyte (Original annotation is inaccurate, this cluster should be Schwann Cell instead. More evidence can be found in the paper)
2. Plasma cells
3. CD8-positive, alpha-beta T cell
4. Transit amplifying cell of large intestine
5. Intestinal enteroendocrine cell
6. Intestinal crypt stem cell

```r
# Load example marker data in both formats
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # Direct Seurat FindAllMarkers output
markers_processed <- loadExampleMarkers(processed = TRUE)     # Processed format

# Preview both data formats
head(markers_unprocessed)
head(markers_processed)
```

## 3. Running Basic Annotation

### 3.1 Fast Mode (All-In-One)

For a quick, comprehensive analysis that performs all steps at once:

```r
# Run the complete CASSIA pipeline in fast mode
fast_results <- runCASSIA_pipeline(
    output_file_name = "CASSIA_Results",
    tissue = "large intestine",
    species = "human",
    marker = markers_unprocessed
)
```

For details about each parameter, please refer to the [**Fast Mode**](/docs/fast-mode) documentation.

### 3.2 Batch Analysis (Faster)

For faster annotation without quality scoring, merging, and annotation boost. This is sufficient for most use cases:

```r
output_name="CASSIA_analysis"

# Run batch analysis with OpenRouter
batch_results <- runCASSIA_batch(
    marker = markers_unprocessed,
    output_name = output_name,
    tissue = "large intestine",
    species = "human",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## 4. Interpreting Results

The HTML report provides a visual summary including:
- Quality scores and reasoning across clusters
- Detailed annotation information for each cluster

For clusters with low scores (<75), consider:
1. Checking if the tissue type was correctly specified
2. Verifying that the species is correct
3. Examining cluster quality through basic QC metrics, doublet detection, and ambient RNA removal
4. Using additional analysis methods like those described in the next section

## 5. Next Steps

After basic annotation, you can explore additional CASSIA capabilities.