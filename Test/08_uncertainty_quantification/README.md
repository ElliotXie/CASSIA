# Test 08: Uncertainty Quantification

## Purpose
Tests the uncertainty quantification functions for measuring annotation robustness by running multiple analyses and calculating consensus/similarity scores.

## What it Tests
- Running multiple single-cluster analyses in parallel
- Running batch analyses multiple times for multiple clusters
- Calculating similarity scores across runs
- Consensus detection and cell type unification
- Ontology-based standardization (optional)

## Functions Tested

### 1. runCASSIA_n_times_similarity_score()
Runs CASSIA annotation n times on a **single cluster** and calculates:
- **Similarity score**: How consistent are the results across runs
- **Consensus types**: Most common general and sub cell types
- **LLM-generated consensus**: AI-determined final annotation
- **Possible mixed cell types**: Evidence of mixed populations

**Test Parameters:**
- **Cluster tested**: plasma cell
- **N iterations**: 3 (reduced for testing; typically 5-10 for production)
- **Main weight**: 0.5 (weight for general cell type in similarity)
- **Sub weight**: 0.5 (weight for sub cell type in similarity)

### 2. runCASSIA_batch_n_times()
Runs CASSIA batch annotation n times on **multiple clusters** to measure reproducibility:
- Processes multiple clusters in parallel
- Generates per-iteration output files (CSV + HTML reports)
- Useful for measuring consistency across repeated batch runs

**Test Parameters:**
- **Clusters tested**: monocyte, plasma cell (first 2 from marker file)
- **N iterations**: 2 (reduced for testing)
- **max_workers**: 3
- **batch_max_workers**: 2

## Expected Output

### From runCASSIA_n_times_similarity_score:
```json
{
  "general_celltype_llm": "plasma cell",
  "sub_celltype_llm": "...",
  "similarity_score": 0.85,
  "consensus_types": ["plasma cell", "..."],
  "mixed_celltypes": []
}
```

### From runCASSIA_batch_n_times:
Per iteration, generates 3 files:
- `batch_results_{n}_summary.csv` - Summary with key predictions
- `batch_results_{n}_conversations.json` - Conversation history
- `batch_results_{n}_report.html` - Interactive HTML report

Each CSV contains columns:
| Column | Description |
|--------|-------------|
| Cluster ID | e.g., "monocyte", "plasma cell" |
| Predicted General Cell Type | Main cell type prediction |
| Predicted Detailed Cell Type | Sub-type prediction |
| Possible Mixed Cell Types | Evidence of mixed populations |
| Marker Number | Number of markers used |
| Model, Provider, Tissue, Species | Analysis parameters |

## Interpreting Results

| Similarity Score | Interpretation |
|-----------------|----------------|
| > 0.9 | High confidence - consistent annotation |
| 0.7 - 0.9 | Moderate confidence - some variation |
| 0.5 - 0.7 | Low confidence - significant disagreement |
| < 0.5 | Very low confidence - possible mixed population |

## Running the Test

### Python
```bash
python test_uncertainty_quantification.py
```

### R
```bash
Rscript test_uncertainty_quantification.R
```

## Results
Results are saved to `results/<timestamp>/outputs/` containing:
- `test_metadata.json`: Test configuration and status for both tests
- `results.json`: Combined results from both tests
- `batch_results_1_summary.csv`, `batch_results_1_conversations.json`: First batch iteration
- `batch_results_2_summary.csv`, `batch_results_2_conversations.json`: Second batch iteration
- HTML reports for each iteration

## Notes
- Test uses reduced iterations for speed (n=3 for similarity score, n=2 for batch)
- For production use, n=5-10 provides better statistical power
- Temperature of 0.3 introduces some variation for robustness testing
- Ontology standardization uses EBI OLS API for cell type normalization
- Both tests must pass for overall test to pass
