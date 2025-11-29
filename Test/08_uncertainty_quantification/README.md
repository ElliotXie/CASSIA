# Test 08: Uncertainty Quantification

## Purpose
Tests the uncertainty quantification functions for measuring annotation robustness by running multiple analyses and calculating consensus/similarity scores.

## What it Tests
- Running multiple single-cluster analyses in parallel
- Calculating similarity scores across runs
- Consensus detection and cell type unification
- Ontology-based standardization (optional)

## Functions Tested

### runCASSIA_n_times_similarity_score()
Runs CASSIA annotation n times on the same cluster and calculates:
- **Similarity score**: How consistent are the results across runs
- **Consensus types**: Most common general and sub cell types
- **LLM-generated consensus**: AI-determined final annotation
- **Possible mixed cell types**: Evidence of mixed populations

## Test Parameters
- **Cluster tested**: plasma cell
- **N iterations**: 3 (reduced for testing; typically 5-10 for production)
- **Main weight**: 0.5 (weight for general cell type in similarity)
- **Sub weight**: 0.5 (weight for sub cell type in similarity)

## Expected Output
```json
{
  "general_celltype_llm": "plasma cell",
  "sub_celltype_llm": "...",
  "similarity_score": 0.85,
  "consensus_types": ["plasma cell", "..."],
  "mixed_celltypes": []
}
```

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
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Uncertainty quantification results

## Notes
- Test uses reduced iterations (n=3) for speed
- For production use, n=5-10 provides better statistical power
- Temperature of 0.3 introduces some variation for robustness testing
- Ontology standardization uses EBI OLS API for cell type normalization
