# Test 09: Subclustering

## Purpose
Tests the subclustering annotation functions for annotating subclusters derived from a major cell type cluster.

## What it Tests
- Subcluster annotation with context from major cluster
- Multi-output extraction (top 2 cell type candidates)
- Key marker identification
- Explanation/reasoning extraction
- HTML report generation

## Use Case
When you have a large cluster identified as a general cell type (e.g., "T cells"), subclustering reveals finer populations. This function annotates each subcluster considering the parent cluster context.

## Functions Tested

### runCASSIA_subclusters()
Main function for subcluster annotation:
- Takes marker data for subclusters
- Major cluster context (e.g., "human lung T cells")
- Returns CSV with annotations for each subcluster

### Output Format
| Column | Description |
|--------|-------------|
| Result ID | Subcluster identifier |
| main_cell_type | Primary cell type prediction |
| sub_cell_type | Secondary cell type prediction |
| key_markers | Most informative markers |
| reason | Biological explanation |

## Test Parameters
- **Major cluster**: human large intestine immune cells
- **Subclusters**: 3 (simulated from marker data)
- **N genes**: 30 markers per subcluster

## Expected Output
```
Subclustering Results:
  Output file: subcluster_results.csv
  Subclusters annotated: 3

  Sample annotations:
    1: Monocyte / Classical monocyte
    2: Plasma cell / IgA plasma cell
    3: CD8+ T cell / Cytotoxic T cell
```

## Running the Test

### Python
```bash
python test_subclustering.py
```

### R
```bash
Rscript test_subclustering.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Subclustering results summary
- `subcluster_results.csv`: Detailed annotations
- `subcluster_results.html`: Visual report (if generated)

## Notes
- Test uses existing marker data as simulated subclusters
- Real usage involves subclusters from single-cell clustering
- Context from major cluster improves annotation specificity
- Temperature of 0.3 allows some flexibility in responses
