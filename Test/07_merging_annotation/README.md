# Test 07: Merging Annotation

## Purpose
Tests the `merge_annotations()` and `merge_annotations_all()` functions for grouping cell type annotations at different detail levels.

## What it Tests
- Single detail level merging (broad, detailed, very_detailed)
- Parallel processing of all detail levels
- Consistency of groupings across detail levels
- Output file generation

## How Merging Works
The merging functions take batch annotation results and group cell types at three levels:

| Level | Column | Description |
|-------|--------|-------------|
| Broad | Merged_Grouping_1 | General categories (e.g., "Myeloid cells", "T cells") |
| Detailed | Merged_Grouping_2 | Intermediate specificity (e.g., "Macrophages", "CD4 T cells") |
| Very Detailed | Merged_Grouping_3 | Normalized specific names (e.g., "Inflammatory macrophages") |

## Prerequisites
- Requires batch annotation results (from Test 02)
- If none exist, the test will run batch annotation first

## Test Scenarios

### Test 1: Broad Merge
Tests `merge_annotations()` with `detail_level="broad"`:
- Groups similar cell types into general categories
- Example: macrophages + dendritic cells -> "Myeloid cells"

### Test 2: Detailed Merge
Tests `merge_annotations()` with `detail_level="detailed"`:
- Provides intermediate-level groupings
- Example: CD4+ naive T cells + CD4+ memory T cells -> "CD4 T cells"

### Test 3: All Levels (Parallel)
Tests `merge_annotations_all()`:
- Processes all three detail levels in parallel
- Outputs a combined CSV with all grouping columns

## Expected Output
- `merge_broad.csv`: Broad groupings
- `merge_detailed.csv`: Detailed groupings
- `merge_all.csv`: All three grouping levels combined

## Running the Test

### Python
```bash
python test_merging_annotation.py
```

### R
```bash
Rscript test_merging_annotation.R
```

## Results
Results are saved to `results/<timestamp>/` containing:
- `test_metadata.json`: Test configuration and status
- `results.json`: Merge results summary
- `merge_*.csv`: Generated merge output files

## Example Output Comparison
For cell type "plasma cell":
| Detail Level | Grouping |
|--------------|----------|
| Broad | B cells / Plasma cells |
| Detailed | Plasma cells |
| Very Detailed | Plasma cells |

For cell type "cd8-positive, alpha-beta t cell":
| Detail Level | Grouping |
|--------------|----------|
| Broad | T cells |
| Detailed | CD8 T cells |
| Very Detailed | CD8+ alpha-beta T cells |

## Notes
- Uses batch_size of 10 for efficient LLM calls
- Parallel processing uses ThreadPoolExecutor with 3 workers
- Lower temperature (0.3) ensures consistent groupings
