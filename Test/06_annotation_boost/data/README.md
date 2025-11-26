# Test 06: Annotation Boost - Input Data

## Files

### findallmarkers_output.csv
Seurat FindAllMarkers output file containing differential expression results for each cluster.

**Columns:**
- `p_val`: P-value from differential expression test
- `avg_log2FC`: Average log2 fold change
- `pct.1`: Percentage of cells in cluster expressing gene
- `pct.2`: Percentage of cells outside cluster expressing gene
- `p_val_adj`: Adjusted p-value (Bonferroni corrected)
- `cluster`: Cell type cluster name
- `gene`: Gene symbol

**Cell Types Included:**
1. monocyte
2. plasma cell
3. cd8-positive, alpha-beta t cell
4. transit amplifying cell of large intestine
5. intestinal enteroendocrine cell
6. intestinal crypt stem cell

**Source:** Human large intestine single-cell RNA-seq data processed with Seurat FindAllMarkers

## How the Test Uses This Data

1. The test loads `findallmarkers_output.csv` which contains the raw marker gene results
2. It either uses existing batch annotation results (from test 02) or runs a new batch annotation
3. Then it runs `runCASSIA_annotationboost()` on a specific cluster (default: "plasma cell")
4. The annotation boost uses the marker gene statistics to perform iterative refinement using LLM analysis
