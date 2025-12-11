# Test 15: CASSIA Pipeline - Input Data

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

**Note:** The test scripts filter this data to use only 2 clusters (monocyte and plasma cell) for faster testing.

**Source:** Human large intestine single-cell RNA-seq data processed with Seurat FindAllMarkers
(Copied from Test 06: Annotation Boost)

## How the Test Uses This Data

1. Each test script loads `findallmarkers_output.csv` from this data folder
2. The raw marker data is filtered to include only 2 clusters: monocyte and plasma cell
3. This filtered data is passed to `runCASSIA_pipeline()` which:
   - Uses it for initial annotation via `runCASSIA_batch`
   - Uses it for annotation boosting when scores fall below threshold (99)
4. With a threshold of 99, annotation boost is guaranteed to trigger for all clusters
5. The annotation boost uses the marker gene statistics (p_val, avg_log2FC, pct.1, pct.2) for iterative refinement

