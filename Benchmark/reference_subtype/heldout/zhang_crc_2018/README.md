# Zhang CRC T-cell held-out benchmark dataset

Source: Zhang et al., "Lineage tracking reveals dynamic relationships of T cells in colorectal cancer", Nature 2018, DOI: 10.1038/s41586-018-0694-x.

Why this is a good held-out test:

- It is a human colorectal cancer T-cell Smart-seq2 dataset, not one of the marker panels used to build the current CASSIA immune T-cell reference files.
- GEO GSE108989 reports 11,138 single T cells from 12 CRC patients and 20 identified T-cell subsets.
- The benchmark input hides the paper cluster labels. Model input should use only neutral `crc_tcell_XX` IDs plus marker genes.
- The ground truth keeps the original Supplementary Table 5 sheet names and expected biological labels for scoring.

Files:

- `supplementary_tables.zip`: Springer/Nature supplementary tables archive downloaded from the paper page.
- `2018-02-02820B-s2/Supplementary Table 5.xlsx`: source signature genes for CRC T-cell clusters.
- `zhang_crc_tcell_marker_inputs.csv`: neutral held-out model input with `cluster_id,marker_genes` only.
- `zhang_crc_tcell_ground_truth.csv`: source sheet, expected label, and scoring terms; do not pass this to the model.

Generation rule:

- Extracted from Supplementary Table 5 sheet order.
- Removed generic ribosomal/housekeeping genes from marker input.
- Forced the sheet marker gene, for example `LEF1`, `LAYN`, `FOXP3`, or `CTLA4`, into the marker list when needed.
