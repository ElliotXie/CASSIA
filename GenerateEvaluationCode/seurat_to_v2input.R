library(data.table)
library(Seurat)
library(dplyr)
# example on colon cancer
# preprocessing code is from nature method paper
a <- fread('../CELLGPT/colon cancer/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz',data.table=F)
d <- fread('../CELLGPT/colon cancer/GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz',data.table=F)
m <- fread('../CELLGPT/colon cancer/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz',data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
rownames(m) <- m[,1]
m <- as.matrix(m[,-1])
ct <- a$Cell_type
ct_subtype <- paste(a$Cell_type, a$Cell_subtype, sep = "_") # when we need subtype, modification needed
ct[a$Cell_subtype%in%c('CMS1','CMS2','CMS3','CMS4')] <- 'colon cancer cell'
names(ct) <- a[,1]
# saveRDS(ct,file='../CELLGPT/colon cancer/ct.rds')
# saveRDS(d,file='../CELLGPT/colon cancer/norm.rds')
# saveRDS(m,file='../CELLGPT/colon cancer/count.rds')

ct <- readRDS('../CELLGPT/colon cancer/ct.rds')
d <- readRDS('../CELLGPT/colon cancer/norm.rds')
k <- CreateSeuratObject(d)
k@meta.data$orig.ident <- ct
Idents(k) <- ct
g <- FindAllMarkers(k)

sorted_markers <- g %>%
  arrange(cluster, p_val, desc(avg_log2FC))

n_genes <- 10 #number of top markers

top_markers <- sorted_markers %>%
  group_by(cluster) %>%
  slice_head(n = n_genes) %>%
  summarise(Top_Markers = paste(gene, collapse = ", "))

colnames(top_markers) <- c("Broad cell type", "Top 10 Markers")

write.csv(top_markers, "../CELLGPT/colon cancer/modified_markers_colon_cancer_celltype_general.csv", row.names = FALSE)
