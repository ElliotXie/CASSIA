# scCATCH
# install.packages("scCATCH")
library(scCATCH)
library(Seurat)
library(data.table)

# tests based on inputs with normalized count matrix and cell clusters with labels
run_sccatch <- function(annotation_file, expression_file, geo_preprocess = FALSE, tissue, species = 'Human', 
                        cancer_types = 'Normal', gene_info = geneinfo, cell_match = cellmatch, # check cellmatch and geneinfo if needed
                        cell_min_pct=0.1, # scCATCH can have small set of markers requiring a even smaller pct to pass
                        logfc=0.25, use_method = '2') {
  
  # Load annotation, expression
  anno <- fread(annotation_file, data.table = FALSE)
  expr <- fread(expression_file, data.table = FALSE)
  
  # Set row names
  if (geo_preprocess){
    rownames(expr) <- expr[, 1]
    expr <- as.matrix(expr[, -1])
  }
  
  # ensure NCBI symbols (based on 2022 version)
  # however it can delete many genes..
  expr_mat <- rev_gene(data = expr, data_type = 'data', species = species, geneinfo = gene_info)
  
  # Process cell types and subtypes
  cell_types <- anno$Cell_type
  # cell_types_subtype <- paste(anno$Cell_type, anno$Cell_subtype, sep = "_")
  names(cell_types) <- anno[, 1]
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(expr_mat)
  seurat_object@meta.data$orig.ident <- cell_types
  Idents(seurat_object) <- cell_types
  
  # Create scCATCH object
  scCATCH_object <- createscCATCH(data = seurat_object[['RNA']]$counts, cluster = as.character(Idents(seurat_object)))
  
  # Find marker genes and cell types
  scCATCH_object <- findmarkergene(object = scCATCH_object, species = species, cancer = cancer_types, 
                                   marker = cell_match, tissue = tissue, use_method = use_method, 
                                   cell_min_pct = cell_min_pct, logfc = logfc)
  
  scCATCH_object <- findcelltype(scCATCH_object)
  
  return(scCATCH_object)
}

# Example usage:
# result <- process_cancer_data(annotation_file = '../CELLGPT/colon cancer/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz',
#                               expression_file = '../CELLGPT/colon cancer/GSE132465_GEO_processed_CRC_10X_natural_log_TPM_matrix.txt.gz',
#                               geo_preprocess = TRUE, tissue, species = 'Human', 
#                               cancer_types = 'Colon Cancer', gene_info = geneinfo, cell_match = cellmatch, 
#                               cell_min_pct=0.05, logfc=0.25, use_method = '2')




