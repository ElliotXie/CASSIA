##if start with raw data use the code below

# Define a function to process single-cell data
process_data_sc5k <- function(data_sc5k, known_tissue_type = "Immune system", 
                              custom_marker_file = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", 
                              resolution = 0.8, dims = 1:10) {
  
  # Normalize the data using log-normalization
  data_sc5k <- NormalizeData(data_sc5k, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identify highly variable features
  data_sc5k <- FindVariableFeatures(data_sc5k, selection.method = "vst", nfeatures = 2000)
  
  # Scale the data and perform principal component analysis (PCA)
  data_sc5k <- ScaleData(data_sc5k)
  data_sc5k <- RunPCA(data_sc5k)
  
  # Find neighbors and identify clusters
  data_sc5k <- FindNeighbors(data_sc5k, dims = dims)
  data_sc5k <- FindClusters(data_sc5k, resolution = resolution)
  
  # Perform UMAP dimensionality reduction
  data_sc5k <- RunUMAP(data_sc5k, dims = dims)
  
  # Generate a UMAP plot
  DimPlot(data_sc5k, reduction = "umap")
  
  # Run ScType for cell type classification
  data_sc5k <- run_sctype(data_sc5k, known_tissue_type = known_tissue_type, 
                          custom_marker_file = custom_marker_file, 
                          name = "sctype_classification", plot = TRUE)
  
  # Display the ScType classification results
  sctype_table <- data_sc5k$sctype_classification %>% table
  print(sctype_table)
  
  # Return the processed Seurat object
  return(data_sc5k)
}


# Example usage of the process_data_sc5k function on multiple datasets
data_sc5k <- process_data_sc5k(data_sc5k)


##if you want to use existing cluster... or determine the tissue type first

#determine the tissue type

1. redefine a function, tehre is bug in the original version

auto_detect_tissue_type=function(path_to_db_file, seuratObject, scaled, assay = "RNA"){
    
    # get all tissue types in DB
    db_read = openxlsx::read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()
    
    for(tissue in tissues_){ print(paste0("Checking...", tissue));
        
        # prepare gene sets
        gs_list = gene_sets_prepare(path_to_db_file, tissue);

        # check Seurat version
        package_type <- substr(packageVersion("Seurat"), 1, 1)
        data_type <- if (scaled) "scale.data" else "counts"
         obj <- if (package_type == "5") {
          as.matrix(seuratObject[[assay]][data_type])
        } else {
          as.matrix(slot(seuratObject[[assay]], data_type))
        } 
        
        es.max = sctype_score(scRNAseqData = obj, scaled = scaled, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                              marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
        
        cL_resutls = do.call("rbind", lapply(unique(seuratObject@meta.data$clustering_monocle3), function(cl){
            
            es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObject@meta.data[seuratObject@meta.data$clustering_monocle3==cl, ])]), decreasing = !0)
            head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
        }))
        
        dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
        
        # return mean score for tissue
        result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
    }
    
    # order by mean score
    result_ = result_[order(-result_$score),]
    
    # plot 
    barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
            xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")
    
    result_
}


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")


2.gusess the tissue

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(data_merge[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))


# guess a tissue type
tissue_guess <- auto_detect_tissue_type(path_to_db_file = db_, seuratObject = data_merge, scaled = TRUE, assay = "RNA")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used 





##now the analysis for your own clustering,data_merge is my seurat object for demonstration. clustering_monocle3 is the customized clustering column...


# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

tissue <- "Liver" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets

gs_list <- gene_sets_prepare(db_, tissue)

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(data_merge[["RNA"]]$scale.data) else as.matrix(data_merge[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(data_merge@meta.data$clustering_monocle3), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data_merge@meta.data[data_merge@meta.data$clustering_monocle3==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data_merge@meta.data$clustering_monocle3==cl)), 10)
}))


sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)


# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"




data_merge@meta.data$sctype_classification_liver = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data_merge@meta.data$sctype_classification_liver[data_merge@meta.data$clustering_monocle3 == j] = as.character(cl_type$type[1])
}






