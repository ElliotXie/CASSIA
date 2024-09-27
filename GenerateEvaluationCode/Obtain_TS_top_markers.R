library(Seurat)
library(data.table)
library(dplyr)

tissues = c("bonemarrow", "eye", "pancreas", "skin", "vasculature")

for (i in tissues){
    mat <- fread(paste0("./TS_datasets/",i,"/exprMatrix.tsv.gz"))
    meta <- read.table(paste0("./TS_datasets/",i,"/meta.tsv"), header=T, sep="\t", as.is=T, row.names=1) 

    genes = mat[,1][[1]] 
    genes = gsub(".+[|]", "", genes) 
    mat = data.frame(mat[,-1], row.names=genes) 

    k = CreateSeuratObject(counts = mat, project = i, meta.data=meta)

    k = NormalizeData(object = k)
    k = FindVariableFeatures(object = k)
    k = ScaleData(object = k)
    k = RunPCA(object = k)
    k = FindNeighbors(object = k, dims = 1:30)
    #k = FindClusters(object = k)
    #k = RunUMAP(object = k, dims = 1:30)

    Idents(k) <- k@meta.data$orig.ident <-  meta$cell_ontology_class
    saveRDS(k,file=paste0("./TS_datasets/",i,"/seurat_object"))
    print(paste0(i,"is done saving Seurat Object."))

    g <- FindAllMarkers(k)
    saveRDS(g,file=paste0("./TS_datasets/",i,"/seurat_object_findallmarkers"))
    print(paste0(i,"is done saving Seurat FindAllMarkers."))


    sorted_markers <- g %>%
        arrange(cluster, p_val, desc(avg_log2FC))

    for (n_genes in c(10,30,50)) {
        top_markers <- sorted_markers %>%
            group_by(cluster) %>%
            slice_head(n = n_genes) %>%
            summarise(Top_Markers = paste(gene, collapse = ", "))

        colnames(top_markers) <- c("Broad cell type", "Top Markers")

        write.csv(top_markers, paste0("./TS_datasets/",i,"/TS_top",n_genes,"_markers.csv"), row.names = FALSE)
        print(n_genes)
	}
}
