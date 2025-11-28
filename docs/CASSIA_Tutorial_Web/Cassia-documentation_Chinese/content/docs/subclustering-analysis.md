---
title: 子聚类分析（可选）
---

子聚类分析是一种强大的技术，用于更详细地研究特定细胞群体。本教程将指导您使用Cassia和Seurat分析子聚类群体（如T细胞或成纤维细胞）的过程。

## 先决条件
- 包含单细胞数据的Seurat对象
- 已安装并加载Cassia包
- 对R和单细胞分析有基本了解

## 工作流程摘要
1. 初始Cassia分析
2. 子簇提取和处理
3. 标记基因识别
4. Cassia子聚类分析
5. 不确定性评估（可选）

## 详细步骤

### 1. 初始分析
首先，在完整数据集上运行默认的Cassia流程，以识别主要细胞群体。

### 2. 子簇处理
使用Seurat提取和处理目标簇：

```r
# 提取目标群体（以CD8+ T细胞为例）
cd8_cells <- subset(large, cell_ontology_class == "cd8-positive, alpha-beta t cell")

# 数据标准化
cd8_cells <- NormalizeData(cd8_cells)

# 识别可变特征
cd8_cells <- FindVariableFeatures(cd8_cells, 
    selection.method = "vst", 
    nfeatures = 2000)

# 缩放数据
all.genes <- rownames(cd8_cells)
cd8_cells <- ScaleData(cd8_cells, features = all.genes)

# 运行PCA
cd8_cells <- RunPCA(cd8_cells, 
    features = VariableFeatures(object = cd8_cells),
    npcs = 30)

# 执行聚类
cd8_cells <- FindNeighbors(cd8_cells, dims = 1:20)
cd8_cells <- FindClusters(cd8_cells, resolution = 0.3)

# 生成UMAP可视化
cd8_cells <- RunUMAP(cd8_cells, dims = 1:20)
```

### 3. 标记基因识别
识别每个子簇的标记基因：

```r
# 寻找标记基因
cd8_markers <- FindAllMarkers(cd8_cells,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25)

# 筛选显著标记
cd8_markers <- cd8_markers %>% filter(p_val_adj < 0.05)

# 保存结果
write.csv(cd8_markers, "cd8_subcluster_markers.csv")
```

### 4. Cassia子聚类分析
对子聚类运行Cassia分析：

```r
# 基本分析
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)

# 对于混合群体
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell mixed with other celltypes",
    output_name = "subclustering_results2",
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)
```

### 5. 不确定性评估（可选）
为获得更可靠的结果，计算CS分数：

```r
# 运行多次迭代
runCASSIA_n_subcluster(
    n = 5,
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results_n",
    model = "anthropic/claude-3.5-sonnet",
    temperature = 0,
    provider = "openrouter",
    max_workers = 5,
    n_genes = 50L
)

# 计算相似性分数
similarity_scores <- runCASSIA_similarity_score_batch(
    marker = marker_sub,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "claude-3-5-sonnet-20241022",
    provider = "anthropic",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

## 提示和建议
- 在子聚类之前，始终先运行默认的Cassia分析
- 根据数据的复杂性调整聚类分辨率
- 处理混合群体时，在`major_cluster_info`参数中指定这一点
- 使用不确定性评估获得更稳健的结果

## 输出文件
分析生成几个输出文件：
- `cd8_subcluster_markers.csv`：每个子簇的标记基因
- `subclustering_results.csv`：基本Cassia分析结果
- `subclustering_uncertainty.csv`：相似性分数（如果执行了不确定性评估）
