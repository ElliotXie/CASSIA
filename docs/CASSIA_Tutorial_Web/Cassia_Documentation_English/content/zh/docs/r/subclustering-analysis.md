---
title: 子聚类分析（可选）
---

## 概述

子聚类分析是一种强大的技术，用于更详细地研究特定细胞群体。此功能允许您使用 Cassia 和 Seurat 分析子聚类群体（如 T 细胞或成纤维细胞）。

## 快速开始

```r
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## 输入

### 先决条件
- 包含单细胞数据的 Seurat 对象
- 已安装并加载 Cassia 包
- 对 R 和单细胞分析有基本了解

### 工作流程摘要
1. 在完整数据集上进行初始 Cassia 分析
2. 子簇提取和处理
3. 子簇标记基因识别
4. Cassia 子聚类分析
5. 不确定性评估（可选）

### 准备子聚类

首先，在完整数据集上运行默认的 Cassia 流程，以识别主要细胞群体。然后使用 Seurat 提取和处理目标簇：

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

### 标记基因识别

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

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `marker` | 子聚类的标记基因（数据框或文件路径） |
| `major_cluster_info` | 父簇的描述或上下文（例如，"CD8+ T 细胞"或"与其他细胞类型混合的 cd8 t 细胞"） |
| `output_name` | 输出 CSV 文件的基本名称 |
| `model` | 使用的 LLM 模型 |
| `provider` | API 提供商 |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `temperature` | 0 | 采样温度（0-1） |
| `n_genes` | 50 | 使用的前 N 个标记基因数量 |

### 示例：混合群体

```r
# 对于混合群体
runCASSIA_subclusters(
    marker = marker_sub,
    major_cluster_info = "cd8 t cell mixed with other celltypes",
    output_name = "subclustering_results2",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

### 不确定性评估函数

为获得更可靠的结果，使用多次迭代计算 CS 分数：

**`runCASSIA_n_subcluster()`** - 运行多次注释迭代：

```r
runCASSIA_n_subcluster(
    n = 5,
    marker = marker_sub,
    major_cluster_info = "cd8 t cell",
    base_output_name = "subclustering_results_n",
    model = "anthropic/claude-sonnet-4.5",
    temperature = 0,
    provider = "openrouter",
    max_workers = 5,
    n_genes = 50
)
```

| 参数 | 描述 |
|------|------|
| `n` | 运行的迭代次数 |
| `base_output_name` | 输出文件的基本名称（附加迭代编号） |
| `max_workers` | 并行工作线程数 |

**`runCASSIA_similarity_score_batch()`** - 计算跨迭代的相似性分数：

```r
similarity_scores <- runCASSIA_similarity_score_batch(
    marker = marker_sub,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `file_pattern` | - | 匹配迭代结果文件的 Glob 模式 |
| `main_weight` | 0.5 | 主要细胞类型相似性的权重 |
| `sub_weight` | 0.5 | 亚型相似性的权重 |

## 输出

| 文件 | 描述 |
|------|------|
| `cd8_subcluster_markers.csv` | 每个子簇的标记基因 |
| `{output_name}.csv` | 基本 Cassia 分析结果 |
| `{output_name}.html` | 包含可视化的 HTML 报告 |
| `{output_name}_uncertainty.csv` | 相似性分数（如果执行了不确定性评估） |

> **自动报告生成**：在 CSV 输出的同时会自动生成 HTML 报告，便于可视化子聚类结果。

### 提示和建议
- 在子聚类之前，始终先运行默认的 Cassia 分析
- 根据数据的复杂性调整聚类分辨率
- 处理混合群体时，在 `major_cluster_info` 参数中指定这一点
- 使用不确定性评估获得更稳健的结果
