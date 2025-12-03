---
title: "使用 Scanpy 进行聚类和注释"
---

本教程演示了如何使用 Scanpy 对预处理后的 AnnData 对象进行聚类，然后使用 CASSIA 对生成的聚类进行注释。我们假设质量控制步骤已经完成。

## 1. 安装和设置

### 1.1 所需的包

```bash
# 安装所需的包
pip install scanpy leidenalg CASSIA
```

### 1.2 导入包

```python
import scanpy as sc
import pandas as pd
import CASSIA
import os

# 设置 API 密钥 (CASSIA 需要)
# 替换为您实际的密钥
os.environ["OPENROUTER_API_KEY"] = "your_openrouter_key"
# os.environ["OPENAI_API_KEY"] = "your_openai_key"
# os.environ["ANTHROPIC_API_KEY"] = "your_anthropic_key"
```

## 2. 探索您的预处理 AnnData 对象

在本教程中，我们将使用来自 GTEX 项目的乳腺组织数据集作为示例。该数据集为乳腺组织细胞类型提供了全面的参考。

```python
# 加载 GTEX 乳腺数据集 (假设为 .h5ad 格式)
# 您可以根据您的数据将其调整为 read_csv 或 read_10x
adata = sc.read("gtex_ref.h5ad")

# 检查元数据
print(adata.obs.columns)
```

我们假设数据集在 `Broad_cell_type` 和 `Granular_cell_type` 等列中包含黄金标准细胞类型标签以供参考。

## 3. 降维和聚类

### 3.1 预处理和聚类

我们将遵循标准的 Scanpy 工作流程：归一化、对数变换、特征选择、缩放、PCA、邻居图构建和聚类。

```python
# 归一化和对数变换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 识别高变基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_raw = adata  # 如果需要，保留原始数据用于可视化
adata = adata[:, adata.var.highly_variable]

# 缩放和 PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 邻居和 UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
sc.tl.umap(adata)

# 聚类 (Leiden)
# 使用分辨率 0.4 以匹配 R 示例
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_res_0.4')

# 可视化
sc.pl.umap(adata, color=['leiden_res_0.4'])
```

![GTEX 乳腺组织聚类的 UMAP 可视化](/images/gtex-umap-clusters.png)

### 3.2 寻找标记基因

CASSIA 需要每个聚类的标记基因列表。我们可以使用 `scanpy.tl.rank_genes_groups` 生成此列表。

```python
# 寻找聚类的标记
sc.tl.rank_genes_groups(adata, 'leiden_res_0.4', method='wilcoxon')

# 将标记提取到 DataFrame 中
markers = sc.get.rank_genes_groups_df(adata, group=None)

# 重命名列以匹配 CASSIA 的预期格式 (类似 Seurat)
# Scanpy 输出: names, scores, logfoldchanges, pvals, pvals_adj
markers = markers.rename(columns={
    'names': 'gene',
    'logfoldchanges': 'avg_log2FC', 
    'pvals_adj': 'p_val_adj', 
    'group': 'cluster'
})

# 筛选正标记 (可选，CASSIA 会处理，但这是个好习惯)
markers = markers[markers['avg_log2FC'] > 0]

print(markers.head())
```

## 4. 使用 CASSIA 注释聚类

### 4.1 基本注释

现在我们有了聚类和标记基因，我们可以使用 CASSIA 进行注释：

```python
# 运行 CASSIA 流程
results = CASSIA.runCASSIA_pipeline(
    output_file_name = "gtex_breast_annotation",
    tissue = "Breast",
    species = "Human",
    marker_path = markers, # 直接传递 DataFrame
    max_workers = 6,  # 匹配数据集中的聚类数
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model="anthropic/claude-sonnet-4.5",
    annotationboost_provider="openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

默认提供商是 OpenRouter，默认模型经过选择以优化注释质量。要查看每个新模型的表现，请访问我们的基准测试网站：[sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

输出文件保存在以组织和物种命名的文件夹中。在文件夹内，您将看到以下文件：

- `gtex_breast_annotation_summary.csv`: 注释结果摘要
- `gtex_breast_annotation_full.csv`: 完整的注释结果，包括完整的对话历史记录和不同级别的合并聚类
- `gtex_breast_annotation_scored.csv`: 评分后的注释结果
- `gtex_breast_annotation_report.html`: 生成的注释结果报告，包含指向所有聚类报告的路由。

`gtex_breast_annotation_scored.csv` 文件：

![CASSIA 注释报告](/images/gtex-breast-annotation-report.png)

### 4.2 将注释整合到 AnnData 对象中

运行 CASSIA 后，您可以将注释整合回您的 AnnData 对象中。我们将把 CSV 结果中的注释映射到聚类观察值。

```python
# 加载结果
cassia_results = pd.read_csv("gtex_breast_annotation_scored.csv")

# 创建映射字典：聚类 ID -> 注释
# 您可以选择 'celltype_1' (最详细) 或 'CASSIA_merged_grouping_1' (最广泛)
annotation_map = dict(zip(cassia_results['cluster'].astype(str), cassia_results['celltype_1']))

# 将注释映射到新的观察列
adata.obs['CASSIA_annotation'] = adata.obs['leiden_res_0.4'].map(annotation_map)

# 检查新注释
print(adata.obs[['leiden_res_0.4', 'CASSIA_annotation']].head())
```

### 4.3 可视化注释

现在我们可以在 UMAP 上可视化 CASSIA 注释。

```python
sc.pl.umap(adata, color=['Broad_cell_type', 'CASSIA_annotation'], legend_loc='on data')
```

![图 1: 黄金标准聚类](/images/Figure1_GoldStandard.png)

![图 2: 注释结果](/images/Figure2_CASSIA.png)

CASSIA 会自动将注释总结为不同粒度级别（从一般到详细），您可以在输出 CSV 中找到这些级别并进行类似的可视化。

## 5. 下一步

完成聚类和注释后：

- 对主要聚类进行子集化并重复步骤
- 对注释质量低的聚类使用 `CASSIA.runCASSIA_annotationboost()`
- 尝试 `CASSIA.symphonyCompare()` 来区分相似的细胞类型
- 探索 `CASSIA.runCASSIA_subclusters()` 以对特定群体进行更详细的分析
- 实施 `CASSIA.runCASSIA_batch_n_times()` 进行不确定性量化

有关这些高级技术的详细信息，请参阅“扩展分析”教程。

