---
title: "使用 Scanpy 进行聚类和注释"
---

本教程演示了如何使用 Scanpy 对预处理后的 AnnData 对象进行聚类，然后使用 CASSIA 对生成的聚类进行注释。我们假设质量控制步骤已经完成。

## 1. 安装和设置

### 1.1 所需的包

```bash
pip install scanpy leidenalg CASSIA
```

### 1.2 导入包

```python
import scanpy as sc
import pandas as pd
import CASSIA
import os
```

### 1.3 设置 API 密钥

**您只需选择一个提供商。** 推荐使用 OpenRouter，因为它提供多种模型的访问。

```python
os.environ["OPENROUTER_API_KEY"] = "your_openrouter_key"
# os.environ["OPENAI_API_KEY"] = "your_openai_key"
# os.environ["ANTHROPIC_API_KEY"] = "your_anthropic_key"
```

## 2. 加载数据

在本教程中，我们将使用来自 GTEX 项目的乳腺组织数据集作为示例。该数据集为乳腺组织细胞类型提供了全面的参考。

```python
# 加载 GTEX 乳腺数据集 (假设为 .h5ad 格式)
adata = sc.read("gtex_ref.h5ad")
```

探索元数据以了解您的数据集：

```python
print(adata.obs.columns)
```

我们假设数据集在 `Broad_cell_type` 和 `Granular_cell_type` 等列中包含黄金标准细胞类型标签以供参考。

## 3. 降维和聚类

我们将遵循标准的 Scanpy 工作流程，并逐步解释每个步骤。

### 3.1 归一化

首先，对数据进行归一化以消除细胞间测序深度的差异：

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

这将每个细胞缩放到 10,000 总计数，并应用对数变换以减少高表达基因的影响。

### 3.2 特征选择

识别用于聚类的高变基因：

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
```

仅保留高变基因（保留原始数据用于后续可视化）：

```python
adata_raw = adata
adata = adata[:, adata.var.highly_variable]
```

### 3.3 降维

缩放数据并执行 PCA：

```python
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
```

缩放确保每个基因的贡献相等。PCA 将数据降维到捕获最大方差的主成分。

### 3.4 邻居图和 UMAP

构建 k-近邻图并计算 UMAP 用于可视化：

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
sc.tl.umap(adata)
```

邻居图连接相似的细胞，并构成聚类的基础。

### 3.5 聚类

执行 Leiden 聚类以识别细胞群体：

```python
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_res_0.4')
```

`resolution` 参数控制聚类的粒度。较低的值 = 更少、更宽泛的聚类。

### 3.6 可视化聚类

```python
sc.pl.umap(adata, color=['leiden_res_0.4'])
```

![GTEX 乳腺组织聚类的 UMAP 可视化](/images/gtex-umap-clusters.png)

## 4. 寻找标记基因

CASSIA 需要每个聚类的标记基因列表。我们将使用 Wilcoxon 秩和检验来识别差异表达基因。

### 4.1 运行差异表达分析

```python
sc.tl.rank_genes_groups(adata, 'leiden_res_0.4', method='wilcoxon')
```

### 4.2 提取结果

将结果提取到 DataFrame：

```python
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head())
```

### 4.3 格式化为 CASSIA 格式

重命名列以匹配 CASSIA 的预期格式：

```python
markers = markers.rename(columns={
    'names': 'gene',
    'logfoldchanges': 'avg_log2FC',
    'pvals_adj': 'p_val_adj',
    'group': 'cluster'
})
```

筛选正标记（上调基因）：

```python
markers = markers[markers['avg_log2FC'] > 0]
```

## 5. 使用 CASSIA 注释聚类

### 5.1 运行 CASSIA 流程

现在我们可以使用 CASSIA 注释聚类：

```python
results = CASSIA.runCASSIA_pipeline(
    output_file_name = "gtex_breast_annotation",
    tissue = "Breast",
    species = "Human",
    marker = markers,
    max_workers = 6,
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model = "anthropic/claude-sonnet-4.5",
    annotationboost_provider = "openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

要查看每个模型的表现，请访问我们的基准测试网站：[sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

### 5.2 输出文件

流程会创建一个包含以下文件的文件夹：

| 文件 | 描述 |
|------|------|
| `gtex_breast_annotation_summary.csv` | 注释结果摘要 |
| `gtex_breast_annotation_full.csv` | 包含对话历史的完整结果 |
| `gtex_breast_annotation_scored.csv` | 评分后的注释结果 |
| `gtex_breast_annotation_report.html` | 交互式 HTML 报告 |

![CASSIA 注释报告](/images/gtex-breast-annotation-report.webp)

### 5.3 加载结果

```python
cassia_results = pd.read_csv("gtex_breast_annotation_scored.csv")
```

### 5.4 创建注释映射

创建将聚类 ID 映射到细胞类型注释的字典：

```python
annotation_map = dict(zip(
    cassia_results['cluster'].astype(str),
    cassia_results['celltype_1']
))
```

您可以选择不同的粒度级别：
- `celltype_1`：最详细的注释
- `CASSIA_merged_grouping_1`：最宽泛的类别

### 5.5 添加到 AnnData

将注释映射到您的 AnnData 对象：

```python
adata.obs['CASSIA_annotation'] = adata.obs['leiden_res_0.4'].map(annotation_map)
```

验证映射：

```python
print(adata.obs[['leiden_res_0.4', 'CASSIA_annotation']].head())
```

## 6. 可视化注释

将 CASSIA 注释与黄金标准标签进行比较：

```python
sc.pl.umap(adata, color=['Broad_cell_type', 'CASSIA_annotation'], legend_loc='on data')
```

![图 1: 黄金标准聚类](/images/Figure1_GoldStandard.webp)

![图 2: 注释结果](/images/Figure2_CASSIA.webp)

CASSIA 会自动将注释总结为不同的粒度级别，您可以在输出 CSV 中找到这些级别。

## 7. 下一步

完成聚类和注释后：

- 对主要聚类进行子集化并重复分析以获得更精细的分辨率
- 对注释质量低的聚类使用 `CASSIA.runCASSIA_annotationboost()`
- 尝试 `CASSIA.symphonyCompare()` 来区分相似的细胞类型
- 探索 `CASSIA.runCASSIA_subclusters()` 以对特定群体进行更详细的分析
- 实施 `CASSIA.runCASSIA_batch_n_times()` 进行不确定性量化

有关这些高级技术的详细信息，请参阅"扩展分析"教程。
