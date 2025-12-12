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
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
```

## 2. 加载数据

在本教程中，我们将使用来自 GTEX 项目的乳腺组织数据集作为示例。该数据集为乳腺组织细胞类型提供了全面的参考。

下载数据集：[GTEx_breast_minimal.h5ad](https://drive.google.com/file/d/1HhGX0AD6tfUYzuYc2LToxghTGk8LZeUp/view?usp=sharing)

```python
# 加载 GTEX 乳腺数据集 (假设为 .h5ad 格式)
adata = sc.read("GTEx_breast_minimal.h5ad")
```

探索元数据以了解您的数据集：

```python
print(adata.obs.columns)
```

该数据集在 `Broad cell type` 和 `Granular cell type` 列中包含黄金标准细胞类型标签以供参考。

## 3. 降维和聚类

我们将遵循标准的 Scanpy 工作流程，并逐步解释每个步骤。

### 3.1 归一化

首先，对数据进行归一化以消除细胞间测序深度的差异：

```python
# 保存原始计数数据
adata.layers["counts"] = adata.X.copy()

# 归一化到中位数总计数
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

### 3.2 特征选择

识别用于聚类的高变基因：

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
```

### 3.3 降维

```python
sc.tl.pca(adata)
```

PCA 将数据降维到捕获最大方差的主成分。

### 3.4 邻居图和 UMAP

构建 k-近邻图并计算 UMAP 用于可视化：

```python
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

```python
sc.pl.umap(
    adata,
    color="Broad cell type",
    size=2,
)
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

## 4. 寻找标记基因

CASSIA 需要每个聚类的标记基因列表。我们将使用 Wilcoxon 秩和检验来识别差异表达基因。

### 4.1 运行差异表达分析

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.4", method="wilcoxon",use_raw=False)
```

### 4.2 提取增强的标记基因

使用 CASSIA 的 `enhance_scanpy_markers()` 提取带有表达百分比值(`pct.1` 和 `pct.2`)的标记基因。此函数直接从 `adata.uns['rank_genes_groups']` 读取并添加重要的元数据：

```python
# 提取带有 pct.1/pct.2 值的增强标记基因
markers = CASSIA.enhance_scanpy_markers(adata, cluster_col="leiden_res_0.4", n_genes=50)
print(markers.head())
```

返回的 DataFrame 包含：
- `pct.1`: 簇内表达每个基因的细胞百分比
- `pct.2`: 簇外表达每个基因的细胞百分比
- 标准统计信息: `avg_log2FC`, `p_val_adj`, `scores`

这些百分比值帮助 CASSIA 更好地理解标记基因的特异性，对于注释增强尤其有用。

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
    annotation_model = "openai/gpt-5.1",
    annotation_provider = "openrouter",
    score_model = "anthropic/claude-sonnet-4.5",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model = "openai/gpt-5.1",
    annotationboost_provider = "openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

要查看每个模型的表现，请访问我们的基准测试网站：[sc-llm-benchmark.com/methods/cassia](https://sc-llm-benchmark.com/methods/cassia)

### 5.2 输出文件

流程会创建一个名为 `CASSIA_Pipeline_{tissue}_{species}_{timestamp}/` 的输出文件夹，包含三个子文件夹：

- `01_annotation_report/` - 分析的交互式 HTML 报告
- `02_annotation_boost/` - 低评分聚类的注释增强结果
- `03_csv_files/` - 汇总 CSV 文件，包括最终结果

### 5.3 加载结果

```python
# 将文件夹名称替换为实际的输出文件夹（包含时间戳）
cassia_results = pd.read_csv("CASSIA_Pipeline_Breast_Human_XXXXXX/03_csv_files/gtex_breast_annotation_FINAL_RESULTS.csv")
```

### 5.4 将注释添加到 AnnData

使用 CASSIA 的 `add_cassia_to_anndata()` 自动将所有注释结果整合到 AnnData 对象中：

```python
# 自动进行簇匹配并添加 CASSIA 注释
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="leiden_res_0.4",
    prefix="CASSIA_"
)
```

此函数自动：
- 匹配 CASSIA 结果和 AnnData 之间的簇 ID（支持模糊匹配）
- 向 `adata.obs` 添加多个注释列
- 在 `adata.uns['CASSIA']` 中存储簇级摘要

### 5.5 查看注释

该函数向 `adata.obs` 添加了多个有用的列：

```python
# 查看添加的列
cassia_cols = [col for col in adata.obs.columns if col.startswith('CASSIA_')]
print("添加的列:", cassia_cols)

# 查看样本注释
print(adata.obs[['leiden_res_0.4', 'CASSIA_general_celltype', 'CASSIA_sub_celltype', 'CASSIA_score']].head(10))

# 查看簇级摘要
print("\n簇级摘要:")
print(adata.uns['CASSIA'])
```

可用的注释列包括：
- `CASSIA_general_celltype`: 宽泛的细胞类型注释
- `CASSIA_sub_celltype`: 详细的亚细胞类型注释
- `CASSIA_score`: 置信度评分 (0-100)
- `CASSIA_combined_celltype`: 格式为"一般类型 :: 亚类型"
- 合并分组和替代预测的其他列

## 6. 可视化注释

将 CASSIA 注释与黄金标准标签进行比较：

```python
# 将 CASSIA 注释与原始聚类一起可视化
sc.pl.umap(adata, color=['leiden_res_0.4', 'CASSIA_general_celltype'], ncols=2, size=2)

# 与黄金标准标签进行比较（如果可用）
sc.pl.umap(adata, color=['CASSIA_general_celltype', 'Broad cell type'], ncols=2, size=2)
```

您还可以创建汇总统计信息：

```python
# 统计每种细胞类型的细胞数量
celltype_counts = adata.obs['CASSIA_general_celltype'].value_counts()
print(celltype_counts)

# 比较聚类分配与 CASSIA 注释
comparison = adata.obs.groupby(['leiden_res_0.4', 'CASSIA_general_celltype']).size().unstack(fill_value=0)
print(comparison)
```

CASSIA 会自动将注释总结为不同的粒度级别，您可以在输出 CSV 和 `adata.obs` 中的各种 `CASSIA_*` 列中找到这些级别。

## 7. 下一步

完成聚类和注释后：

- 对主要聚类进行子集化并重复分析以获得更精细的分辨率
- 对注释质量低的聚类使用 `CASSIA.runCASSIA_annotationboost()`
- 尝试 `CASSIA.symphonyCompare()` 来区分相似的细胞类型
- 探索 `CASSIA.runCASSIA_subclusters()` 以对特定群体进行更详细的分析
- 实施 `CASSIA.runCASSIA_batch_n_times()` 进行不确定性量化

有关这些高级技术的详细信息，请参阅"扩展分析"教程。
