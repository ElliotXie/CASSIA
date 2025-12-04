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

### 4.2 提取结果

将结果提取到 DataFrame：

```python
markers = sc.get.rank_genes_groups_df(adata, group=None)
print(markers.head())
```

```
  group     names     scores  logfoldchanges          pvals      pvals_adj
0     0     KRT15  48.403599        5.089063   0.000000e+00   0.000000e+00
1     0   TFCP2L1  36.434315        4.719573  1.218868e-290  1.078393e-286
2     0       KIT  35.415207        4.628579  9.961057e-275  5.875364e-271
3     0      NFIB  34.773083        2.135615  6.208985e-265  2.746700e-261
4     0  ANKRD36C  34.635426        3.357123  7.404516e-263  2.620458e-259
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

流程会创建一个包含三个子文件夹的输出文件夹：

- `01_html_reports/` - 每个聚类的交互式 HTML 报告
- `02_cluster_annotations/` - 每个聚类的单独注释结果
- `03_csv_files/` - 汇总 CSV 文件，包括最终结果

### 5.3 加载结果

```python
cassia_results = pd.read_csv("CASSIA_Pipeline_output/gtex_breast_annotation_FINAL_RESULTS.csv")
```

### 5.4 创建注释映射

创建将聚类 ID 映射到细胞类型注释的字典：

```python
annotation_map = dict(zip(
    cassia_results['Cluster ID'].astype(str),
    cassia_results['Predicted General Cell Type']
))
```

您可以选择不同的粒度级别：
- `Predicted General Cell Type`：宽泛的注释

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

CASSIA 会自动将注释总结为不同的粒度级别，您可以在输出 CSV 中找到这些级别。

## 7. 下一步

完成聚类和注释后：

- 对主要聚类进行子集化并重复分析以获得更精细的分辨率
- 对注释质量低的聚类使用 `CASSIA.runCASSIA_annotationboost()`
- 尝试 `CASSIA.symphonyCompare()` 来区分相似的细胞类型
- 探索 `CASSIA.runCASSIA_subclusters()` 以对特定群体进行更详细的分析
- 实施 `CASSIA.runCASSIA_batch_n_times()` 进行不确定性量化

有关这些高级技术的详细信息，请参阅"扩展分析"教程。
