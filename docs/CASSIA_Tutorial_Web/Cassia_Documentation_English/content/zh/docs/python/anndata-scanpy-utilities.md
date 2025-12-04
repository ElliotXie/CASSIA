---
title: "AnnData 和 Scanpy 集成工具"
---

CASSIA 提供实用函数来简化与 Scanpy 工作流程和 AnnData 对象的集成。这些函数简化了标记基因提取和注释整合,使单细胞分析流程更加便捷。

## 概述

提供两个用于 Scanpy/AnnData 集成的实用函数:

1. **`enhance_scanpy_markers()`**: 从 Scanpy 的差异表达结果中提取标记基因,并添加 `pct.1` 和 `pct.2` 值(每个基因在簇内/簇外表达的细胞百分比)。

2. **`add_cassia_to_anndata()`**: 自动将 CASSIA 注释结果添加到 AnnData 对象,支持模糊簇匹配和全面的元数据。

## 安装

这些函数需要可选依赖项:

```bash
pip install scanpy anndata
```

在 Python 脚本中导入:

```python
import scanpy as sc
import CASSIA
```

## enhance_scanpy_markers()

### 描述

提取并增强 Scanpy 标记基因,添加表达百分比值(`pct.1` 和 `pct.2`)。此函数专为使用 Scanpy 的 `rank_genes_groups` 进行差异表达分析的工作流程设计。

**重要提示:** 此函数接受 **AnnData 对象**作为输入,而非 DataFrame。它直接从 `adata.uns['rank_genes_groups']` 读取数据(由 `sc.tl.rank_genes_groups()` 存储的结果),并返回一个新的增强 DataFrame。它是 `sc.get.rank_genes_groups_df()` 的**替代品**,而非其输出的后处理器。

### 参数

| 参数 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `adata` | AnnData | *必需* | 带有 `rank_genes_groups` 结果的注释数据矩阵(存储在 `.uns` 中) |
| `cluster_col` | str or None | None | `adata.obs` 中包含簇分配的列名。如果为 None,自动检测(先尝试 'leiden',再尝试 'louvain') |
| `n_genes` | int or None | None | 每个簇包含的顶部基因数量。如果为 None,包含所有基因 |
| `min_expression` | float | 0.0 | 细胞被视为"表达"基因的阈值 |
| `include_stats` | bool | True | 是否包含额外统计信息(logfoldchanges, pvals, scores) |
| `key` | str | "rank_genes_groups" | `adata.uns` 中存储 `rank_genes_groups` 结果的键 |

### 返回值

返回包含以下列的 pandas DataFrame:

| 列名 | 描述 |
|------|------|
| `cluster` | 簇 ID |
| `gene` | 基因名称 |
| `pct.1` | 簇内表达基因的细胞比例 (0.0-1.0) |
| `pct.2` | 簇外表达基因的细胞比例 (0.0-1.0) |
| `avg_log2FC` | 对数折叠变化(如果 `include_stats=True` 且可用) |
| `p_val_adj` | 调整后的 p 值(如果 `include_stats=True` 且可用) |
| `scores` | Scanpy 评分(如果 `include_stats=True` 且可用) |

### 基本用法

```python
import scanpy as sc
import CASSIA

# 加载数据并进行聚类
adata = sc.read_h5ad("your_data.h5ad")
sc.tl.leiden(adata, resolution=0.5)

# 运行差异表达分析
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# 提取增强的标记基因
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)
print(markers.head())
```

### 高级用法

#### 自定义簇列和表达阈值

```python
# 使用自定义簇列和表达阈值
markers = CASSIA.enhance_scanpy_markers(
    adata,
    cluster_col="my_custom_clusters",
    n_genes=100,
    min_expression=0.1,  # 仅统计表达量 > 0.1 的细胞
    include_stats=True
)
```

#### 与 CASSIA 注释集成

```python
# 提取标记基因并运行 CASSIA 注释
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)

# 运行 CASSIA 批量注释
results = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="brain",
    species="mouse",
    output_name="brain_annotation"
)
```

## add_cassia_to_anndata()

### 描述

将 CASSIA 注释结果整合到 AnnData 对象中。此函数自动:
- 匹配 CASSIA 结果和 AnnData 之间的簇 ID(支持模糊匹配)
- 向 `adata.obs` 添加多个注释列
- 在 `adata.uns['CASSIA']` 中存储簇级摘要
- 优雅地处理缺失或不匹配的簇名称

### 参数

| 参数 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `adata` | AnnData | *必需* | 注释数据矩阵 |
| `cassia_results` | str or DataFrame | *必需* | CASSIA 结果 CSV 文件路径或 `runCASSIA_batch()` 的 DataFrame |
| `cluster_col` | str or None | None | `adata.obs` 中包含簇分配的列名。如果为 None,自动检测(先尝试 'leiden',再尝试 'louvain') |
| `cassia_cluster_col` | str | "Cluster ID" | CASSIA 结果中簇 ID 的列名 |
| `prefix` | str | "CASSIA_" | 新列名的前缀 |
| `replace_existing` | bool | False | 是否覆盖现有的 CASSIA 列 |
| `fuzzy_match` | bool | True | 启用簇 ID 对齐的模糊匹配(处理"0"与"cluster_0"等变体) |
| `columns_to_include` | int | 2 | 1 = 仅合并分组, 2 = 所有指标 |
| `inplace` | bool | True | 如果为 True,就地修改 adata。如果为 False,返回副本 |

### 输出

向 `adata.obs` 添加以下列(带指定前缀):

| 列名 | 描述 |
|------|------|
| `CASSIA_general_celltype` | 主要细胞类型预测 |
| `CASSIA_sub_celltype` | 主要亚细胞类型(排名列表中的第一个) |
| `CASSIA_sub_celltype_all` | 亚细胞类型的完整逗号分隔字符串 |
| `CASSIA_sub_celltype_1/2/3` | 拆分的排名候选亚细胞类型 |
| `CASSIA_mixed_celltype` | 可能的混合细胞类型 |
| `CASSIA_score` | 质量/共识评分 (0-100) |
| `CASSIA_merged_grouping_1/2/3` | 层次分组(如果可用) |
| `CASSIA_combined_celltype` | 格式为"一般类型 :: 亚类型" |

同时在 `adata.uns['CASSIA']` 中存储簇级摘要。

### 基本用法

```python
import scanpy as sc
import CASSIA

# 加载 AnnData 和 CASSIA 结果
adata = sc.read_h5ad("your_data.h5ad")
cassia_results = "brain_annotation_FINAL_RESULTS.csv"

# 添加 CASSIA 注释
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="leiden"
)

# 查看注释
print(adata.obs[['leiden', 'CASSIA_general_celltype', 'CASSIA_sub_celltype']].head())

# 查看簇级摘要
print(adata.uns['CASSIA'])
```

### 高级用法

#### 直接使用 runCASSIA_batch 的 DataFrame

```python
# 在一个工作流程中运行 CASSIA 并添加结果
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)

results_df = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="lung",
    species="human",
    output_name="lung_annotation"
)

# 直接从 DataFrame 添加注释(无需读取 CSV)
CASSIA.add_cassia_to_anndata(
    adata,
    results_df,
    cluster_col="leiden"
)
```

#### 自定义前缀和列选择

```python
# 仅添加合并分组,使用自定义前缀
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    cluster_col="seurat_clusters",
    prefix="CellType_",
    columns_to_include=1  # 仅合并分组
)
```

#### 处理不匹配的簇名称

```python
# 模糊匹配自动处理变体
# 例如,adata 中的"0"匹配 CASSIA 结果中的"cluster_0"
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    fuzzy_match=True  # 默认为 True
)
```

## 完整工作流程示例

以下是从聚类到注释的完整工作流程:

```python
import scanpy as sc
import CASSIA

# 设置 API 密钥
CASSIA.set_api_key("your-api-key", provider="openrouter")

# 加载和预处理数据
adata = sc.read_h5ad("pbmc3k.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# 降维和聚类
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 查找标记基因(存储在 adata.uns 中)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# 提取增强的标记基因,包含 pct.1/pct.2
markers = CASSIA.enhance_scanpy_markers(adata, n_genes=50)
print(f"提取了 {len(markers)} 个标记基因")

# 运行 CASSIA 注释
results = CASSIA.runCASSIA_batch(
    marker=markers,
    tissue="blood",
    species="human",
    output_name="pbmc_annotation",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter"
)

# 将注释添加到 AnnData
CASSIA.add_cassia_to_anndata(adata, results)

# 可视化
sc.pl.umap(adata, color=['leiden', 'CASSIA_general_celltype'], ncols=2)

# 访问注释
print(adata.obs[['leiden', 'CASSIA_general_celltype', 'CASSIA_score']].head(20))
print("\n簇级摘要:")
print(adata.uns['CASSIA'])
```

## 疑难解答

### 错误: "rank_genes_groups results not found in adata.uns"

**原因:** 在运行 `sc.tl.rank_genes_groups()` 之前调用了 `enhance_scanpy_markers()`。

**解决方案:** 确保先运行差异表达分析:

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
markers = CASSIA.enhance_scanpy_markers(adata)
```

### 错误: "Could not auto-detect cluster column"

**原因:** `adata.obs` 中既不存在 'leiden' 也不存在 'louvain' 列。

**解决方案:** 明确指定簇列:

```python
markers = CASSIA.enhance_scanpy_markers(adata, cluster_col="my_clusters")
```

### 警告: "Could not find matches for clusters"

**原因:** AnnData 中的簇名称与 CASSIA 结果中的不匹配。

**解决方案:** 检查簇命名并启用模糊匹配(默认):

```python
# 查看簇名称
print("AnnData 簇:", adata.obs['leiden'].unique())
print("CASSIA 簇:", cassia_results['Cluster ID'].unique())

# 默认启用模糊匹配
CASSIA.add_cassia_to_anndata(adata, cassia_results, fuzzy_match=True)
```

### 错误: "CASSIA columns already exist"

**原因:** 当列已存在时尝试添加注释。

**解决方案:** 设置 `replace_existing=True` 以覆盖:

```python
CASSIA.add_cassia_to_anndata(
    adata,
    cassia_results,
    replace_existing=True
)
```

### pct.1 或 pct.2 值缺失(NaN)

**原因:** 在 AnnData 对象中未找到该基因。

**解决方案:** 验证标记结果和 AnnData 之间的基因名称是否匹配:

```python
# 检查基因是否存在
missing_genes = [g for g in markers['gene'].unique() if g not in adata.var_names]
print(f"缺失的基因: {len(missing_genes)}")
```

## 最佳结果提示

1. **使用足够的标记基因**: 每个簇通常使用 50-100 个基因可提供良好的注释质量。
2. **检查表达阈值**: 如果数据有不同的缩放,请调整 `min_expression`。
3. **验证簇匹配**: 始终检查分析和 CASSIA 结果之间的簇 ID 是否对齐。
4. **保存工作**: 添加注释后,保存 AnnData 对象:
   ```python
   adata.write_h5ad("annotated_data.h5ad")
   ```
