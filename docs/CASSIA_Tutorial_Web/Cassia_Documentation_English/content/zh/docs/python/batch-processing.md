---
title: 批量处理
---

CASSIA 支持批量处理以同时分析多个聚类。本指南说明如何高效地运行批量分析。

### 准备标记数据

您有四种提供标记数据的选项：

1. 创建包含簇和标记基因的 DataFrame 或 CSV 文件
2. 使用 Seurat 的 `FindAllMarkers` 输出（导出为 CSV）
3. 直接使用 Scanpy 的 `rank_genes_groups` 输出
4. 使用 CASSIA 的示例标记数据

```python
import CASSIA
import scanpy as sc
import pandas as pd

# 选项1：加载您自己的标记数据
markers = pd.read_csv("path/to/your/markers.csv")

# 选项2：加载 Seurat FindAllMarkers 输出（导出为 CSV）
markers = pd.read_csv("seurat_markers.csv")

# 选项3：直接使用 Scanpy rank_genes_groups 输出
#（假设您已经有一个计算了 rank_genes_groups 的 AnnData 对象）
markers = sc.get.rank_genes_groups_df(adata, group=None)  # 获取所有组

# 选项4：加载示例标记数据
markers = CASSIA.loadmarker(marker_type="unprocessed")

# 预览数据
print(markers.head())
```

#### 标记数据格式
CASSIA 接受三种格式：

**1. Seurat FindAllMarkers 输出**

Seurat 的 `FindAllMarkers` 函数的标准输出，包含差异表达统计信息：

```
p_val  avg_log2FC  pct.1  pct.2  p_val_adj  cluster  gene
0      3.02        0.973  0.152  0          0        CD79A
0      2.74        0.938  0.125  0          0        MS4A1
0      2.54        0.935  0.138  0          0        CD79B
0      1.89        0.812  0.089  0          1        IL7R
0      1.76        0.756  0.112  0          1        CCR7
```

**2. Scanpy rank_genes_groups 输出（Python 推荐）**

Scanpy 的 `sc.tl.rank_genes_groups()` 函数的输出，通常使用 `sc.get.rank_genes_groups_df()` 导出：

```
group  names   scores  pvals  pvals_adj  logfoldchanges
0      CD79A   28.53   0      0          3.02
0      MS4A1   25.41   0      0          2.74
0      CD79B   24.89   0      0          2.54
1      IL7R    22.15   0      0          1.89
1      CCR7    20.87   0      0          1.76
```

**3. 简化格式**

包含簇 ID 和逗号分隔标记基因的两列 DataFrame：

```
cluster  marker_genes
0        CD79A,MS4A1,CD79B,HLA-DRA,TCL1A
1        IL7R,CCR7,LEF1,TCF7,FHIT,MAL
2        CD8A,CD8B,GZMK,CCL5,NKG7
```

### 运行批量分析

```python
output_name="intestine_detailed"

# 运行批量分析
CASSIA.runCASSIA_batch(
    marker = markers,
    output_name = output_name,
    model = "anthropic/claude-sonnet-4.5",
    tissue = "large intestine",
    species = "human",
    max_workers = 6,  # 匹配聚类数
    n_genes = 50,
    additional_info = None,
    provider = "openrouter",
    reasoning = "medium"  # 可选: "high", "medium", "low" 用于兼容模型
)
```

### 参数详情

- **`marker`**: 标记数据（DataFrame 或文件路径）。数据应包含聚类 ID 和基因名称。
- **`output_name`**: 输出文件的基本名称。
- **`model`**: 要使用的模型（例如，"anthropic/claude-sonnet-4.5"）。
- **`tissue`**: 组织类型（例如，"brain"）。
- **`species`**: 物种（例如，"human"）。
- **`max_workers`**: 并行工作者数。建议使用大约 75-80% 的可用 CPU 核心。
- **`n_genes`**: 每个聚类使用的顶部标记基因数（默认 50）。
- **`additional_info`**: 关于实验的额外背景信息（可选）。
- **`provider`**: API 提供商 ("openrouter", "openai", "anthropic")。
- **`ranking_method`**: 基因排序方法 ("avg_log2FC", "p_val_adj", "pct_diff", "Score")。
- **`ascending`**: 排序方向 (布尔值)。
- **`temperature`**: 模型温度 (0-1)。
- **`celltype_column`**: 聚类 ID 的列名。
- **`gene_column_name`**: 基因符号的列名。
- **`max_retries`**: 失败 API 调用的最大重试次数。
- **`validator_involvement`**: 验证级别 ("v1" 或 "v0")。
- **`reasoning`**: （可选）控制推理深度（"high"、"medium"、"low"）。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。

### 输出文件

分析生成以下文件：
1. `{output_name}_full.csv`: 完整的对话历史和推理。
2. `{output_name}_summary.csv`: 精简的结果摘要。
3. `{output_name}_report.html`: 可视化结果的交互式 HTML 报告（单独生成）。

### 获得最佳结果的提示

1. **资源管理**：设置 `max_workers` 时监控系统资源。如果不确定，从保守的数字开始。
2. **上下文优化**：使用 `additional_info` 提供关键的实验细节（例如，“肿瘤样本”，“处理条件”）以指导模型。

