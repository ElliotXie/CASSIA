---
title: 批量处理
---

CASSIA 支持批量处理以同时分析多个聚类。本指南说明如何高效地运行批量分析。

### 准备标记数据

对于 Python，标记数据通常从 CSV 或处理过的数据帧加载。

```python
# 加载您的标记数据
# 可以是 Seurat 的 FindAllMarkers 输出格式
markers = CASSIA.loadmarker(marker_type="unprocessed")
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
    provider = "openrouter"
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

### 输出文件

分析生成以下文件：
1. `{output_name}_full.csv`: 完整的对话历史和推理。
2. `{output_name}_summary.csv`: 精简的结果摘要。
3. `{output_name}_report.html`: 可视化结果的交互式 HTML 报告（单独生成）。

### 获得最佳结果的提示

1. **资源管理**：设置 `max_workers` 时监控系统资源。如果不确定，从保守的数字开始。
2. **上下文优化**：使用 `additional_info` 提供关键的实验细节（例如，“肿瘤样本”，“处理条件”）以指导模型。

