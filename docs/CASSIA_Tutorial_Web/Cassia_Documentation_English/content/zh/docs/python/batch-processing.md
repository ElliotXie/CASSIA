---
title: 批量处理
---

## 概述

批量处理允许您一次性注释单细胞数据集中的所有聚类。CASSIA 并行处理多个聚类，为每个聚类生成带有详细推理的细胞类型预测。

---

## 快速开始

```python
CASSIA.runCASSIA_batch(
    marker = markers,
    output_name = "my_annotation",
    model = "anthropic/claude-sonnet-4.5",
    tissue = "brain",
    species = "human",
    provider = "openrouter"
)
```

---

## 输入

### 标记数据格式

CASSIA 接受三种输入格式：

**1. Seurat FindAllMarkers 输出**

Seurat 的 `FindAllMarkers` 函数的标准输出，包含差异表达统计信息：

| p_val | avg_log2FC | pct.1 | pct.2 | p_val_adj | cluster | gene |
|-------|------------|-------|-------|-----------|---------|------|
| 0 | 3.02 | 0.973 | 0.152 | 0 | 0 | CD79A |
| 0 | 2.74 | 0.938 | 0.125 | 0 | 0 | MS4A1 |
| 0 | 2.54 | 0.935 | 0.138 | 0 | 0 | CD79B |
| 0 | 1.89 | 0.812 | 0.089 | 0 | 1 | IL7R |
| 0 | 1.76 | 0.756 | 0.112 | 0 | 1 | CCR7 |

**2. Scanpy rank_genes_groups 输出（Python 推荐）**

Scanpy 的 `sc.tl.rank_genes_groups()` 函数的输出，通常使用 `sc.get.rank_genes_groups_df()` 导出：

| group | names | scores | pvals | pvals_adj | logfoldchanges |
|-------|-------|--------|-------|-----------|----------------|
| 0 | CD79A | 28.53 | 0 | 0 | 3.02 |
| 0 | MS4A1 | 25.41 | 0 | 0 | 2.74 |
| 0 | CD79B | 24.89 | 0 | 0 | 2.54 |
| 1 | IL7R | 22.15 | 0 | 0 | 1.89 |
| 1 | CCR7 | 20.87 | 0 | 0 | 1.76 |

**3. 简化格式**

包含聚类 ID 和逗号分隔标记基因的两列 DataFrame：

| cluster | marker_genes |
|---------|--------------|
| 0 | CD79A,MS4A1,CD79B,HLA-DRA,TCL1A |
| 1 | IL7R,CCR7,LEF1,TCF7,FHIT,MAL |
| 2 | CD8A,CD8B,GZMK,CCL5,NKG7 |

### 加载标记数据

```python
import CASSIA
import scanpy as sc
import pandas as pd

# 选项1：加载 Seurat FindAllMarkers 输出（导出为 CSV）
markers = pd.read_csv("seurat_markers.csv")

# 选项2：直接使用 Scanpy rank_genes_groups 输出（Python 推荐）
#（假设您已经有一个计算了 rank_genes_groups 的 AnnData 对象）
markers = sc.get.rank_genes_groups_df(adata, group=None)  # 获取所有组

# 选项3：加载您自己的简化格式标记数据
markers = pd.read_csv("path/to/your/markers.csv")

# 加载示例标记数据进行测试
markers = CASSIA.loadmarker(marker_type="unprocessed")

# 预览数据
print(markers.head())
```

---

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `marker` | 标记数据（DataFrame 或文件路径） |
| `output_name` | 输出文件的基本名称 |
| `model` | LLM 模型 ID（见下方模型选择） |
| `tissue` | 组织类型（例如 `"brain"`、`"blood"`） |
| `species` | 物种（例如 `"human"`、`"mouse"`） |
| `provider` | API 提供商（`"openrouter"`、`"openai"`、`"anthropic"`）或自定义 base URL |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `max_workers` | 4 | 并行工作者数。建议使用约 75% 的 CPU 核心。 |
| `n_genes` | 50 | 每个聚类使用的顶部标记基因数 |
| `additional_info` | `None` | 额外的实验上下文（见下方） |
| `temperature` | 0 | 输出随机性（0=确定性，1=创造性）。保持为 0 以获得可重复的结果。 |
| `validator_involvement` | `"v1"` | 验证强度：`"v1"`（中等）或 `"v0"`（高，较慢） |
| `reasoning` | `None` | 仅适用于通过 OpenRouter 使用的 GPT-5 系列模型的推理深度（`"low"`、`"medium"`、`"high"`）。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。 |

### 高级参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `ranking_method` | `"avg_log2FC"` | 基因排序：`"avg_log2FC"`、`"p_val_adj"`、`"pct_diff"`、`"Score"` |
| `ascending` | `None` | 排序方向（使用方法默认值） |
| `celltype_column` | `None` | 聚类 ID 的列名 |
| `gene_column_name` | `None` | 基因符号的列名 |
| `max_retries` | 1 | 失败 API 调用的最大重试次数 |

### 参数详情

**模型选择**
- 默认为 `anthropic/claude-sonnet-4.5` 以获得最佳性能
- 使用 `google/gemini-2.5-flash` 进行更快的初步分析
- 详细的模型推荐请参阅 [如何选择模型和提供商](setting-up-cassia.md#如何选择模型和提供商)

**标记基因选择**
- 默认：每个聚类的前 50 个基因
- `ranking_method` 控制标记基因的排序和选择方式：
  - `"avg_log2FC"`（默认）：按平均 log2 倍数变化排序
  - `"p_val_adj"`：按调整后的 p 值排序
  - `"pct_diff"`：按表达百分比差异排序
  - `"Score"`：按自定义分数排序
- 筛选标准：调整后的 p 值 < 0.05，avg_log2FC > 0.25，最小百分比 > 0.1
- 如果通过筛选的基因少于 50 个，则使用所有通过的基因

**附加上下文**
- 使用 `additional_info` 提供实验上下文
- 示例：
  - 处理条件：`"Samples were antibody-treated"`
  - 分析重点：`"Please carefully distinguish between cancer and non-cancer cells"`
- 提示：比较有无附加上下文的结果

---

## 输出

### 生成的文件

| 文件 | 描述 |
|------|------|
| `{output_name}_summary.csv` | 包含细胞类型、标记和元数据的注释结果 |
| `{output_name}_conversations.json` | 用于调试的完整对话历史 |
| `{output_name}_report.html` | 带有可视化的交互式 HTML 报告 |
