---
title: 亚群聚类分析（可选）
---

## 概述

亚群聚类分析是一种更详细地研究特定细胞群体的强大技术。此功能允许您在初始 CASSIA 注释后分析亚群（如 T 细胞或成纤维细胞）。

## 快速开始

```python
CASSIA.runCASSIA_subclusters(
    marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## 输入

### 工作流程摘要
1. 在完整数据集上进行初始 CASSIA 分析
2. 亚群提取和处理（使用 Seurat 或 Scanpy）
3. 亚群标记基因识别
4. CASSIA 亚群聚类分析
5. 不确定性评估（可选）

### 必需输入
- **标记基因**：包含每个亚群标记基因的 DataFrame 或文件路径（来自 Seurat 中的 `FindAllMarkers` 或 Scanpy 中的 `sc.tl.rank_genes_groups` 的输出）

我们建议首先应用默认的 CASSIA。然后，在目标聚类上，应用标准流程（Seurat/Scanpy）进行亚群聚类并获取标记结果。

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `marker` | 亚群的标记基因（DataFrame 或文件路径） |
| `major_cluster_info` | 父聚类或背景的描述（例如，"CD8+ T 细胞"或"与其他细胞类型混合的 cd8 t 细胞"） |
| `output_name` | 输出 CSV 文件的基本名称 |
| `model` | 要使用的 LLM 模型 |
| `provider` | API 提供商 |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `temperature` | 0 | 采样温度 (0-1) |
| `n_genes` | 50 | 要使用的顶部标记基因数 |

### 不确定性评估函数

为了获得更可信的结果，使用多次迭代计算一致性分数 (CS)：

**`runCASSIA_n_subcluster()`** - 运行多次注释迭代：

```python
CASSIA.runCASSIA_n_subcluster(
    n=5,
    marker=subcluster_results,
    major_cluster_info="cd8 t cell",
    base_output_name="subclustering_results_n",
    model="anthropic/claude-sonnet-4.5",
    temperature=0,
    provider="openrouter",
    max_workers=5,
    n_genes=50
)
```

| 参数 | 描述 |
|------|------|
| `n` | 运行的迭代次数 |
| `base_output_name` | 输出文件的基本名称（附加迭代编号） |
| `max_workers` | 并行工作线程数 |

**`runCASSIA_similarity_score_batch()`** - 计算跨迭代的相似性评分：

```python
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
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
| `{output_name}.csv` | 基本 CASSIA 分析结果 |
| `{output_name}.html` | 包含可视化的 HTML 报告 |
| `{output_name}_uncertainty.csv` | 相似性评分（如果进行了不确定性评估） |

> **自动报告生成**：在 CSV 输出的同时会自动生成 HTML 报告，便于可视化亚群聚类结果。
