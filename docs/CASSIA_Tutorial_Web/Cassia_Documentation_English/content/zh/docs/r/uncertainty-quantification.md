---
title: 不确定性量化（可选）
---

## 概述

CASSIA 中的不确定性量化通过多次分析迭代和相似性评分来评估注释可靠性。此过程对于以下方面至关重要：
- 识别稳健的细胞类型分配
- 检测混合或模糊的簇
- 量化注释置信度
- 理解预测变异性

> **成本警告**：使用 LLM 模型运行多次迭代可能会产生显著成本。每次迭代都会进行单独的 API 调用，因此总成本将大约是单次运行成本的 n 倍。

## 快速开始

```R
library(CASSIA)

# 步骤 1：运行多次迭代
runCASSIA_batch_n_times(
    n = 5,
    marker = marker_data,
    output_name = "my_annotation",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    tissue = "brain",
    species = "human",
    reasoning = "medium"
)

# 步骤 2：计算相似性分数
runCASSIA_similarity_score_batch(
    marker = marker_data,
    file_pattern = "my_annotation_*_summary.csv",
    output_name = "similarity_results",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    reasoning = "medium"
)
```

## 输入

| 输入 | 描述 | 格式 |
|------|------|------|
| `marker` | 标记基因数据 | 数据框或文件路径 |
| `tissue` | 组织类型上下文 | 字符串（如 "brain"、"large intestine"） |
| `species` | 物种上下文 | 字符串（如 "human"、"mouse"） |
| `file_pattern` | 匹配迭代结果的模式 | 使用 `*` 通配符的 Glob 模式 |

## 参数

### 批量迭代 (`runCASSIA_batch_n_times`)

| 参数 | 必需 | 默认值 | 描述 |
|------|------|--------|------|
| `n` | 是 | - | 分析迭代次数（推荐：5） |
| `marker` | 是 | - | 标记基因数据（数据框或路径） |
| `output_name` | 是 | - | 输出文件的基本名称 |
| `model` | 是 | - | 使用的 LLM 模型 |
| `provider` | 是 | - | API 提供商 |
| `tissue` | 是 | - | 组织类型 |
| `species` | 是 | - | 物种 |
| `max_workers` | 否 | 4 | 整体并行处理限制 |
| `batch_max_workers` | 否 | 2 | 每次迭代的工作进程（max_workers * batch_max_workers 应与核心数匹配） |
| `reasoning` | 否 | NULL | 推理深度级别（"low"、"medium"、"high"）- 仅适用于 GPT-5 模型 |

### 相似性评分 (`runCASSIA_similarity_score_batch`)

| 参数 | 必需 | 默认值 | 描述 |
|------|------|--------|------|
| `marker` | 是 | - | 标记基因数据 |
| `file_pattern` | 是 | - | 匹配迭代结果的模式（如 `"output_*_summary.csv"`） |
| `output_name` | 是 | - | 结果的基本名称 |
| `model` | 是 | - | 用于评分的 LLM 模型 |
| `provider` | 是 | - | API 提供商 |
| `max_workers` | 否 | 4 | 并行工作进程数 |
| `main_weight` | 否 | 0.5 | 主细胞类型匹配的重要性（0-1） |
| `sub_weight` | 否 | 0.5 | 亚类型匹配的重要性（0-1） |
| `generate_report` | 否 | TRUE | 生成 HTML 报告 |
| `report_output_path` | 否 | "uq_batch_report.html" | HTML 报告路径 |
| `reasoning` | 否 | NULL | 推理深度级别（"low"、"medium"、"high"）- 仅适用于 GPT-5 模型 |

## 输出

### 生成的文件

| 文件 | 描述 |
|------|------|
| `{output_name}_{n}_summary.csv` | 每次迭代的结果 |
| `{output_name}_similarity.csv` | 跨迭代的相似性分数 |
| `uq_batch_report.html` | HTML 可视化报告 |

### 相似性分数解读

| 分数范围 | 解读 | 操作 |
|----------|------|------|
| > 0.9 | 高一致性 | 稳健的注释 |
| 0.75 - 0.9 | 中等一致性 | 建议审查 |
| < 0.75 | 低一致性 | 使用 [注释增强智能体](annotation-boost.md) 或 [子聚类](subclustering-analysis.md) |

### 低分故障排除

1. **审查数据**：检查标记基因质量和簇异质性
2. **尝试高级智能体**：使用 [注释增强智能体](annotation-boost.md) 或 [子聚类](subclustering-analysis.md)
3. **调整参数**：增加迭代次数以获得更可靠的共识
