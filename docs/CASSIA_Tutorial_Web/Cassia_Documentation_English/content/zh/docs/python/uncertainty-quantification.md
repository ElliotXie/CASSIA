---
title: 不确定性量化（可选）
---

## 概述

CASSIA 中的不确定性量化通过多次分析迭代和相似性评分来评估注释的可靠性。此过程对于以下方面至关重要：
- 识别稳健的细胞类型分配
- 检测混合或模棱两可的聚类
- 量化注释置信度
- 了解预测的可变性

> **成本警告**：使用 LLM 模型运行多次迭代可能会产生大量费用。每次迭代都会进行单独的 API 调用，因此总成本约为单次运行成本的 n 倍。

## 快速开始

### 单聚类分析

```python
from CASSIA import runCASSIA_n_times_similarity_score

result = runCASSIA_n_times_similarity_score(
    tissue="large intestine",
    species="human",
    marker_list=["CD38", "CD138", "JCHAIN", "MZB1", "SDC1"],
    model="openai/gpt-5.1",
    provider="openrouter",
    n=5,
    reasoning="medium"
)

print(f"主要细胞类型: {result['general_celltype_llm']}")
print(f"相似性评分: {result['similarity_score']}")
```

### 批量分析

```python
import CASSIA

# 步骤 1：运行多次迭代
CASSIA.runCASSIA_batch_n_times(
    n=5,
    marker=marker_data,
    output_name="my_annotation",
    model="openai/gpt-5.1",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    reasoning="medium"
)

# 步骤 2：计算相似性评分
CASSIA.runCASSIA_similarity_score_batch(
    marker=marker_data,
    file_pattern="my_annotation_*_full.csv",
    output_name="similarity_results",
    model="openai/gpt-5.1",
    provider="openrouter",
    reasoning="medium"
)
```

## 输入

| 输入 | 描述 | 格式 |
|------|------|------|
| `marker_list` | 单聚类的标记基因 | 基因名称列表 |
| `marker` | 批量处理的标记基因数据 | DataFrame 或文件路径 |
| `tissue` | 组织类型上下文 | 字符串（如 "brain"、"large intestine"） |
| `species` | 物种上下文 | 字符串（如 "human"、"mouse"） |
| `file_pattern` | 匹配迭代结果的模式 | 使用 `*` 通配符的 Glob 模式 |

## 参数

### 单聚类 (`runCASSIA_n_times_similarity_score`)

| 参数 | 必需 | 默认值 | 描述 |
|------|------|--------|------|
| `tissue` | 是 | - | 用于上下文的组织类型 |
| `species` | 是 | - | 用于上下文的物种 |
| `marker_list` | 是 | - | 标记基因列表 |
| `model` | 是 | - | 要使用的 LLM 模型 |
| `provider` | 是 | - | API 提供商（"openrouter"、"openai"、"anthropic"） |
| `n` | 否 | 5 | 分析迭代次数 |
| `temperature` | 否 | 0.3 | LLM 温度（较低 = 更一致） |
| `max_workers` | 否 | 3 | 并行处理工作者数 |
| `main_weight` | 否 | 0.5 | 相似性中主要细胞类型的权重 (0-1) |
| `sub_weight` | 否 | 0.5 | 相似性中亚型的权重 (0-1) |
| `validator_involvement` | 否 | "v1" | 验证器模式（"v0" 严格，"v1" 中等） |
| `additional_info` | 否 | None | 额外上下文字符串 |
| `generate_report` | 否 | True | 生成 HTML 报告 |
| `report_output_path` | 否 | "uq_report.html" | HTML 报告路径 |
| `reasoning` | 否 | None | 推理深度级别（"low"、"medium"、"high"）- 仅适用于 GPT-5 模型 |

### 批量迭代 (`runCASSIA_batch_n_times`)

| 参数 | 必需 | 默认值 | 描述 |
|------|------|--------|------|
| `n` | 是 | - | 分析迭代次数（建议：5） |
| `marker` | 是 | - | 标记基因数据（DataFrame 或路径） |
| `output_name` | 是 | - | 输出文件的基本名称 |
| `model` | 是 | - | 要使用的 LLM 模型 |
| `provider` | 是 | - | API 提供商 |
| `tissue` | 是 | - | 组织类型 |
| `species` | 是 | - | 物种 |
| `max_workers` | 否 | 4 | 总体并行处理限制 |
| `batch_max_workers` | 否 | 2 | 每次迭代的工作者数 |
| `reasoning` | 否 | None | 推理深度级别（"low"、"medium"、"high"）- 仅适用于 GPT-5 模型 |

### 相似性评分 (`runCASSIA_similarity_score_batch`)

| 参数 | 必需 | 默认值 | 描述 |
|------|------|--------|------|
| `marker` | 是 | - | 标记基因数据 |
| `file_pattern` | 是 | - | 匹配迭代结果的模式（如 `"output_*_full.csv"`） |
| `output_name` | 是 | - | 结果的基本名称 |
| `model` | 是 | - | 用于评分的 LLM 模型 |
| `provider` | 是 | - | API 提供商 |
| `max_workers` | 否 | 4 | 并行工作者数 |
| `main_weight` | 否 | 0.5 | 主要细胞类型匹配的重要性 (0-1) |
| `sub_weight` | 否 | 0.5 | 亚型匹配的重要性 (0-1) |
| `generate_report` | 否 | True | 生成 HTML 报告 |
| `report_output_path` | 否 | "uq_batch_report.html" | HTML 报告路径 |
| `reasoning` | 否 | None | 推理深度级别（"low"、"medium"、"high"）- 仅适用于 GPT-5 模型 |

## 输出

### 生成的文件

| 文件 | 描述 |
|------|------|
| `{output_name}_{n}_full.csv` | 每次迭代的结果 |
| `{output_name}_similarity.csv` | 跨迭代的相似性评分 |
| `uq_report.html` / `uq_batch_report.html` | HTML 可视化报告 |

### 返回值（单聚类）

| 键 | 描述 |
|----|------|
| `general_celltype_llm` | 共识主要细胞类型 |
| `sub_celltype_llm` | 共识亚细胞类型 |
| `similarity_score` | 跨迭代的总体相似性 (0-1) |
| `consensus_types` | 出现频率最高的细胞类型 |
| `Possible_mixed_celltypes_llm` | 检测到的混合细胞类型群体 |
| `original_results` | 每次迭代的原始结果 |

### 相似性评分解读

| 分数范围 | 解读 | 操作 |
|----------|------|------|
| > 0.9 | 高一致性 | 稳健的注释 |
| 0.75 - 0.9 | 中等一致性 | 建议审查 |
| < 0.75 | 低一致性 | 使用 [注释增强智能体](annotation-boost.md) 或 [亚群聚类](subclustering-analysis.md) |

### 低分故障排除

1. **检查数据**：检查标记基因质量和聚类异质性
2. **尝试高级智能体**：使用 [注释增强智能体](annotation-boost.md) 或 [亚群聚类](subclustering-analysis.md)
3. **调整参数**：增加迭代次数以获得更可靠的共识
