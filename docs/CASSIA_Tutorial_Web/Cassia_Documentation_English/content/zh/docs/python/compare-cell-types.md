---
title: Symphony Compare (Optional)
---

## Overview

> **注意:** Symphony Compare 需要使用 OpenRouter 作为提供商。这是此模块唯一支持的提供商。

`symphonyCompare` 智能体充当虚拟专家小组来解决模棱两可的细胞类型注释。它协调多个 AI 模型（智能体"交响乐团"）来比较潜在的细胞类型，在讨论轮次中辩论它们的发现，并根据标记基因证据达成共识。

此智能体在运行默认 CASSIA 流程后特别有用，如果您不确定特定聚类的身份。例如，区分浆细胞的不同亚型。

## Quick Start

```python
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells"],
    marker_set = "IGLL5, IGLV6-57, JCHAIN, IGKC, TNFRSF17, IGHG1, MZB1",
    species = "human"
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

## Input

- **标记基因**: 来自 CASSIA 之前结果或您自己分析的逗号分隔的标记基因字符串
- **候选细胞类型**: 要比较的 2-4 种细胞类型的列表
- **组织背景**: 正在分析的组织类型
- **物种**: 样本的物种（默认：人类）

## Parameters

### Required Parameters

| 参数 | 类型 | 描述 |
|-----------|------|-------------|
| `tissue` | 字符串 | 正在分析的组织类型（例如，"large intestine"） |
| `celltypes` | 列表 | 要比较的 2-4 种细胞类型的列表 |
| `marker_set` | 字符串 | 逗号分隔的标记基因字符串 |

### Optional Parameters

| 参数 | 类型 | 默认值 | 描述 |
|-----------|------|---------|-------------|
| `species` | 字符串 | `"human"` | 样本的物种 |
| `model_preset` | 字符串 | `"premium"` | 要使用的模型配置。`"premium"`: 高性能组合 (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4)。`"budget"`: 具成本效益的模型 (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash) |
| `output_basename` | 字符串 | - | 输出文件的基本名称 |
| `enable_discussion` | 布尔值 | `True` | 是否启用模型之间的多轮辩论 |

### Advanced Parameters

| 参数 | 类型 | 默认值 | 描述 |
|-----------|------|---------|-------------|
| `max_discussion_rounds` | 整数 | `2` | 最大讨论轮数 |
| `consensus_threshold` | 浮点数 | `0.8` | 达成共识所需的模型比例（0-1） |

## Output

该函数返回包含共识结果的字典，并生成输出文件。

**返回值:**
- `consensus`: 模型小组达成的共识细胞类型
- `confidence`: 共识的置信度水平（0-1）

**生成的文件:**
- `{output_basename}.csv`: 来自所有模型和轮次的详细比较结果、推理和分数
- `{output_basename}_report.html`: 可视化辩论和共识过程的交互式 HTML 报告
