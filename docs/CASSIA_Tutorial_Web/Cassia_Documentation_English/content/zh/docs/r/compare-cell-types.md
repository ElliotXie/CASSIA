---
title: Symphony Compare (可选)
---

## Overview

> **注意:** Symphony Compare 需要使用 OpenRouter 作为提供商。这是此模块唯一支持的提供商。

`symphonyCompare` 是一个高级模块，它协调多个 AI 模型以比较细胞类型并自动建立共识。它并行使用多个 AI 模型进行全面的细胞类型比较，当模型对最佳匹配细胞类型存在分歧时，自动触发讨论轮次。可以把它想象成一个专家生物学家小组在辩论并达成共识。

## Quick Start

```r
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

cat("Consensus:", results$consensus, "\n")
cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
```

## Input

- **标记基因**: 来自 CASSIA 之前结果或您自己分析的字符向量或逗号分隔的标记基因字符串
- **候选细胞类型**: 要比较的 2-4 种细胞类型的向量
- **组织背景**: 正在分析的组织类型
- **物种**: 样本的物种（默认：人类）

## Parameters

### Required Parameters

| 参数 | 类型 | 描述 |
|-----------|------|-------------|
| `tissue` | 字符串 | 被分析的组织类型（例如，"blood", "brain", "liver"） |
| `celltypes` | 字符向量 | 要比较的 2-4 种细胞类型的向量 |
| `marker_set` | 字符向量/字符串 | 标记基因的列表或字符串 |

### Optional Parameters

| 参数 | 类型 | 默认值 | 描述 |
|-----------|------|---------|-------------|
| `species` | 字符串 | `"human"` | 样本的物种 |
| `model_preset` | 字符串 | `"premium"` | 要使用的预配置模型集合。`"premium"`: 高性能组合 (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4)。`"budget"`: 具成本效益的模型 (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash) |
| `enable_discussion` | 逻辑值 | `TRUE` | 如果为 TRUE，若未达成初始共识，模型将"讨论"并重新考虑其投票 |
| `generate_report` | 逻辑值 | `TRUE` | 是否生成分析的 HTML 报告 |
| `verbose` | 逻辑值 | `TRUE` | 是否向控制台打印进度消息 |

### Advanced Parameters

| 参数 | 类型 | 默认值 | 描述 |
|-----------|------|---------|-------------|
| `max_discussion_rounds` | 整数 | `2` | 允许的最大讨论轮数 |
| `consensus_threshold` | 数值 | `0.8` | 声明共识所需的置信度阈值（0-1） |

## Output

函数返回一个包含共识结果的列表，并生成输出文件。

**返回值:**
- `results`: 所有模型响应和分数的列表
- `consensus`: 共识细胞类型（如果达成）
- `confidence`: 共识的置信度水平（0-1）
- `csv_file`: 生成的包含详细结果的 CSV 文件路径
- `html_file`: 生成的交互式 HTML 报告路径（如果启用）
- `summary`: 比较的汇总统计信息
