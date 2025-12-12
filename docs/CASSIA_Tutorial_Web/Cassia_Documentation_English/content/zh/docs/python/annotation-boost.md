---
title: 注释增强（可选）
---

## 概述

注释增强是一种高级验证工具，通过多次分析迭代增强注释置信度。它特别适用于：
- 验证低置信度注释
- 深入了解特定细胞聚类
- 解决模棱两可的细胞类型分配
- 生成全面的验证报告

## 快速开始

```python
CASSIA.runCASSIA_annotationboost(
    full_result_path = "batch_results_summary.csv",
    marker = marker_data,
    cluster_name = "CD4+ T cell",
    major_cluster_info = "Human PBMC",
    output_name = "Cluster1_report",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter",
)
```

## 输入

- 来自 CASSIA 批量分析的完整结果 CSV（`_summary.csv`）
- 原始标记基因文件（***注意：标记文件不应被过滤！***）
- 聚类背景信息
- 特定聚类标识符
- （可选）批量注释的对话 JSON 文件（`_conversations.json`）

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `full_result_path` | CASSIA 结果 CSV 文件的路径（`_summary.csv`） |
| `marker` | 标记基因数据（数据帧或路径）。使用与初始分析相同的标记数据（不要过滤） |
| `cluster_name` | 要验证的目标聚类的确切名称 |
| `major_cluster_info` | 关于数据集的背景信息（例如，"Human PBMC"，"Mouse Brain"） |
| `output_name` | 输出验证报告的基本名称 |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `num_iterations` | 5 | 验证轮数 |
| `model` | - | 使用的 LLM 模型。推荐：`anthropic/claude-sonnet-4.5` 或更好的模型 |
| `provider` | - | 模型的 API 提供商 |
| `conversations_json_path` | - | 批量注释的对话 JSON 文件路径 |
| `conversation_history_mode` | `"full"` | 如何使用先前的对话历史：`"full"`、`"final"` 或 `"none"` |

### 高级参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `search_strategy` | `"breadth"` | 探索假设的策略：`"breadth"` 或 `"depth"` |
| `report_style` | `"per_iteration"` | 最终报告的格式：`"per_iteration"` 或 `"total_summary"` |
| `validator_involvement` | `"v1"` | 验证严格程度：`"v1"`（中等）或 `"v0"`（高） |
| `reasoning` | - | 推理深度级别：`"low"`、`"medium"`、`"high"`。仅支持 OpenAI GPT-5 系列模型 |

## 输出

分析生成以下输出文件：
- `{output_name}_summary.html`：包含详细分析结果和可视化的 HTML 报告。
- `{output_name}_raw_conversation.txt`：包含完整分析对话的原始对话文本。
