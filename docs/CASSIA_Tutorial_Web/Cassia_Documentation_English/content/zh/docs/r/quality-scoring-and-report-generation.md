---
title: 质量评分
---

## 概述

质量评分通过分析每个预测背后的推理和证据来评估细胞类型注释的可靠性。评分智能体为每个注释分配置信度分数（0-100），帮助识别可能需要通过注释增强智能体或比较智能体进一步验证的聚类。

---

## 快速开始

```R
runCASSIA_score_batch(
    input_file = "my_annotation_summary.csv",  # JSON 自动检测
    output_file = "my_annotation_scored.csv",
    model = "openai/gpt-5.1",
    provider = "openrouter"
)
```

---

## 输入

质量评分需要由 `runCASSIA_batch` 生成的**两个文件**：

1. **摘要 CSV** (`*_summary.csv`) - 包含聚类标识符、预测的细胞类型、标记基因列表和元数据
2. **对话 JSON** (`*_conversations.json`) - 包含完整的对话历史和推理过程

**自动检测：** 当您提供摘要 CSV 作为输入时，CASSIA 会自动在同一目录中查找对应的对话 JSON 文件。无需手动指定两个文件。

---

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `input_file` | 摘要 CSV 的路径（来自 `runCASSIA_batch`）。对话 JSON 会自动检测。 |
| `model` | 用于评分的 LLM 模型 ID |
| `provider` | API 提供商（`"openrouter"`、`"openai"`、`"anthropic"`）或自定义 base URL |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `output_file` | `{input}_scored.csv` | 评分结果的输出文件名 |
| `conversations_json_path` | 自动检测 | 对话 JSON 的路径。默认从输入 CSV 文件名自动检测。 |
| `max_workers` | 4 | 并行评分线程数 |
| `reasoning` | `NULL` | 仅适用于通过 OpenRouter 使用的 GPT-5 系列模型的推理深度（`"low"`、`"medium"`、`"high"`）。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。 |

### 参数详情

**模型选择**
- 强烈推荐：`openai/gpt-5.1` 或 `anthropic/claude-sonnet-4.5` 以获得最佳准确性
- 评分需要强大的推理能力来评估注释质量
- 使用与注释运行相同的提供商设置

---

## 输出

### 生成的文件

| 文件 | 描述 |
|------|------|
| `{output_file}` | 包含质量分数和推理的评分结果 CSV |
| `{output_file}_report.html` | 包含所有 CASSIA 输出的交互式 HTML 报告 |

### 输出内容

评分 CSV 文件包括：
- 原始注释数据
- 质量分数（0-100）
- 每个分数的详细推理

### 解读分数

| 分数范围 | 置信度 | 建议操作 |
|----------|--------|----------|
| 90-100 | 高 | 证据充分，注释可靠 |
| 76-89 | 良好 | 证据充足，一般可信 |
| <75 | 低 | 通过注释增强智能体和比较智能体运行 |
