---
title: 注释增强（可选）
---


注释增强是一种高级验证工具，通过多次迭代分析增强注释的置信度。它特别适用于：
- 验证低置信度注释
- 获取特定细胞簇的详细见解
- 解决模糊的细胞类型分配
- 生成全面的验证报告

### 所需组件

1. **输入数据**：
   - CASSIA批量分析的完整结果CSV
   - 原始标记基因文件（Seurat输出或自定义标记文件）
   - 簇上下文信息
   - 特定簇标识符

2. **模型配置**：
   - 使用至少比 `gemini-2.5-flash` 更先进的模型
   - 推荐：通过 OpenRouter 或 Anthropic 使用 `anthropic/claude-sonnet-4.5`

### 运行注释增强

```R
# 设置参数
validation_config <- list(
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)

# 定义簇信息
cluster_info <- "Human PBMC"

# 指定要验证的簇
target_cluster = "CD4+ T cell"

# 运行验证
runCASSIA_annotationboost(
    # 必需参数
    full_result_path = "batch_results_summary.csv",
    marker = marker_data,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "Cluster1_report",

    # 可选参数
    num_iterations = 5,             # 验证轮数
    model = validation_config$model,
    provider = validation_config$provider,

    # 对话历史（来自批量注释）
    conversations_json_path = "batch_results_conversations.json",
    conversation_history_mode = "full",  # "full"、"final"（摘要）或 "none"

    # 高级选项
    search_strategy = "breadth",    # "breadth"（广度）或 "depth"（深度）
    report_style = "per_iteration", # "per_iteration"（每次迭代）或 "total_summary"（总体摘要）
    validator_involvement = "v1",   # "v0"（高）或 "v1"（中等）
    reasoning = "low"               # 可选: "low", "medium", "high"
)
```

### 参数详情

- **`full_result_path`**: CASSIA 结果 CSV 文件的路径（`_summary.csv`）。
- **`marker`**: 标记基因数据（数据框或路径）。**重要**：请使用与初始分析相同的标记数据（不要过滤）。
- **`cluster_name`**: 要验证的目标簇的确切名称。
- **`major_cluster_info`**: 数据集的上下文信息（例如，"Human PBMC", "Mouse Brain"）。
- **`output_name`**: 输出验证报告的基本名称。
- **`num_iterations`**: 验证轮数（默认：5）。
- **`model`**: 用于验证的 LLM 模型。
- **`provider`**: 模型的 API 提供商。
- **`conversations_json_path`**: 批量注释的对话 JSON 文件路径（`_conversations.json`）。这为深入分析提供了先前的注释上下文。
- **`conversation_history_mode`**: 如何使用先前的对话历史。
    - `"full"`（默认）：使用完整的先前对话历史作为上下文。
    - `"final"`：在使用之前先用 LLM 对历史进行摘要。
    - `"none"`：不使用任何先前的对话历史。
- **`search_strategy`**: 探索假设的策略。
    - `"breadth"`（默认）：并行测试多个假设。
    - `"depth"`：深入关注验证单个假设。
- **`report_style`**: 最终报告的格式。
    - `"per_iteration"`：每轮的详细分解。
    - `"total_summary"`：简明概述。
- **`validator_involvement`**: 验证严格程度级别。
    - `"v1"`（默认）：中等参与度。
    - `"v0"`：高参与度（更严格的检查）。
- **`reasoning`** (字符串, 可选): 推理深度级别（"low"、"medium"、"high"）。控制模型在响应前"思考"的程度。仅支持 OpenAI GPT-5 系列模型（如 `gpt-5.1`）。通过 OpenRouter 使用无需额外验证。通过 OpenAI API 直接使用需要身份验证（KYC）。

### 输出文件

分析生成以下输出文件：
- `{output_name}_summary.html`：包含详细分析结果和可视化的 HTML 报告。
- `{output_name}_raw_conversation.txt`：包含完整分析对话的原始对话文本。

### 故障排除

1. **低置信度结果**：
   - 增加 `num_iterations` 以收集更多证据。
   - 检查标记基因质量。
   - 检查是否存在双细胞、批次效应或背景污染等潜在问题。

2. **不一致结果**：
   - 验证输入数据的质量和一致性。
   - 考虑生物变异性或混合群体（可能需要子聚类）。
