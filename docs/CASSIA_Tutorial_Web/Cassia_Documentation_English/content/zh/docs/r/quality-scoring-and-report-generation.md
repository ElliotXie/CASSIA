---
title: 质量评分
---


质量评分有助于评估细胞类型注释的可靠性。CASSIA通过`runCASSIA_score_batch`函数提供自动评分功能，该功能分析每个注释背后的推理和证据。

## 运行质量评分

### 基本用法
```R
runCASSIA_score_batch(
    input_file = "my_annotation_full.csv",
    output_file = "my_annotation_scored.csv",
    max_workers = 4,
    model = "openai/gpt-5.1", # 推荐用于精确评分
    provider = "openrouter",
    reasoning = "low"  # 可选: "low", "medium", "high"
)
```

### 参数详情

- **`input_file`** (字符串): 完整注释结果 CSV 文件的路径（由 `runCASSIA_batch` 生成）。
- **`output_file`** (字符串): 输出评分结果 CSV 文件的名称。
- **`max_workers`** (整数): 并行评分线程数。
- **`model`** (字符串): 用于质量评分的 LLM。强烈推荐使用 `openai/gpt-5.1` 或 `anthropic/claude-4.5-sonnet` 以获得最佳准确性。
- **`provider`** (字符串): 模型的 API 提供商（例如 "openrouter"）。
- **`reasoning`** (字符串, 可选): 推理深度级别（"low"、"medium"、"high"）。控制模型在响应前"思考"的程度。仅支持 OpenAI GPT-5 系列模型（如 `gpt-5.1`）。通过 OpenRouter 使用无需额外验证。通过 OpenAI API 直接使用需要身份验证（KYC）。

### 输出格式
评分输出文件包含：
- 原始注释数据
- 质量分数（0-100）
- 置信度指标
- 分数的详细推理

### 解读分数

- **90-100**：高置信度，强有力的证据
- **76-89**：良好的置信度，充分的证据
- **<75**：低置信度，需要通过注释增强智能体和比较智能体运行

HTML 报告会自动生成在 `{output_file}_report.html`，包含所有 CASSIA 输出，包括结构化结果、对话历史和质量分数。
