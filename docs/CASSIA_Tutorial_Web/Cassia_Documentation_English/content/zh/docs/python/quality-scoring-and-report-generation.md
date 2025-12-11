---
title: 质量评分
---

质量评分有助于评估细胞类型注释的可靠性。CASSIA 通过 `runCASSIA_score_batch` 函数提供自动评分功能，该函数分析每个注释背后的推理和证据。

### 运行质量评分

```python
# 运行质量评分
CASSIA.runCASSIA_score_batch(
    input_file = output_name + "_full.csv",
    output_file = output_name + "_scored.csv",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter",
    generate_report = True,  # 默认值，自动生成 HTML 报告
    reasoning = "medium"  # 可选: "low", "medium", "high"
)
```

### 参数详情

- **`input_file`**: 完整注释结果 CSV 文件的路径（由 `runCASSIA_batch` 生成）。
- **`output_file`**: 输出评分结果 CSV 文件名。
- **`max_workers`**: 并行评分线程数。
- **`model`**: 用于质量评分的 LLM。建议使用 `claude-4.5-sonnet` 或 `gpt-4o` 等高能力模型以获得准确评分。
- **`provider`**: 模型的 API 提供商（例如，"openrouter"）。
- **`generate_report`**: 是否自动生成 HTML 报告（默认：`True`）。报告保存为 `{output_file}_report.html`。
- **`reasoning`**: （可选）推理深度级别（"low"、"medium"、"high"）。控制模型在响应前"思考"的程度。仅支持 OpenAI GPT-5 系列模型（如 `gpt-5.1`）。通过 OpenRouter 使用无需额外验证。通过 OpenAI API 直接使用需要身份验证（KYC）。

### 解读评分

- **90-100**: 高置信度，证据充分。
- **76-89**: 良好置信度，证据充足。
- **<75**: 低置信度。这些聚类是使用注释增强智能体或比较智能体进行进一步分析的候选者。

HTML 报告会自动生成在 `{output_file}_report.html`，包含所有 CASSIA 输出，包括结构化结果、对话历史记录和质量评分。

