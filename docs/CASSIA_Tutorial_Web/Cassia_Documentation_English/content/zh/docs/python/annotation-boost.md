---
title: 注释增强（可选）
---

注释增强是一种高级验证工具，通过多次分析迭代增强注释置信度。它特别适用于：
- 验证低置信度注释
- 深入了解特定细胞聚类
- 解决模棱两可的细胞类型分配
- 生成全面的验证报告

### 所需组件

1. **输入数据**：
   - 来自 CASSIA 批量分析的完整结果 CSV
   - 原始标记基因文件（***注意：标记文件不应被过滤！***）
   - 聚类背景信息
   - 特定聚类标识符

2. **模型配置**：
   - 使用至少比 `gemini-2.5-flash` 更先进的模型
   - 推荐：通过 OpenRouter 或 Anthropic 使用 `anthropic/claude-sonnet-4.5`

### 运行注释增强

单核细胞聚类有时被注释为免疫细胞和神经元/胶质细胞的混合群体。这里我们使用注释增强智能体来更详细地测试这些假设。

```python
# 对高线粒体含量聚类运行增强验证
CASSIA.runCASSIA_annotationboost(
    full_result_path = output_name + "_full.csv",
    marker = unprocessed_markers,
    output_name = "monocyte_annotationboost",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

### 参数详情

- **`full_result_path`**: 原始 CASSIA 结果 CSV 文件的路径。
- **`marker`**: 标记基因数据（数据帧或路径）。**重要**：使用与初始分析相同的标记数据（不要过滤）。
- **`cluster_name`**: 要验证的目标聚类的确切名称。
- **`major_cluster_info`**: 关于数据集的背景信息（例如，“人类 PBMC”，“小鼠脑”）。
- **`output_name`**: 输出验证报告的基本名称。
- **`num_iterations`**: 验证轮数（默认：5）。
- **`model`**: 用于验证的 LLM 模型。
- **`provider`**: 模型的 API 提供商。
- **`search_strategy`**: 探索假设的策略（"breadth" 或 "depth"）。
- **`report_style`**: 最终报告的格式（"per_iteration" 或 "total_summary"）。
- **`validator_involvement`**: 验证严格程度（"v1" 或 "v0"）。

### 输出文件

分析生成以下输出文件：
- `{output_name}_summary.html`：包含详细分析结果和可视化的 HTML 报告。
- `{output_name}_raw_conversation.txt`：包含完整分析对话的原始对话文本。

### 故障排除

1. **低置信度结果**：
   - 增加 `num_iterations` 以收集更多证据。
   - 检查标记基因质量。
   - 检查是否存在双峰、批次效应或背景污染等潜在问题。

2. **结果不一致**：
   - 验证输入数据的质量和一致性。
   - 考虑生物变异性或混合群体（可能需要亚群聚类）。

