---
title: 快速模式
---

CASSIA 的快速模式提供了一个精简的、单行解决方案来运行完整的分析流程。此模式在单个函数调用中结合了注释、评分和用于纠正低质量注释的注释增强，使用优化的默认参数。

### 基本用法

```python
# 以快速模式运行 CASSIA 流程
CASSIA.runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker_path = unprocessed_markers,
    max_workers = 6,  # 与数据集中的聚类数匹配
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model="anthropic/claude-sonnet-4.5",
    annotationboost_provider="openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

### 参数详情

- **`output_file_name`**: 输出文件夹和文件的基本名称。
- **`tissue`**: 样本的组织类型（例如，“脑”）。
- **`species`**: 样本的物种（例如，“人类”）。
- **`marker_path`**: 标记基因数据（数据帧或 CSV 路径）。
- **`max_workers`**: 使用的并行进程数。
- **`annotation_model`**: 用于初始细胞类型注释步骤的模型。
- **`annotation_provider`**: 注释模型的提供商。
- **`score_model`**: 用于质量评分的模型。**建议**：使用像 `claude-4.5-sonnet` 这样的高能力模型以获得准确的评分。
- **`score_provider`**: 评分模型的提供商。
- **`score_threshold`**: 质量评分低于此阈值 (0-100) 的注释将触发注释增强过程。默认值为 75。
- **`annotationboost_model`**: 用于改进低置信度注释的模型。
- **`annotationboost_provider`**: 注释增强模型的提供商。
- **`do_merge_annotations`**: 逻辑值。如果为 `True`，则将详细的细胞类型合并为更广泛的类别。
- **`merge_model`**: 用于合并步骤的模型。
- **`merge_provider`**: 合并模型的提供商。
- **`additional_info`**: 可选的实验背景（例如，“用药物 X 处理”）。
- **`validator_involvement`**: 控制验证的严格程度（"v1" = 中等，"v0" = 高）。

### 主要功能说明

#### 合并注释
该流程包括一个自动合并步骤 (`do_merge_annotations = True`)，将注释的聚类分组为更广泛的类别（例如，将“CD4+ T 细胞”和“CD8+ T 细胞”分组为“T 细胞”）。这提供了细胞类型的分层视图，使您更容易理解数据集中的主要群体。

#### 验证者参与度
`validator_involvement` 参数控制验证过程的强度：
- `"v0"`: 高参与度。应用更强的验证检查，可能较慢但更严格。
- `"v1"`: 中等参与度（默认）。适合大多数标准分析的平衡验证。

### 输出文件
该流程生成一个包含以下文件的文件夹：
1. 注释结果 csv 文件
2. 评分结果 csv 文件
3. 合并的注释结果（如果启用）
4. 基本 CASSIA 报告
5. 注释增强报告

### 性能提示
- 为了获得最佳性能，请根据系统的 CPU 核心数调整 `max_workers`
- 使用 `additional_info` 提供相关的实验背景
- 监控 `score_threshold` 以平衡严格性和吞吐量

