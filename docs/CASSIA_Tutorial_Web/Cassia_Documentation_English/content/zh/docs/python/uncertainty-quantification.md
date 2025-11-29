---
title: 不确定性量化（可选）
---

CASSIA 中的不确定性量化有助于通过多次分析迭代和相似性评分来评估注释的可靠性。此过程对于以下方面至关重要：
- 识别稳健的细胞类型分配
- 检测混合或模棱两可的聚类
- 量化注释置信度
- 了解预测的可变性

### 多重迭代分析

```python
# 运行多次迭代
iteration_results = CASSIA.runCASSIA_batch_n_times(
    n=2,
    marker=unprocessed_markers,
    output_name=output_name + "_Uncertainty",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    max_workers=6,
    batch_max_workers=3  # API 速率限制的保守设置
)
```

> **⚠️ 成本警告**：使用 LLM 模型运行多次迭代可能会产生大量费用。每次迭代都会进行单独的 API 调用，因此总成本约为单次运行成本的 n 倍。建议从较小的迭代次数开始进行测试。

#### 参数详情

- **`n`**: 分析迭代次数（建议：5）。
- **`marker`**: 标记基因数据（数据帧或路径）。
- **`output_name`**: 输出文件的基本名称。
- **`model`**: 要使用的 LLM 模型。
- **`provider`**: API 提供商。
- **`tissue`**: 组织类型。
- **`species`**: 物种。
- **`max_workers`**: 总体并行处理限制。
- **`batch_max_workers`**: 每次迭代的工作者数。

### 相似性评分计算

```python
# 计算相似性评分
similarity_scores = CASSIA.runCASSIA_similarity_score_batch(
    marker=unprocessed_markers,
    file_pattern=output_name + "_Uncertainty_*_full.csv",
    output_name="intestine_uncertainty",
    max_workers=6,
    model="openai/gpt-5.1",
    provider="openrouter",
    main_weight=0.5,
    sub_weight=0.5
)
```

#### 参数详情

- **`marker`**: 标记基因数据。
- **`file_pattern`**: 匹配迭代结果的模式（使用 `*` 通配符）。
    - 示例：`"my_annotation_repeat_*_full.csv"` 匹配 `my_annotation_repeat_1_full.csv`, `my_annotation_repeat_2_full.csv` 等。
- **`output_name`**: 结果的基本名称。
- **`max_workers`**: 并行工作者数。
- **`model`**: 用于评分的 LLM 模型。
- **`provider`**: API 提供商。
- **`main_weight`**: 主要细胞类型匹配的重要性 (0-1)。
- **`sub_weight`**: 亚型匹配的重要性 (0-1)。（权重总和应为 1.0）。

### 输出解读与故障排除

#### 相似性评分
- **范围**: 0（完全不同）到 1（完全相同）。
- **> 0.9**: 高一致性（稳健的注释）。
- **0.75 - 0.9**: 中等一致性。
- **< 0.75**: 低一致性（模棱两可）。

#### 低分故障排除
1. **检查数据**：检查标记基因质量和聚类异质性。
2. **尝试高级智能体**：使用 [注释增强智能体](annotation-boost.md) 或 [亚群聚类](subclustering-analysis.md)。
3. **参数**：增加迭代次数。

