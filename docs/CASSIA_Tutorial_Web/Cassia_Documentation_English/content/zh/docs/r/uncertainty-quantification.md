---
title: 不确定性量化（可选）
---


CASSIA中的不确定性量化通过多次分析迭代和相似性评分来评估注释可靠性。这个过程对以下方面至关重要：
- 识别稳健的细胞类型分配
- 检测混合或模糊的簇
- 量化注释置信度
- 理解预测变异性

### 多次迭代分析

#### 基本用法
```R
# 运行多次分析
runCASSIA_batch_n_times(
    n = 5,
    marker = marker_data,
    output_name = "my_annotation_repeat",
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    tissue = "brain",
    species = "human",
    max_workers = 4,
    batch_max_workers = 2
)
```
> **⚠️ 成本警告**：使用LLM模型运行多次迭代可能会产生显著成本。每次迭代都会进行单独的API调用，因此总成本将大约是单次运行成本的n倍。考虑从较小的迭代次数开始进行测试。

#### 参数详情

- **`n`**：分析迭代次数（推荐：5）。
- **`marker`**：标记基因数据（数据框或路径）。
- **`output_name`**：输出文件的基本名称。
- **`model`**：使用的 LLM 模型。
- **`provider`**：API 提供商。
- **`tissue`**：组织类型。
- **`species`**：物种。
- **`max_workers`**：整体并行处理限制。
- **`batch_max_workers`**：每次迭代的工作进程（max_workers * batch_max_workers 应与您的核心数量匹配）。

### 相似性分数计算

#### 运行相似性分析
```R
# 计算相似性分数
runCASSIA_similarity_score_batch(
    marker = marker_data,
    file_pattern = "my_annotation_repeat_*_full.csv",
    output_name = "similarity_results",
    max_workers = 4,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

#### 参数详情

- **`marker`**：标记基因数据。
- **`file_pattern`**：匹配迭代结果的模式（使用 `*` 通配符）。
    - 示例：`"my_annotation_repeat_*_full.csv"` 匹配 `my_annotation_repeat_1_full.csv`, `my_annotation_repeat_2_full.csv` 等。
- **`output_name`**：结果的基本名称。
- **`max_workers`**：并行工作进程数。
- **`model`**：用于评分的 LLM 模型。
- **`provider`**：API 提供商。
- **`main_weight`**：主细胞类型匹配的重要性（0-1）。
- **`sub_weight`**：亚类型匹配的重要性（0-1）。（权重总和应为 1.0）。

### 输出解释与故障排除

#### 相似性分数
- **范围**：0（完全不同）到 1（完全相同）。
- **> 0.9**：高一致性（稳健的注释）。
- **0.75 - 0.9**：中等一致性。
- **< 0.75**：低一致性（模糊）。

#### 低分故障排除
1. **审查数据**：检查标记基因质量和簇异质性。
2. **尝试高级智能体**：使用 [注释增强智能体](annotation-boost.md) 或 [子聚类](subclustering-analysis.md)。
3. **参数**：增加迭代次数。
