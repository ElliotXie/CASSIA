---
title: 不确定性量化（可选）
---

CASSIA 中的不确定性量化有助于通过多次分析迭代和相似性评分来评估注释的可靠性。此过程对于以下方面至关重要：
- 识别稳健的细胞类型分配
- 检测混合或模棱两可的聚类
- 量化注释置信度
- 了解预测的可变性

### 单聚类不确定性分析

对于分析单个聚类的不确定性，使用 `runCASSIA_n_times_similarity_score()`：

```python
from CASSIA import runCASSIA_n_times_similarity_score

# 对单个聚类运行多次迭代并计算相似性评分
result = runCASSIA_n_times_similarity_score(
    tissue="large intestine",
    species="human",
    marker_list=["CD38", "CD138", "JCHAIN", "MZB1", "SDC1"],
    model="google/gemini-2.5-flash",
    provider="openrouter",
    n=5,  # 迭代次数
    temperature=0.3,
    max_workers=3,
    main_weight=0.5,
    sub_weight=0.5,
    validator_involvement="v1"
)

# 访问结果
print(f"主要细胞类型: {result['general_celltype_llm']}")
print(f"亚细胞类型: {result['sub_celltype_llm']}")
print(f"相似性评分: {result['similarity_score']}")
print(f"共识类型: {result['consensus_types']}")

# 检查混合细胞类型
if result.get('Possible_mixed_celltypes_llm'):
    print(f"可能的混合类型: {result['Possible_mixed_celltypes_llm']}")
```

#### 参数详情（单聚类）

- **`tissue`**: 用于上下文的组织类型。
- **`species`**: 用于上下文的物种。
- **`marker_list`**: 聚类的标记基因列表。
- **`model`**: 要使用的 LLM 模型。
- **`provider`**: API 提供商（"openrouter"、"openai"、"anthropic"）。
- **`n`**: 分析迭代次数（默认：5）。
- **`temperature`**: LLM 温度（较低 = 更一致）。
- **`max_workers`**: 并行处理工作者数。
- **`main_weight`**: 相似性中主要细胞类型的权重 (0-1)。
- **`sub_weight`**: 相似性中亚型的权重 (0-1)。
- **`validator_involvement`**: 验证器模式（"v0" 严格，"v1" 中等）。
- **`additional_info`**: 可选的额外上下文字符串。

#### 返回值（单聚类）

函数返回包含以下内容的字典：
- **`general_celltype_llm`**: 共识主要细胞类型。
- **`sub_celltype_llm`**: 共识亚细胞类型。
- **`similarity_score`**: 跨迭代的总体相似性 (0-1)。
- **`consensus_types`**: 出现频率最高的细胞类型。
- **`Possible_mixed_celltypes_llm`**: 检测到的混合细胞类型群体。
- **`original_results`**: 每次迭代的原始结果。

### 批量迭代分析

对于跨多个聚类的批处理，使用 `runCASSIA_batch_n_times`：

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

