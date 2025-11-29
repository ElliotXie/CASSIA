---
title: 亚群聚类分析（可选）
---

亚群聚类分析是一种更详细地研究特定细胞群体的强大技术。本教程将带您了解分析亚群（如 T 细胞或成纤维细胞）的过程。

### 工作流程摘要
1. 初始 CASSIA 分析
2. 亚群提取和处理（使用 Seurat 或 Scanpy）
3. 标记识别
4. CASSIA 亚群聚类分析
5. 不确定性评估（可选）

### 运行亚群聚类分析

我们建议首先应用默认的 CASSIA。然后，在目标聚类上，应用标准流程（Seurat/Scanpy）进行亚群聚类并获取标记结果。

```python
CASSIA.runCASSIA_subclusters(
    marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

#### 参数详情

- **`marker`**: 亚群的标记基因（数据帧或文件路径）。
- **`major_cluster_info`**: 父聚类或背景的描述（例如，“CD8+ T 细胞”或“与其他细胞类型混合的 cd8 t 细胞”）。
- **`output_name`**: 输出 CSV 文件的基本名称。
- **`model`**: 要使用的 LLM 模型。
- **`provider`**: API 提供商。
- **`temperature`**: 采样温度 (0-1)。
- **`n_genes`**: 要使用的顶部标记基因数。

### 不确定性评估

为了获得更可信的结果，计算一致性分数 (CS)：

```python
# 运行多次迭代
CASSIA.runCASSIA_n_subcluster(
    n=5, 
    marker=subcluster_results,
    major_cluster_info="cd8 t cell", 
    base_output_name="subclustering_results_n",
    model="anthropic/claude-sonnet-4.5",
    temperature=0,
    provider="openrouter",
    max_workers=5,
    n_genes=50
)

# 计算相似性评分
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

### 输出文件
- `{output_name}.csv`: 基本 Cassia 分析结果。
- `{output_name}_uncertainty.csv`: 相似性评分（如果进行了不确定性评估）。

