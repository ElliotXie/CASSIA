---
title: 快速模式
---

CASSIA的快速模式提供了一个简化的、单行的解决方案，用于运行完整的分析流程。此模式在单个函数调用中结合了注释、评分和注释增强，使用优化的默认参数。

### 基本用法
```R
runCASSIA_pipeline(
    output_file_name = "my_analysis",
    tissue = "brain",
    species = "human",
    marker = marker_data,
    max_workers = 4
)
```


### 将CASSIA结果添加回Seurat对象
```R
seurat_corrected <- add_cassia_to_seurat(
  seurat_obj = seurat_corrected, # 您想要添加CASSIA结果的Seurat对象
  cassia_results_path = "/FastAnalysisResults_scored.csv", #保存评分结果的位置，指定路径
  cluster_col = "celltype", # Seurat对象中包含细胞类型的列
  cassia_cluster_col="True Cell Type" # 评分结果中包含真实细胞类型的列
)

# 这将向Seurat对象添加六个新列：一般细胞类型、所有三个亚细胞类型、最可能的细胞类型、第二可能的细胞类型、第三可能的细胞类型和混合细胞类型，以及每种细胞类型的质量评分。
```




### 完整参数选项
```R
runCASSIA_pipeline(
    # 必需参数
    output_file_name,
    tissue,
    species,
    marker,
    
    # 带默认值的可选参数
    max_workers = 4,
    
    # 模型配置
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    annotation_reasoning = "medium",  # 可选: "high", "medium", "low"
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    annotationboost_model="anthropic/claude-sonnet-4.5",
    annotationboost_provider="openrouter",
    
    # 合并参数
    do_merge_annotations = TRUE,
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter",
    
    # 分析参数
    score_threshold = 75,
    additional_info = NULL,
    validator_involvement = "v1"
)
```

### 参数详情

- **`output_file_name`**: 输出文件夹和文件的基本名称。
- **`tissue`**: 样本的组织类型（例如，"brain"）。
- **`species`**: 样本的物种（例如，"human"）。
- **`marker`**: 标记基因数据（数据框或 CSV 路径）。
- **`max_workers`**: 并行处理的工作进程数。
- **`annotation_model`**: 用于初始细胞类型注释步骤的模型。
- **`annotation_reasoning`**: （可选）控制注释的推理深度（"high"、"medium"、"low"）。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。
- **`score_model`**: 用于质量评分的模型。**推荐**：使用像 `claude-4.5-sonnet` 这样高性能的模型以获得准确的评分。
- **`annotationboost_model`**: 用于优化低置信度注释的模型。
- **`do_merge_annotations`**: 逻辑值。如果为 `TRUE`，将详细的细胞类型合并为更广泛的类别。
- **`merge_model`**: 用于合并步骤的模型。
- **`score_threshold`**: 质量分数低于此阈值（0-100）的注释将触发注释增强过程。默认为 75。
- **`additional_info`**: 可选的实验上下文信息（例如，“药物 X 处理”）。
- **`validator_involvement`**: 控制验证严格程度（"v1" = 中等，"v0" = 高）。

### 主要功能说明

#### 合并注释
该流程包含一个自动合并步骤 (`do_merge_annotations = TRUE`)，将注释的簇分组为更广泛的类别（例如，将“CD4+ T 细胞”和“CD8+ T 细胞”分组为“T 细胞”）。这提供了细胞类型的分层视图，使您更容易了解数据集中的主要群体。

#### 验证器参与
`validator_involvement` 参数控制验证过程的强度：
- `"v0"`: 高参与度。应用更强的验证检查，可能较慢但更严格。
- `"v1"`: 中等参与度（默认）。适合大多数标准分析的平衡验证。

### 输出文件
流程生成一个文件夹，其中包含以下文件：
1. 注释结果 csv 文件
2. 评分结果 csv 文件
3. 合并注释结果（如果启用）
4. 基础 CASSIA 报告
5. 注释增强报告

### 性能提示
- 为获得最佳性能，根据系统的CPU核心调整`max_workers`
- 使用`additional_info`提供相关的实验上下文
- 监控`score_threshold`以平衡严格性和处理量



接下来我们详细介绍每个功能...
