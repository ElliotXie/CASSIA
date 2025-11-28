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
    output_file_name,     # 输出文件的基本名称
    tissue,               # 组织类型（例如，"brain"脑，"blood"血液）
    species,              # 物种（例如，"human"人类，"mouse"小鼠）
    marker,               # 来自findallmarker的标记文件，路径或数据对象
    
    # 带默认值的可选参数
    max_workers = 4,      # 并行工作进程数
    
    # 模型配置
    annotation_model = "gpt-4o",             # 用于注释的模型
    annotation_provider = "openai",         # 注释的提供商
    score_model = "anthropic/claude-3.5-sonnet",  # 用于评分的模型
    score_provider = "openrouter",         # 评分的提供商
    annotationboost_model="anthropic/claude-3.5-sonnet", #用于注释增强的模型
    annotationboost_provider="openrouter", #注释增强的提供商
    
    # 分析参数
    score_threshold = 75,     # 最低可接受分数
    additional_info = NULL    # 附加上下文信息
)
```


### 输出文件
流程生成：
1. 初始注释结果
2. 质量分数和推理
3. 摘要报告
4. 注释增强报告

### 性能提示
- 为获得最佳性能，根据系统的CPU核心调整`max_workers`
- 使用`additional_info`提供相关的实验上下文
- 监控`score_threshold`以平衡严格性和处理量



接下来我们详细介绍每个功能...
