---
title: 注释增强（可选）
---


注释增强是一种高级验证工具，通过多次迭代分析增强注释的置信度。它特别适用于：
- 验证低置信度注释
- 获取特定细胞簇的详细见解
- 解决模糊的细胞类型分配
- 生成全面的验证报告

### 所需组件

1. **输入数据**：
   - CASSIA批量分析的完整结果CSV
   - 原始标记基因文件（Seurat输出或自定义标记文件）
   - 簇上下文信息
   - 特定簇标识符

2. **模型配置**：
   - 推荐：通过OpenRouter使用`anthropic/claude-3.5-sonnet`
   - 可以使用其他模型，但可能提供不同级别的详细信息

### 运行注释增强

```R
# 设置参数
validation_config <- list(
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter"
)

# 定义簇信息
cluster_info <- "Human PBMC"

# 指定要验证的簇
target_cluster = "CD4+ T cell"

# 运行验证
runCASSIA_annotationboost(

    # 必需参数
    full_result_path = "cell_type_analysis_results.csv",
    marker = marker_data,
    cluster_name = target_cluster,
    major_cluster_info = cluster_info,
    output_name = "Cluster1_report",
    num_iterations = 5, # 验证轮数

    # 模型配置
    model = validation_config$model,
    provider = validation_config$provider,
)
```

### 参数详情

   - `full_result_path`：原始CASSIA结果的路径
   - `marker`：标记基因数据（与初始分析中使用的相同）
   - `cluster_name`：目标簇名称
   - `major_cluster_info`：数据集上下文
   - `num_iterations`：验证轮数（默认：5）


### 故障排除

1. **低置信度结果**：
   - 增加`num_iterations`
   - 检查标记基因质量
   - 检查簇的质量，包括去双胞，去批次，去背景污染等

3. **不一致结果**：
   - 检查标记基因一致性
   - 验证输入数据质量
   - 考虑生物变异性
