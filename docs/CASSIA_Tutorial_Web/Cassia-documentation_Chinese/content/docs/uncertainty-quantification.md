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
    # 核心参数
    n = 5, #迭代次数
    marker = marker_data,
    output_name = "my_annottaion_repeat",
    
    # 模型设置
    model = "gpt-4o,
    provider = "openai",
    
    # 上下文信息
    tissue = "brain",
    species = "human",
    additional_info = NULL,


    # 处理控制
    max_workers = 4,        # 总并行工作进程
    batch_max_workers = 2   # 每批次的工作进程
)
```

> **⚠️ 成本警告**：使用LLM模型运行多次迭代可能会产生显著成本。每次迭代都会进行单独的API调用，因此总成本将大约是单次运行成本的n倍。考虑从较小的迭代次数开始进行测试。

#### 参数详情

1. **迭代控制**：
   - `n`：分析迭代次数
   - 推荐：标准分析使用5次迭代
   - 对于关键应用考虑更多迭代

2. **资源管理**：
   - `max_workers`：整体并行处理限制
   - `batch_max_workers`：每次迭代的工作进程
   - max_workers * batch_max_workers 应与您的核心数量匹配。


### 相似性分数计算

#### 运行相似性分析
```R

# 计算相似性分数
runCASSIA_similarity_score_batch(
    # 输入参数
    marker = marker_data,
    file_pattern = "my_annottaion_repeat_*_full.csv",
    output_name = "similarity_results",
    
    
    # 处理参数
    max_workers = 4,
    model = "anthropic/claude-3.5-sonnet",
    provider = "openrouter",
    
    # 评分权重
    main_weight = 0.5, # 主细胞类型匹配的权重
    sub_weight = 0.5  # 亚类型匹配的权重
)
```

#### 评分参数

1. **权重配置**：
   - `main_weight`：主细胞类型匹配的重要性（0-1）
   - `sub_weight`：亚类型匹配的重要性（0-1）
   - 权重总和应为1.0

2. **文件管理**：
   - `file_pattern`：匹配迭代结果的模式
   - 使用*匹配迭代数字
   - 示例：如果您有"my_annottaion_repeat_1_full.csv"、"my_annottaion_repeat_2_full.csv"和"my_annottaion_repeat_3_full.csv"，使用"my_annottaion_repeat__full.csv"来匹配模式。

#### 输出解释

1. **相似性分数**：
   - 范围：0（完全不同）到1（完全相同）
   - 解释指南：
     - 0.9：高一致性
     - 0.75-0.9：中等一致性
     - <0.75：低一致性

#### 故障排除

1. **性能问题**：
   - 减少工作进程数量
   - 以较小批次处理

2. **低相似性分数**：
   - 审查标记基因质量
   - 使用注释增强功能
   - 审查簇异质性
   - 考虑生物变异性
   - 增加迭代次数
   - 尝试子聚类
