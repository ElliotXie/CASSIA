---
title: 比较细胞类型（可选）
---


此功能允许您使用多个LLM确定哪种细胞类型最可能是簇的真实细胞类型。在默认设置下，会使用3个最先进的LLM基于标记基因对候选细胞类型进行评分。


## 函数参数

```r
compareCelltypes(
    tissue,        # 被分析的组织类型（例如，"large intestine"大肠）
    celltypes,     # 要比较的细胞类型向量（例如，c("Plasma Cells", "IgA-secreting Plasma Cells")）
    marker,        # 由逗号分隔的标记基因字符串
    species,       # 来源物种（"human"人类或"mouse"小鼠）
    output_file,   # 输出文件的名称（例如，"plasma_cell_subtype"）
    model_list     # 可选：要使用的LLM模型列表（有默认值）
)
```

## 参数详情

- `tissue`：指定数据的组织来源（例如，"large intestine"大肠，"small intestine"小肠，"brain"大脑）

- `celltypes`：您想要比较的细胞类型列表（最大推荐：4-5个）。示例：`c("Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells")`

- `marker`：逗号分隔的标记基因列表（例如，"IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3"）

- `species`：指定数据的物种来源

- `output_file`：输出文件的名称（不含扩展名）

- `model_list`：可选参数。默认模型（如果未提供）是最先进的LLM：
  ```r
  model_list = c(
      "anthropic/claude-3.5-sonnet",  # Anthropic的最新模型
      "openai/o1-mini",              # OpenAI的模型
      "google/gemini-pro-1.5"        # Google的模型
  )
  ```



## 输出格式

1. **控制台输出**：
   - 每个LLM对每种细胞类型的相似性评分
   - 共识结果（如果达成）
   - 警告消息（如果有）

2. **输出文件**（保存为"[output_file].txt"）：
   - 每个LLM的详细比较结果
   - 标记基因分析
   - 最终共识（如果达成）

## 解释指南

### 高置信度结果
- 当所有LLM对同一细胞类型给出80%以上的评分时，获得高置信度结果
- 这表明细胞类型识别清晰、明确

### 未达成共识
如果未达成明确共识，请考虑以下可能情况：

1. **低质量簇**
   - 症状：LLM之间的评分不一致或较低
   - 解决方案：增加分析中的标记基因数量

2. **混合簇**
   - 症状：不同的LLM强烈偏好不同的细胞类型
   - 解决方案：执行子聚类以分离潜在的不同群体

3. **最后手段**
   - 如果尝试上述解决方案后问题仍然存在
   - 咨询领域专家进行手动审查
