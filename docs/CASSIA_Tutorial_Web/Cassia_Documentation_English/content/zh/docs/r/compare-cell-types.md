---
title: Symphony Compare (可选)
---


`symphonyCompare` 是一个高级模块，它协调多个 AI 模型以比较细胞类型并自动建立共识。它并行使用多个 AI 模型进行全面的细胞类型比较，当模型对最佳匹配细胞类型存在分歧时，自动触发讨论轮次。可以把它想象成一个专家生物学家小组在辩论并达成共识。

### 使用示例

```r
# 基本用法 - 让 Symphony Compare 处理一切
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

# 访问结果
cat("Consensus:", results$consensus, "\n")
cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
```

### 函数参数

```r
symphonyCompare(
    tissue,             # 组织类型
    celltypes,          # 要比较的2-4种细胞类型的向量
    marker_set,         # 标记基因的向量或字符串
    species = "human",  # 物种
    model_preset = "premium", # 预设模型配置
    enable_discussion = TRUE,  # 启用自动讨论轮次
    max_discussion_rounds = 2, # 最大讨论轮次
    consensus_threshold = 0.8, # 共识阈值 (0-1)
    generate_report = TRUE,    # 生成 HTML 报告
    verbose = TRUE             # 打印进度消息
)
```

### 参数详情

- **`tissue`** (字符串): 被分析的组织类型（例如，"blood", "brain", "liver"）。
- **`celltypes`** (字符向量): 要比较的 2-4 种细胞类型的向量。
- **`marker_set`** (字符向量或字符串): 标记基因的列表或字符串。
- **`species`** (字符串): 样本的物种（默认："human"）。
- **`model_preset`** (字符串): 要使用的预配置模型集合。
  - `"premium"`: 高性能组合 (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4)
  - `"budget"`: 具成本效益的模型 (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash)
- **`enable_discussion`** (逻辑值): 如果为 `TRUE`，若未达成初始共识，模型将“讨论”并重新考虑其投票。
- **`max_discussion_rounds`** (整数): 允许的最大讨论轮数（默认：2）。
- **`consensus_threshold`** (数值): 声明共识所需的置信度阈值（0-1，默认：0.8）。
- **`generate_report`** (逻辑值): 是否生成分析的 HTML 报告（默认：`TRUE`）。
- **`verbose`** (逻辑值): 是否向控制台打印进度消息（默认：`TRUE`）。

### 输出格式

函数返回一个包含以下内容的列表：
- **`results`**: 所有模型响应和分数的列表。
- **`consensus`**: 共识细胞类型（如果达成）。
- **`confidence`**: 共识的置信度水平 (0-1)。
- **`csv_file`**: 生成的包含详细结果的 CSV 文件路径。
- **`html_file`**: 生成的交互式 HTML 报告路径（如果启用）。
- **`summary`**: 比较的汇总统计信息。
