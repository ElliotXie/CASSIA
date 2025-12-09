---
title: 批量处理
---

CASSIA支持批量处理，可同时分析多个簇。本篇我们将介绍如何准备数据并高效运行批量分析。


### 模型推荐
如果您使用OpenRouter作为提供商，可以指定模型如`"openai/gpt-4o-2024-11-20"`或`"anthropic/claude-3.5-sonnet"`。以下是一些模型推荐：

- **Claude 3.5 Sonnet**（最佳性能，略贵）
    - 模型ID：`"anthropic/claude-3.5-sonnet"`
- **GPT-4o**（平衡选项）
    - 模型ID：`"openai/gpt-4o-2024-11-20"`
- **Llama 3.2**（开源，经济实惠）
    - 模型ID：`"meta-llama/llama-3.2-90b-vision-instruct"`
- **Deepseek v3**（开源，几乎免费，性能与gpt4o相当，最推荐选项）
    - 模型ID：`"deepseek/deepseek-chat-v3-0324"`
    - 模型ID：`"deepseek/deepseek-chat-v3-0324:free"`

### 准备标记数据
您有四种提供标记数据的选项：

1. 创建包含簇和标记基因的数据框或CSV文件
2. 直接使用Seurat的`FindAllMarkers`函数输出
3. 使用Scanpy的`rank_genes_groups`输出（Python/导出的CSV）
4. 使用CASSIA的示例标记数据

```R
# 选项1：加载您自己的标记数据
markers <- read.csv("path/to/your/markers.csv")

# 选项2：直接使用Seurat的FindAllMarkers输出
#（假设您已经有一个Seurat对象）
markers <- FindAllMarkers(seurat_obj)

# 选项3：加载Scanpy rank_genes_groups输出（导出为CSV）
markers <- read.csv("scanpy_markers.csv")

# 选项4：加载示例标记数据
markers <- loadExampleMarkers()

# 预览数据
head(markers)
```

#### 标记数据格式
CASSIA接受三种格式：

**1. Seurat FindAllMarkers 输出（推荐）**

Seurat的`FindAllMarkers`函数的标准输出，包含差异表达统计信息：

```
p_val  avg_log2FC  pct.1  pct.2  p_val_adj  cluster  gene
0      3.02        0.973  0.152  0          0        CD79A
0      2.74        0.938  0.125  0          0        MS4A1
0      2.54        0.935  0.138  0          0        CD79B
0      1.89        0.812  0.089  0          1        IL7R
0      1.76        0.756  0.112  0          1        CCR7
```

**2. Scanpy rank_genes_groups 输出**

Scanpy的`sc.tl.rank_genes_groups()`函数的输出，通常使用`sc.get.rank_genes_groups_df()`导出：

```
group  names   scores  pvals  pvals_adj  logfoldchanges
0      CD79A   28.53   0      0          3.02
0      MS4A1   25.41   0      0          2.74
0      CD79B   24.89   0      0          2.54
1      IL7R    22.15   0      0          1.89
1      CCR7    20.87   0      0          1.76
```

**3. 简化格式**

包含簇ID和逗号分隔标记基因的两列数据框：

```
cluster  marker_genes
0        CD79A,MS4A1,CD79B,HLA-DRA,TCL1A
1        IL7R,CCR7,LEF1,TCF7,FHIT,MAL
2        CD8A,CD8B,GZMK,CCL5,NKG7
```

### 运行批量分析

#### 设置参数

```R
# 检测可用CPU核心
available_cores <- parallel::detectCores()

# 计算推荐的工作进程数（可用核心的75%）
recommended_workers <- floor(available_cores * 0.75)

runCASSIA_batch(
    # 必需参数
    marker = markers,                    # 标记数据（数据框或文件路径）
    output_name = "my_annotation",       # 输出文件的基本名称
    model = "anthropic/claude-4.5-sonnet", # 要使用的模型
    tissue = "brain",                    # 组织类型
    species = "human",                   # 物种

    # 可选参数
    max_workers = recommended_workers,    # 并行工作进程数
    n_genes = 50,                        # 要使用的顶部标记基因数量
    additional_info = "",                # 附加上下文
    provider = "openrouter",              # API提供商
    reasoning = "medium",                # 可选: "high", "medium", "low" 用于兼容模型

    # 高级参数
    ranking_method = "avg_log2FC",       # 排序方法："avg_log2FC"（默认）、"p_val_adj"、"pct_diff"或"Score"
    ascending = NULL,                    # 排序方向（NULL使用方法默认值）
    temperature = 0,                     # 模型温度
    celltype_column = NULL,              # 簇列名称（默认：第一列）
    gene_column_name = NULL,             # 基因列名称（默认：第二列）
    max_retries = 1,                     # 失败调用的最大重试次数
    validator_involvement = "v1"         # 验证器级别："v1"（中等）或"v0"（高）
)
```

### 参数详情

1. **标记基因选择**：
   - 默认：每个簇的前50个基因
   - `ranking_method`: 控制标记基因的排序和选择方式
     - `"avg_log2FC"` (默认): 按平均 log2 倍数变化排序
     - `"p_val_adj"`: 按调整后的 p 值排序
     - `"pct_diff"`: 按表达百分比差异排序
     - `"Score"`: 按自定义分数排序
   - 筛选标准：
     - 调整后的p值 < 0.05
     - 平均log2倍数变化 > 0.25
     - 最小百分比 > 0.1
   - 如果通过筛选的基因少于50个，则使用所有通过的基因

2. **并行处理**：
   - `max_workers`：控制并行处理线程
   - 推荐：可用CPU核心的80%
   - 示例：对于16核机器，设置为13

3. **附加上下文**（可选）：
   - 使用`additional_info`提供实验上下文
   - 示例：
     - 处理条件："Samples were antibody-treated"（样本经抗体处理）
     - 分析重点："Please carefully distinguish between cancer and non-cancer cells"（请仔细区分癌细胞和非癌细胞）
   - 提示：比较有无附加上下文的结果

4. **模型选择**：
   - 默认为 `anthropic/claude-4.5-sonnet` 以获得最佳性能。
   - 您可以使用 `google/gemini-2.5-flash` 快速查看数据。

5. **推理深度**（可选）：
   - 使用 `reasoning` 控制模型推理深度（"high"、"medium"、"low"）
   - 推荐 GPT-5.1 使用 `"medium"` - 高深度可能需要很长时间
   - 详见 [推理深度参数](setting-up-cassia.md#推理深度参数)

### 输出文件

分析会生成三个文件：

1. **`my_annotation_summary.csv`**：包含注释结果的清洁 CSV 文件：
   - 聚类 ID、预测的一般/详细细胞类型、可能的混合细胞类型
   - 标记数量和标记列表
   - 迭代次数、模型、提供商、组织、物种元数据

2. **`my_annotation_conversations.json`**：包含每个聚类完整对话历史的结构化 JSON 文件：
   - 注释尝试（如果验证失败，包含所有迭代）
   - 验证响应和格式化代理响应
   - 用于调试和理解推理过程

3. **`my_annotation_report.html`**：交互式 HTML 报告：
   - 摘要统计和可视化
   - 每个聚类的可展开对话历史
   - 可搜索和可过滤的结果表格

### 将 CASSIA 结果添加回 Seurat 对象

您可以使用 `add_cassia_to_seurat` 轻松将注释结果添加回您的 Seurat 对象。此函数根据簇标识符将 CASSIA 结果映射到您的 Seurat 对象。

```R
seurat_corrected <- add_cassia_to_seurat(
  seurat_obj = seurat_corrected,                         # 您想要添加 CASSIA 结果的 Seurat 对象
  cassia_results_path = "my_annotation_summary.csv",     # 摘要结果文件的路径
  cluster_col = "seurat_clusters",                       # Seurat 对象中包含簇 ID 的列
  cassia_cluster_col = "Cluster"                         # CASSIA 结果中包含簇 ID 的列
)
```

这将根据摘要 CSV 的列向您的 Seurat 对象添加多个新列：
- `Cluster ID`：簇标识符
- `Predicted General Cell Type`：一般细胞类型类别
- `Predicted Detailed Cell Type`：具体细胞类型预测
- `Possible Mixed Cell Types`：关于潜在混合群体的信息
- `Marker List`：用于注释的标记基因
- 其他元数据列：`Iterations`、`Model`、`Provider`、`Tissue`、`Species`

### 获取最佳结果的提示

1. **资源管理**：
   - 设置`max_workers`时监控系统资源
   - 从推荐的核心75%开始，根据需要调整

2. **标记基因选择**：
   - 默认的50个基因适用于大多数情况
   - 对于更复杂的细胞类型，可以增加基因数量
   - 如果遇到API速率限制，可以减少

3. **上下文优化**：
   - 测试运行有无附加上下文的情况
   - 保持上下文简洁相关
