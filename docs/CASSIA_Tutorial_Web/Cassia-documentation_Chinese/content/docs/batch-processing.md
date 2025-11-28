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
您有三种提供标记数据的选项：

1. 创建包含簇和标记基因的数据框或CSV文件
2. 直接使用Seurat的`findAllMarkers`函数输出
3. 使用CASSIA的示例标记数据

```R
# 选项1：加载您自己的标记数据
markers <- read.csv("path/to/your/markers.csv")

# 选项2：直接使用Seurat的findAllMarkers输出
#（假设您已经有一个Seurat对象）
markers <- FindAllMarkers(seurat_obj)

# 选项3：加载示例标记数据
markers <- loadExampleMarkers()

# 预览数据
head(markers)
```

#### 标记数据格式
CASSIA接受两种格式：

1. **FindAllMarkers输出**：Seurat的FindAllMarkers函数的标准输出
2. **简化格式**：一个两列数据框，其中：
   - 第一列：簇标识符
   - 第二列：逗号分隔的排序标记基因

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
    model = "gpt-4o",                     # 要使用的模型
    tissue = "brain",                    # 组织类型
    species = "human",                   # 物种
    
    # 可选参数
    max_workers = recommended_workers,    # 并行工作进程数
    n_genes = 50,                        # 要使用的顶部标记基因数量
    additional_info = "",                # 附加上下文
    provider = "openai"                  # API提供商
)
```

### 参数详情

1. **标记基因选择**：
   - 默认：每个簇的前50个基因
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

### 输出文件

分析会生成两个文件：
1. `my_annotation_full.csv`：完整对话历史
2. `my_annotation_summary.csv`：简明结果摘要

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
