---
title: "使用标记文件进行基础注释"
---

本教程介绍使用CASSIA进行细胞类型注释的基本步骤，适用于您已经准备好标记基因列表的情况。这非常适合您已经完成聚类分析并希望对细胞簇进行注释的场景。

## 1. 安装和设置

在使用CASSIA之前，您需要安装软件包并设置必要的环境。

### 1.1 软件包安装

```r
# 安装前置依赖包
install.packages("reticulate")
install.packages("devtools")

# 从GitHub安装CASSIA
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")

# 加载CASSIA包
library(CASSIA)
```

### 1.2 Python环境设置

CASSIA依赖Python进行部分后端操作。当您加载软件包时，环境会自动设置。如果遇到问题，请运行以下命令手动设置Python环境：

```r
# 手动设置Python环境
setup_cassia_env(conda_env = "cassia_env")
```

### 1.3 API密钥配置

CASSIA至少需要一个API密钥才能运行。我们建议设置OpenRouter、OpenAI和Anthropic的API密钥以获得最佳体验：

```r
# 设置API密钥（请替换为您的实际密钥）
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
```

设置 `persist = TRUE` 会将密钥保存到您的 `.Renviron` 文件中，以便在未来的会话中使用。

## 2. 使用标记文件

CASSIA使用标记基因数据，通常由差异表达分析工具（如Seurat的 `FindAllMarkers` 函数）生成。

### 2.1 所需格式

CASSIA接受两种标记文件格式：

**1. 原始FindAllMarkers输出（推荐）：**

直接来自Seurat的 `FindAllMarkers` 函数的输出，应包含以下基本列：
- `cluster`：细胞簇标识符
- `gene`：基因名称/符号
- `avg_log2FC`：对数倍数变化
- `p_val_adj`：校正后的p值
- `pct.1`：细胞簇内表达该基因的细胞百分比
- `pct.2`：细胞簇外表达该基因的细胞百分比

```r
# 原始FindAllMarkers输出格式示例
head(markers)
#   p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene
# 1 0     3.02       0.973 0.152 0         0       CD79A
# 2 0     2.74       0.938 0.125 0         0       MS4A1
# 3 0     2.54       0.935 0.138 0         0       CD79B
# ...
```

**2. 处理后格式：**

简化格式，包含细胞簇ID和逗号分隔的基因列表：

```r
# 处理后标记格式示例
head(markers_processed)
#   cluster marker_genes
# 1 0       CD79A,MS4A1,CD79B,HLA-DRA,TCL1A,HLA-DRB1,HLA-DQB1,HLA-DQA1,...
# 2 1       IL7R,CCR7,LEF1,TCF7,FHIT,MAL,NOSIP,CMTM8,TRABD2A,...
# ...
```

### 2.2 示例数据

在本教程中，我们将使用CASSIA自带的示例数据，其中包含来自大肠数据集的六个不同细胞群体的细胞簇：
1. 单核细胞（原始注释不准确，该细胞簇应为施旺细胞。更多证据可在论文中找到）
2. 浆细胞
3. CD8阳性αβT细胞
4. 大肠过渡扩增细胞
5. 肠道肠内分泌细胞
6. 肠隐窝干细胞

```r
# 加载两种格式的示例标记数据
markers_unprocessed <- loadExampleMarkers(processed = FALSE)  # 直接的Seurat FindAllMarkers输出
markers_processed <- loadExampleMarkers(processed = TRUE)     # 处理后格式

# 预览两种数据格式
head(markers_unprocessed)
head(markers_processed)
```

## 3. 运行基础注释

### 3.1 快速模式（一体化）

如需快速、全面的分析，可一次性执行所有步骤：

```r
# 以快速模式运行完整的CASSIA流程
fast_results <- runCASSIA_pipeline(
    output_file_name = "CASSIA_Results",
    tissue = "large intestine",
    species = "human",
    marker = markers_unprocessed
)
```

有关每个参数的详细信息，请参阅[***快速模式文档***](/zh/docs/fast-mode)。

以下是推荐的成本效益分析默认设置：

```r
fast_results <- runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker = markers_unprocessed,
    max_workers = 6,
    annotation_model = "google/gemini-2.5-flash-preview",
    annotation_provider = "openrouter",
    score_model = "deepseek/deepseek-chat-v3-0324",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model = "google/gemini-2.5-flash-preview",
    annotationboost_provider = "openrouter",
    merge_model = "deepseek/deepseek-chat-v3-0324"
    max_retries = 2
)
```

### 3.2 详细批量分析

如需对注释过程有更多控制：

```r

output_name="CASSIA_analysis"

# 使用OpenRouter运行批量分析
batch_results <- runCASSIA_batch(
    marker = markers_unprocessed,
    output_name = output_name
    tissue = "large intestine",
    species = "human"
)

```

### 3.3 质量评分

评估注释的质量：

```r
# 运行质量评分
quality_scores <- runCASSIA_score_batch(
  input_file = paste0(output_name, "_full.csv"),
  output_file = paste0(output_name, "_scored.csv")
)
```

生成的HTML报告页面如下所示，您可以点击按钮导航到相应的细胞簇。

![CASSIA评分报告](/images/report_score.png)


## 4. 解读结果

HTML报告提供可视化摘要，包括：
- 各细胞簇的质量评分和推理过程
- 每个细胞簇的详细注释信息

对于低分数（<75）的细胞簇，请考虑：
1. 检查组织类型是否正确指定
2. 验证物种是否正确
3. 通过基本QC指标、双细胞检测和环境RNA去除来检查细胞簇质量
4. 使用下一节描述的其他分析方法

## 5. 后续步骤

完成基础注释后，您可以探索CASSIA的其他功能。
