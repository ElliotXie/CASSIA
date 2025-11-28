---
title: "CASSIA入门：基础分析"
---

本教程介绍使用CASSIA进行单细胞RNA测序分析的基本步骤。我们将重点介绍基础细胞类型注释的核心工作流程。

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

CASSIA依赖Python进行部分后端操作。当您加载软件包时，环境会自动设置。如果遇到问题：

```r
# 手动设置Python环境
setup_cassia_env(conda_env = "cassia_env")
```

### 1.3 API密钥配置

CASSIA至少需要一个API密钥才能运行。您可以使用OpenAI、Anthropic或OpenRouter：

```r
# 设置至少一个API密钥（选择您拥有账户的提供商）
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)
# 或
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)
# 或
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

设置 `persist = TRUE` 会将密钥保存到您的 `.Renviron` 文件中，以便在未来的会话中使用。

## 2. 准备数据

CASSIA使用标记基因数据，通常由Seurat中的 `FindAllMarkers` 函数生成。

### 2.1 示例数据

在本教程中，我们将使用CASSIA自带的示例数据：

```r
# 加载示例标记数据（来自Seurat的FindAllMarkers输出）
markers <- loadExampleMarkers(processed = FALSE)

# 预览数据
head(markers)
```

### 2.2 使用您自己的数据

如果您有自己的Seurat对象，以下是准备数据的方法：

```r
# 从Seurat对象准备数据的示例
# seurat_obj <- your_seurat_object
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj)
# seurat_obj <- FindNeighbors(seurat_obj)
# seurat_obj <- FindClusters(seurat_obj)
# markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
```

## 3. 基础分析工作流程

### 3.1 快速分析（一体化）

如需快速执行所有步骤（注释、评分和注释增强）：

```r
# 以快速模式运行完整的CASSIA流程
results <- runCASSIA_pipeline(
    output_file_name = "MyAnalysis",
    tissue = "blood", # 指定样本的组织类型
    species = "human", # 指定"human"或"mouse"
    marker = markers,
    max_workers = 4, # 根据您计算机的性能调整
    annotation_model = "gpt-4o", # 用于初始注释的模型
    annotation_provider = "openai",
    score_model = "claude-3-5-sonnet-20241022", # 用于质量评分的模型
    score_provider = "anthropic",
    score_threshold = 75, # 高质量注释的最低分数
    annotationboost_model = "claude-3-5-sonnet-20241022", # 用于注释增强的模型
    annotationboost_provider = "anthropic"
)
```

### 3.2 分步分析

如果您希望有更多控制，可以分别运行每个步骤：

```r
# 步骤1：初始注释
batch_results <- runCASSIA_batch(
    marker = markers,
    output_name = "StepByStep",
    model = "gpt-4o",
    tissue = "blood",
    species = "human",
    max_workers = 4,
    n_genes = 50, # 使用的顶级标记基因数量
    provider = "openai"
)

# 步骤2：评分注释
quality_scores <- runCASSIA_score_batch(
    input_file = "StepByStep_full.csv",
    output_file = "StepByStep_scored.csv",
    max_workers = 4,
    model = "claude-3-5-sonnet-20241022",
    provider = "anthropic"
)

# 步骤3：生成报告
runCASSIA_generate_score_report(
    csv_path = "StepByStep_scored.csv",
    output_name = "StepByStep_report.html"
)
```

## 4. 解读结果

CASSIA生成多个输出文件：

### 4.1 完整结果CSV

`_full.csv` 文件包含每个细胞簇的详细注释信息：

```r
# 读取完整结果文件
results <- read.csv("MyAnalysis_full.csv")
head(results)
```

主要列包括：
- `cluster`：细胞簇标识符
- `celltype_1`、`celltype_2`、`celltype_3`：预测的前三个细胞类型
- `reasoning`：模型对注释的解释
- `confidence`：模型对注释的置信度

### 4.2 评分结果CSV

`_scored.csv` 文件包含每个注释的质量评分：

```r
# 读取评分结果文件
scored <- read.csv("MyAnalysis_scored.csv")
head(scored)
```

主要附加列：
- `score`：0-100的质量评分
- `score_reasoning`：评分说明
- `score_category`：基于评分阈值的分类

### 4.3 HTML报告

HTML报告提供注释的可视化摘要，包括：
- 各细胞簇的评分分布
- 每个细胞簇的详细注释信息
- 质量评估和建议

## 5. 常见问题及解决方案

### 5.1 API速率限制

如果遇到速率限制错误：

```r
# 调整最大工作进程数以避免速率限制
results <- runCASSIA_pipeline(
    # ... 其他参数 ...
    max_workers = 2, # 减少以避免速率限制
    # ... 其他参数 ...
)
```

### 5.2 内存问题

对于大型数据集：

```r
# 减少每个细胞簇处理的基因数量
results <- runCASSIA_pipeline(
    # ... 其他参数 ...
    n_genes = 30, # 从默认的50减少
    # ... 其他参数 ...
)
```

## 6. 后续步骤

完成基础分析后，可以考虑：

- 使用 `runCASSIA_annotationboost()` 处理低置信度分数的细胞簇
- 尝试 `compareCelltypes()` 解决模糊的注释
- 实施 `runCASSIA_batch_n_times()` 进行不确定性量化
- 探索 `runCASSIA_subclusters()` 对特定细胞群体进行详细分析

有关这些高级技术的详细信息，请参阅"完整CASSIA工作流程"教程。
