---
title: "使用标记文件进行基础注释"
---

本教程介绍使用CASSIA进行细胞类型注释的基本步骤，适用于您已经准备好标记基因列表的情况。这非常适合您已经完成聚类分析并希望对细胞簇进行注释的场景。

## 1. 安装和设置

在开始之前，请确保已安装并配置CASSIA。有关详细说明，请参阅[**设置CASSIA**](/zh/docs/r/setting-up-cassia)文档。

```r
library(CASSIA)

# 设置API密钥（推荐使用OpenRouter）
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

## 2. 使用标记文件

CASSIA使用差异表达分析产生的标记基因数据。它接受Seurat的`FindAllMarkers`输出、Scanpy的`rank_genes_groups`输出或简化格式。有关详细的格式说明，请参阅[**批量处理**](/zh/docs/r/batch-processing#标记数据格式)文档。

### 2.1 示例数据

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

有关每个参数的详细信息，请参阅[**快速模式**](/zh/docs/fast-mode)文档。

### 3.2 批量分析（更快）

如需更快的注释，不运行质量评分、合并和注释增强。这对于大多数情况已经足够：

```r
output_name="CASSIA_analysis"

# 使用OpenRouter运行批量分析
batch_results <- runCASSIA_batch(
    marker = markers_unprocessed,
    output_name = output_name,
    tissue = "large intestine",
    species = "human",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

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
