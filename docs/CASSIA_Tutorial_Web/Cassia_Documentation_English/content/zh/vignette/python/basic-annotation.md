---
title: 基本注释
---

本 Python 教程演示了使用 CASSIA 对单细胞 RNA 测序数据进行细胞类型注释的完整工作流程。我们将分析包含六个不同细胞群体的肠道细胞数据集：

1. 单核细胞
2. 浆细胞
3. cd8 阳性 alpha-beta T 细胞
4. 大肠过渡扩增细胞
5. 肠道肠内分泌细胞
6. 肠隐窝干细胞

## 设置和环境准备

首先，让我们安装并导入所需的包：

```bash
pip install CASSIA
```

```python
import CASSIA
```

### 设置 API 密钥

**您只需选择一个提供商。** 推荐使用 OpenRouter，因为它提供多种模型的访问。

```python
# 设置 API 密钥（选择一个提供商）
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")  # 推荐
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
```

### 加载数据

```python
processed_markers = CASSIA.loadmarker(marker_type="processed")
unprocessed_markers = CASSIA.loadmarker(marker_type="unprocessed")
subcluster_results = CASSIA.loadmarker(marker_type="subcluster_results")

# 列出可用的标记集
available_markers = CASSIA.list_available_markers()
print(available_markers) 
```

## 快速模式

以快速模式运行 CASSIA 流程。这是一个一步到位的过程，可以快速获得结果。

```python
# 以快速模式运行 CASSIA 流程
CASSIA.runCASSIA_pipeline(
    output_file_name = "FastAnalysisResults",
    tissue = "large intestine",
    species = "human",
    marker = unprocessed_markers,
    max_workers = 6,  # 与数据集中的聚类数匹配
    annotation_model = "anthropic/claude-sonnet-4.5",
    annotation_provider = "openrouter",
    score_model = "openai/gpt-5.1",
    score_provider = "openrouter",
    score_threshold = 75,
    annotationboost_model="anthropic/claude-sonnet-4.5",
    annotationboost_provider="openrouter",
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter"
)
```

## 详细的批量分析

为了进行更精细的控制，您可以单独运行各个步骤。

```python
output_name="intestine_detailed"

# 运行批量分析
CASSIA.runCASSIA_batch(
    marker = unprocessed_markers,
    output_name = output_name,
    model = "anthropic/claude-sonnet-4.5",
    tissue = "large intestine",
    species = "human",
    max_workers = 6,  # 匹配聚类数
    n_genes = 50,
    additional_info = None,
    provider = "openrouter",
    reasoning = "medium"  # 可选: 用于 GPT-5 系列模型
)
```

> **提示: 推理深度参数**
>
> 使用 `reasoning="medium"` 配合 GPT-5.1 可获得增强推理而不会耗费过长时间。对于 OpenAI 推理模型，我们推荐使用 OpenRouter 以避免身份验证要求。Claude 模型默认使用最佳推理。详见 [推理深度参数](/docs/python/setting-up-cassia#推理深度参数)。

## 质量评分

注释后，运行质量评分以评估置信度。

```python
# 运行质量评分
CASSIA.runCASSIA_score_batch(
    input_file = output_name + "_summary.csv",  # JSON 自动检测
    output_file = output_name + "_scored.csv",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter"
)
```

