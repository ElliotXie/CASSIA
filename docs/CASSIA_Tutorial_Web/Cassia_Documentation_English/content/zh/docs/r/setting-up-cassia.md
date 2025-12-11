---
title: 设置CASSIA
---


首先，确保您已安装reticulate包和devtools包。

```r
install.packages("reticulate")
install.packages("devtools")
```

接下来，您需要安装CASSIA包。

```r
# 安装CASSIA包
library(devtools)
devtools::install_github("ElliotXie/CASSIA/CASSIA_R")
```

## 设置Python环境

CASSIA依赖Python进行一些后端处理。当您加载CASSIA包时，它会尝试自动设置所需的Python环境。但是，如果您遇到问题，可以使用`setup_cassia_env()`函数自动创建和配置必要的Python环境。

```r
library(CASSIA)

# 如有需要，自动设置Python环境
setup_cassia_env(conda_env = "cassia_env")
```

此函数将：

- 创建一个名为`cassia_env`的新Python虚拟环境（默认）。如果创建虚拟环境失败，它将尝试创建一个Conda环境。
- 安装所需的Python包：`openai`、`pandas`、`numpy`、`requests`、`anthropic`、`matplotlib`和`seaborn`。

## 设置API密钥

要使用如OpenAI的GPT-4、Anthropic的Claude或通过OpenRouter访问的模型，您需要先从提供商获取API密钥，再通过setLLMApiKey()函数设置API密钥。从提供商获取API密钥大约需要3分钟。（推荐优先设置Openrouter API密钥）

**注意：您必须至少设置一个API密钥才能使用CASSIA。** 您也可以使用**[自定义 API 提供商](#自定义-api-提供商)**，如 DeepSeek 或本地 LLM。

```r
# 对于OpenAI
setLLMApiKey("your_openai_api_key", provider = "openai", persist = TRUE)

# 对于Anthropic
setLLMApiKey("your_anthropic_api_key", provider = "anthropic", persist = TRUE)

# 对于OpenRouter
setLLMApiKey("your_openrouter_api_key", provider = "openrouter", persist = TRUE)
```

- 将`"your_api_key"`替换为您的实际API密钥。
- 根据您的提供商将`provider`设置为`"openai"`、`"anthropic"`或`"openrouter"`。
- 设置`persist = TRUE`将密钥保存在您的`.Renviron`文件中，供future会话使用。

## 验证 API 密钥

您可以在运行分析之前验证 API 密钥是否正常工作：

```r
# 验证所有已配置的提供商
validate_api_keys(force_revalidate = TRUE)

# 验证特定提供商
validate_api_keys("openai", force_revalidate = TRUE)
```

## 如何选择模型和提供商

有三个提供商可供选择：`openrouter`、`openai`和`anthropic`。每个提供商都有自己的模型和定价。
请注意，模型名称必须完全按照下面所示设置，否则将找不到模型。

### OpenRouter

OpenRouter是一个平台，提供对主要提供商支持的几乎所有模型的访问。建议使用OpenRouter，因为它具有更高的速率限制，并且可以访问包括开源选项在内的各种模型。

- `anthropic/claude-sonnet-4.5`: 当前默认的高性能模型。
- `openai/gpt-5.1`: 平衡选项。
- `google/gemini-2.5-flash`: 性能最好的低成本模型之一。

### OpenAI

- `gpt-5.1`: 高性能模型。

### Anthropic

- `claude-sonnet-4.5`: 高性能模型。

### 其他提供商

这些模型可以通过其自有 API 使用。设置方法请参阅 **[自定义 API 提供商](#自定义-api-提供商)**。

- `deepseek-chat` (DeepSeek v3.2): 高性能，价格实惠。提供商：`https://api.deepseek.com`
- `glm-4.6` (GLM 4.6): 快速且经济实惠。提供商：`https://api.z.ai/api/paas/v4/`
- `kimi-k2` (Kimi K2): 强大的推理能力。提供商：`https://api.moonshot.cn/v1`

### 本地 LLM

- `gpt-oss:20b`: 可通过 Ollama 在本地运行。适合大批量分析，准确率可接受。设置方法请参阅 **[本地 LLM](#本地-llmollama-lm-studio)**。

## 智能模型设置 (Smart Model Settings)

CASSIA包含一个智能模型选择系统，允许您使用简单的别名或"层级"（tiers）而不是记住确切的模型版本字符串。

### 层级快捷方式 (Tier Shortcuts)
您可以将这些快捷方式与任何提供商一起使用，以获得适合您需求的模型：

- `"best"`: 选择性能最高的模型（例如，`gpt-5.1`，`claude-opus-4.5`）
- `"balanced"`: 选择性能和成本平衡良好的模型（例如，`gpt-4o`，`claude-sonnet-4.5`）
- `"fast"`: 选择最快/最便宜的模型（例如，`gpt-5-mini`，`claude-haiku-4.5`）

示例：
```r
# 这将自动选择OpenAI的最佳模型 (gpt-5.1)
runCASSIA_pipeline(..., model = "best", provider = "openai")
```

### 模糊匹配和别名 (Fuzzy Matching & Aliases)
您还可以使用通用名称，CASSIA会将它们解析为正确的版本：

- `"gpt"` -> 解析为 `gpt-5.1` (对于OpenAI)
- `"claude"` -> 解析为 `claude-sonnet-4.5` (对于Anthropic)
- `"gemini"` -> 解析为 `google/gemini-2.5-flash` (对于OpenRouter)

这使您的代码对模型版本更新更加稳健。

## 自定义 API 提供商

CASSIA 支持任何兼容 OpenAI 的 API 端点，让您可以使用自定义提供商，如 DeepSeek、本地 LLM 服务器或其他第三方服务。

> **中国大陆用户推荐**：由于 Claude 和 GPT 等模型在中国大陆可能存在访问限制，我们推荐使用 **DeepSeek**。DeepSeek 是中国公司开发的高性能大语言模型，性能与 GPT-4o 相当，价格实惠，访问稳定。

### 设置自定义提供商

使用自定义 API 提供商时，将完整的基础 URL 作为 `provider` 参数：

```r
# 为自定义提供商设置 API 密钥
setLLMApiKey("your-api-key", provider = "https://api.your-provider.com", persist = TRUE)

# 在分析中使用
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    provider = "https://api.your-provider.com",
    model = "your-model-name",
    tissue = "brain",
    species = "human"
)
```

### DeepSeek 使用示例（推荐）

DeepSeek 提供高性能模型，价格实惠，特别适合中国用户：

1. 从 [DeepSeek 开放平台](https://platform.deepseek.com/) 获取 API 密钥
2. 在 CASSIA 中配置：

```r
setLLMApiKey("your-deepseek-key", provider = "https://api.deepseek.com", persist = TRUE)

runCASSIA_pipeline(
    output_file_name = "analysis",
    marker = markers,
    annotation_provider = "https://api.deepseek.com",
    annotation_model = "deepseek-chat",
    tissue = "brain",
    species = "human"
)
```

### 本地 LLM（Ollama、LM Studio）

为了完全的数据隐私和零 API 费用，您可以在本地运行 LLM。CASSIA 支持任何兼容 OpenAI 的本地服务器。

**本地 URL 无需 API 密钥。**

#### Ollama 设置

1. 从 [ollama.ai](https://ollama.ai) 安装 Ollama
2. 拉取模型：`ollama pull gpt-oss:20b`
3. Ollama 自动运行在 `http://localhost:11434`

#### 使用方法

```r
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    provider = "http://localhost:11434/v1",
    model = "gpt-oss:20b",
    tissue = "brain",
    species = "human"
)
```

## 推理深度参数

**注意：** 此参数仅对 OpenAI GPT-5 系列模型（如 `gpt-5.1`）有效。推荐通过 OpenRouter 使用，或作为已验证的 OpenAI 用户使用。

`reasoning` 参数控制兼容模型的推理深度。

### 语法

```r
reasoning = "high"    # 最大推理深度
reasoning = "medium"  # 平衡（推荐用于 GPT-5.1）
reasoning = "low"     # 最小推理
reasoning = NULL      # 标准模式（或省略该参数）
```

### 提供商说明

- **OpenAI 通过 OpenRouter**: 完全控制推理参数。**推荐**使用以避免直接 OpenAI API 所需的身份验证。
- **OpenAI 直接访问**: 需要身份验证才能使用推理模型。
- **Anthropic Claude**: 默认自动使用最高推理深度。
- **Gemini**: 动态思考 - 模型自行决定何时以及思考多少。

### 建议

- **GPT-5.1**: 使用 `reasoning = "medium"` - 最高深度可能需要很长时间
- **GPT-4o**: 无需推理参数仍有出色性能
- **Claude**: 无需设置 - 自动使用最佳推理
- **一般**: 更高深度 = 更长处理时间 + 更高成本

示例：
```r
# 通过 OpenRouter 使用 GPT-5.1 的推理功能
runCASSIA_batch(
    marker = markers,
    output_name = "results",
    model = "openai/gpt-5.1",
    provider = "openrouter",
    reasoning = "medium",  # 推荐用于平衡速度/质量
    tissue = "brain",
    species = "human"
)
```
