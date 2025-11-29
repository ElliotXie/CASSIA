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

**注意：您必须至少设置一个API密钥才能使用CASSIA。**

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
