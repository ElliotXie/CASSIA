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

- 如果不存在，创建一个名为`cassia_env`的新Conda环境。
- 安装所需的Python包：`openai`、`pandas`、`numpy`、`scikit-learn`、`requests`和`anthropic`。

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

有三个提供商可供选择：`openai`、`anthropic`和`openrouter`。每个提供商都有自己的模型和定价。
请注意，模型名称必须完全按照下面所示设置，否则将找不到模型。

### OpenAI

- `gpt-4o`

`gpt-4o`是最均衡的模型。在与GPTcelltype对比基准测试CASSIA的性能时，我们将其作为默认选择。

### Anthropic

- `claude-3-5-sonnet-20241022`

`claude-3-5-sonnet-20241022`是最强大的模型。在基准测试中，我们将其用于评分和注释增强。

### OpenRouter

OpenRouter是一个平台，提供对主要提供商支持的几乎所有模型的访问。实际上，建议使用OpenRouter访问claude-3-5-sonnet，因为它具有最高的速率限制。我们还可以使用它访问许多开源模型，如llama-3.2和DeepseekV3。这些开源模型价格更便宜，性能略有下降。

- `anthropic/claude-3.5-sonnet`
- `openai/gpt-4o-2024-11-20`
- `meta-llama/llama-3.2-90b-vision-instruct` 
- `deepseek/deepseek-chat-v3-0324`(最推荐的模型)
- `deepseek/deepseek-chat-v3-0324:free`(DeepseekV3的免费版本，稍慢且稳定性较低)

DeepseekV3是最推荐的模型。它是一个几乎与gpt4o一样好的免费开源模型。
