---
title: 安装 CASSIA
---

首先，使用 pip 安装 CASSIA 包：

```bash
pip install CASSIA
```

## 导入 CASSIA

```python
import CASSIA
```

## 设置 API 密钥

要使用 OpenAI 的 GPT-4、Anthropic 的 Claude 或通过 OpenRouter 使用模型，您首先需要从提供商处获取 API 密钥，然后使用 `set_api_key` 函数设置您的 API 密钥。

**注意：您必须设置至少一个 API 密钥才能使用 CASSIA。**

**您只需选择一个提供商。** 推荐使用 OpenRouter，因为它提供多种模型的访问。

```python
# 设置 API 密钥（选择一个提供商）
CASSIA.set_api_key("your-openrouter-key", provider="openrouter")  # 推荐
# CASSIA.set_api_key("your-openai-key", provider="openai")
# CASSIA.set_api_key("your-anthropic-key", provider="anthropic")
```

- 将 `"your-key"` 替换为您的实际 API 密钥。
- 根据您的提供商，将 `provider` 设置为 `"openai"`、`"anthropic"` 或 `"openrouter"`。

## 验证 API 密钥

您可以在运行分析之前验证 API 密钥是否正常工作：

```python
# 验证所有已配置的提供商
CASSIA.validate_api_keys(force_revalidate=True)

# 验证特定提供商
CASSIA.validate_api_keys("openai", force_revalidate=True)
```

## 如何选择模型和提供商

有三个提供商可供选择：`openrouter`、`openai` 和 `anthropic`。每个提供商都有自己的模型和定价。
**注意，模型名称必须完全按照下面所示设置，否则将找不到模型。**

### OpenRouter

OpenRouter 是一个平台，提供对主要提供商支持的几乎所有模型的访问。建议使用 OpenRouter，因为它具有更高的速率限制，并且可以访问包括开源选项在内的各种模型。

- `anthropic/claude-sonnet-4.5`: 当前默认的高性能模型。
- `openai/gpt-5.1`: 平衡的选项。
- `google/gemini-2.5-flash`: 性能最好的低成本模型之一。

### OpenAI

- `gpt-5.1`: 高性能模型。

### Anthropic

- `claude-sonnet-4.5`: 高性能模型。

## 智能模型设置（推荐）

CASSIA 包含一个智能模型选择系统，允许您使用简单的别名或“层级”代替记住确切的模型版本字符串。这使您的代码对模型版本更新更具鲁棒性。

### 层级快捷方式
您可以将这些快捷方式与任何提供商一起使用，以获取适合您需求的模型：

- `"best"`: 选择性能最高的模型（例如，`gpt-5.1`，`claude-opus-4.5`）
- `"balanced"`: 选择性能和成本平衡良好的模型（例如，`gpt-4o`，`claude-sonnet-4.5`）
- `"fast"`: 选择最快/最便宜的模型（例如，`gpt-5-mini`，`claude-haiku-4.5`）

示例：
```python
# 这将自动为 OpenAI 选择最佳模型
CASSIA.runCASSIA_pipeline(..., model = "best", provider = "openai")
```

### 模糊匹配和别名
您也可以使用常用名称，CASSIA 将把它们解析为正确的版本：

- `"gpt"` -> 解析为 `gpt-5.1` (用于 OpenAI)
- `"claude"` -> 解析为 `claude-sonnet-4.5` (用于 Anthropic)
- `"gemini"` -> 解析为 `google/gemini-2.5-flash` (用于 OpenRouter)

## 推理深度参数

`reasoning` 参数控制兼容模型的推理深度。

### 语法

```python
# 简单字符串格式
reasoning="high"    # 最大推理深度
reasoning="medium"  # 平衡（推荐用于 GPT-5.1）
reasoning="low"     # 最小推理

# 字典格式（仅限 Python）
reasoning={"effort": "high"}

# 标准模式（无扩展推理）
reasoning=None  # 或省略该参数
```

### 提供商说明

- **OpenAI 通过 OpenRouter**: 完全控制推理参数。**推荐**使用以避免直接 OpenAI API 所需的身份验证。
- **OpenAI 直接访问**: 需要身份验证才能使用推理模型。
- **Anthropic Claude**: 默认自动使用最高推理深度。
- **Gemini**: 动态思考 - 模型自行决定何时以及思考多少。

### 建议

- **GPT-5.1**: 使用 `reasoning="medium"` - 最高深度可能需要很长时间
- **GPT-4o**: 无需推理参数仍有出色性能
- **Claude**: 无需设置 - 自动使用最佳推理
- **一般**: 更高深度 = 更长处理时间 + 更高成本

示例：
```python
# 通过 OpenRouter 使用 GPT-5.1 的推理功能
CASSIA.runCASSIA_batch(
    marker=markers,
    output_name="results",
    model="openai/gpt-5.1",
    provider="openrouter",
    reasoning="medium",  # 推荐用于平衡速度/质量
    tissue="brain",
    species="human"
)
```

## 自定义 API 提供商

CASSIA 支持任何兼容 OpenAI 的 API 端点，让您可以使用自定义提供商，如 DeepSeek、本地 LLM 服务器或其他第三方服务。

> **中国大陆用户推荐**：由于 Claude 和 GPT 等模型在中国大陆可能存在访问限制，我们推荐使用 **DeepSeek**。DeepSeek 是中国公司开发的高性能大语言模型，性能与 GPT-4o 相当，价格实惠，访问稳定。

### 设置自定义提供商

使用自定义 API 提供商时，将完整 URL 作为 `provider` 参数：

```python
# 为自定义提供商设置 API 密钥
CASSIA.set_api_key("your-api-key", provider="https://api.deepseek.com")

# 在分析中使用
CASSIA.runCASSIA_batch(
    marker=markers,
    output_name="results",
    provider="https://api.deepseek.com",
    model="deepseek-chat",
    tissue="brain",
    species="human"
)
```

### DeepSeek 使用示例（推荐）

DeepSeek 提供高性能模型，价格实惠，特别适合中国用户：

1. 从 [DeepSeek 开放平台](https://platform.deepseek.com/) 获取 API 密钥
2. 在 CASSIA 中配置：

```python
CASSIA.set_api_key("your-deepseek-key", provider="https://api.deepseek.com")

CASSIA.runCASSIA_pipeline(
    output_file_name="analysis",
    marker=markers,
    annotation_provider="https://api.deepseek.com",
    annotation_model="deepseek-chat",
    tissue="brain",
    species="human"
)
```

### 兼容的提供商

任何遵循 OpenAI 聊天补全格式的 API 都可以使用，包括：
- DeepSeek (`https://api.deepseek.com`) - **推荐中国用户使用**
- 本地 LLM 服务器（如 Ollama、vLLM）
- 其他兼容 OpenAI 的服务
