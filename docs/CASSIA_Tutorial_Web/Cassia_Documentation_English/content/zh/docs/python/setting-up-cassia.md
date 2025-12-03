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

