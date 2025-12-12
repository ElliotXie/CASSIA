---
title: 单簇分析
---

## 概述

`runCASSIA` 函数分析单个标记基因簇以识别细胞类型。此函数专为只需分析单个簇的用户设计。

注意：CASSIA 可以通过[批量处理](batch-processing.md)一次性处理多个簇。当您只需要注释单个簇时，请使用此函数。

---

## 快速开始

```R
result <- runCASSIA(
    marker_list = c("CD3D", "CD3E", "CD2", "TRAC"),
    model = "anthropic/claude-sonnet-4.5",
    tissue = "blood",
    species = "human",
    provider = "openrouter"
)

# 查看注释结果
print(result$structured_output)
```

模型推荐请参阅[如何选择模型和提供商](setting-up-cassia.md#how-to-select-models-and-providers)。

---

## 输入

### 标记基因列表格式

提供一个包含簇标记基因名称的字符向量：

```R
marker_list <- c("CD3D", "CD3E", "CD2", "TRAC", "IL7R")
```

这些应该是表征您感兴趣簇的顶级差异表达基因。

---

## 参数

### 必需参数

| 参数 | 描述 |
|------|------|
| `marker_list` | 簇的标记基因名称字符向量 |
| `model` | LLM 模型 ID（例如 `"anthropic/claude-sonnet-4.5"`） |
| `tissue` | 组织类型（例如 `"blood"`、`"brain"`） |
| `species` | 物种（例如 `"human"`、`"mouse"`） |
| `provider` | API 提供商（`"openrouter"`、`"openai"`、`"anthropic"`） |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `temperature` | 0 | 输出随机性（0=确定性，1=创造性）。保持为 0 以获得可重复结果。 |
| `additional_info` | `NULL` | 关于样本的额外实验上下文 |
| `validator_involvement` | `"v1"` | 验证强度：`"v1"`（中等）或 `"v0"`（高，较慢） |
| `reasoning` | `NULL` | 兼容模型的推理深度（`"low"`、`"medium"`、`"high"`）。见下文。 |

### 参数详情

**模型选择**
- 默认：`anthropic/claude-sonnet-4.5` 以获得最佳性能
- 替代：`google/gemini-2.5-flash` 以获得更快分析
- 使用 OpenRouter 时，请指定完整的模型 ID
- 详细推荐请参阅[如何选择模型和提供商](setting-up-cassia.md#how-to-select-models-and-providers)

**推理参数**
- 控制兼容模型的推理深度（通过 OpenRouter 使用 GPT-5 系列）
- 选项：`"low"`、`"medium"`、`"high"`
- 省略此参数以使用标准模式
- 详见[推理深度参数](setting-up-cassia.md#推理深度参数)

**额外上下文**
- 使用 `additional_info` 提供实验上下文
- 示例：`"来自肿瘤微环境的样本，关注免疫浸润"`

---

## 输出

函数返回包含两个组件的列表：

| 组件 | 描述 |
|------|------|
| `structured_output` | 包含预测细胞类型和推理的注释结果 |
| `conversation_history` | 用于调试和透明度的完整对话日志 |

