---
title: 单簇分析
---

runCASSIA函数分析单个标记基因簇以识别细胞类型。
**CASSIA可以一次性处理多个簇，因此此方法只适用于分析单个感兴趣的簇。**

### 示例

有关模型设置和推荐的详细信息，请参阅 **[如何选择模型和提供商](setting-up-cassia.md#how-to-select-models-and-providers)** 部分。

#### 示例代码

```R
# 参数
model <- "anthropic/claude-4.5-sonnet"  # 要使用的模型
temperature <- 0
marker_list <- c("CD3D", "CD3E", "CD2", "TRAC")
tissue <- "blood"
species <- "human"
additional_info <- NULL
provider <- "openrouter"  # 或"openai"、"anthropic"

# 运行分析
result <- runCASSIA(
  model = model,
  temperature = temperature,
  marker_list = marker_list,
  tissue = tissue,
  species = species,
  additional_info = additional_info,
  provider = provider,
  validator_involvement = "v1",
  reasoning = "medium"  # 可选: "high", "medium", "low" 用于兼容模型
)

# 查看结构化输出
print(result$structured_output)

# 查看对话历史
print(result$conversation_history)
```

### 参数详情

- **`model`**: 用于分析的大语言模型。选项请参阅 [设置 CASSIA](setting-up-cassia.md)。
- **`temperature`**: 控制模型输出的随机性（0 = 确定性，1 = 创造性）。默认为 0。
- **`marker_list`**: 单个簇的标记基因名称字符向量。
- **`tissue`**: 样本的组织来源。
- **`species`**: 样本的物种（例如 "human", "mouse"）。
- **`additional_info`**: (可选) 关于实验或样本的任何额外上下文信息。
- **`provider`**: 要使用的 API 提供商 ("openrouter", "openai", "anthropic")。
- **`validator_involvement`**: 验证严格程度级别（"v1" 为中等，"v0" 为高）。
- **`reasoning`**: （可选）控制兼容模型的推理深度（"high"、"medium"、"low"）。省略则使用标准模式。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。

更多参数详情请参阅 [批量处理参数详情](batch-processing.md#参数详情)。

_注意：_ 使用OpenRouter时，需指定完整的模型ID。
