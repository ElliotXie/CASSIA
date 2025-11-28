---
title: 单簇分析
---

runCASSIA函数分析单个标记基因簇以识别细胞类型。
**CASSIA可以一次性处理多个簇，因此此方法只适用于分析单个感兴趣的簇。**

### 示例

如果您使用OpenRouter作为提供商，可以指定模型如`"openai/gpt-4o-2024-11-20"`或`"anthropic/claude-3.5-sonnet"`。以下是一些模型推荐：

- **Claude 3.5 Sonnet**（最佳性能，略贵）
    - 模型ID：`"anthropic/claude-3.5-sonnet"`
- **GPT-4o**（平衡选项）
    - 模型ID：`"openai/gpt-4o-2024-11-20"`
- **Llama 3.2**（开源，经济实惠）
    - 模型ID：`"meta-llama/llama-3.2-90b-vision-instruct"`
- **Deepseek v3**（开源，几乎免费，性能与gpt4o相当，最推荐选项）
    - 模型ID：`"deepseek/deepseek-chat-v3-0324"`
    - 模型ID：`"deepseek/deepseek-chat-v3-0324:free"`

#### 示例代码

```R
# 参数
model <- "openai/gpt-4o-2024-11-20"  # 使用OpenRouter时的模型ID
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
  provider = provider
)

# 查看结构化输出
print(result$structured_output)

# 查看对话历史
print(result$conversation_history)
```

_注意：_ 使用OpenRouter时，需指定完整的模型ID。
