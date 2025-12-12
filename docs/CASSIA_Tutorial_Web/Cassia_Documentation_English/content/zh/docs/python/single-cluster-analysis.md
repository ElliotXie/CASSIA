---
title: 单聚类分析
---

`runCASSIA` 函数分析单个标记基因聚类以识别细胞类型。
请注意，CASSIA 设计用于同时处理多个聚类，此函数专为只有一个聚类需要分析的用户设计。

### 示例

有关模型设置和建议的详细信息，请参阅 **[如何选择模型和提供商](setting-up-cassia.md#如何选择模型和提供商)** 部分。

#### 示例代码

```python
import CASSIA

# 参数
model = "anthropic/claude-sonnet-4.5"  # 使用的模型
temperature = 0
marker_list = ["CD3D", "CD3E", "CD2", "TRAC"]
tissue = "blood"
species = "human"
additional_info = None
provider = "openrouter"  # 或 "openai", "anthropic"

# 运行分析
result, conversation_history, _ = CASSIA.runCASSIA(
    model=model,
    temperature=temperature,
    marker_list=marker_list,
    tissue=tissue,
    species=species,
    additional_info=additional_info,
    provider=provider,
    validator_involvement="v1",
    reasoning="medium"  # 可选: "high", "medium", "low" 用于兼容模型
)

# 查看结构化输出
print(result['main_cell_type'])
print(result['sub_cell_types'])

# 查看对话历史
print(conversation_history)
```

### 参数详情

- **`model`**: 用于分析的 LLM 模型。有关选项，请参阅 [安装 CASSIA](setting-up-cassia.md)。
- **`temperature`**: 控制模型输出的随机性（0 = 确定性，1 = 创造性）。默认为 0。
- **`marker_list`**: 单个聚类的标记基因名称列表。
- **`tissue`**: 样本的组织来源。
- **`species`**: 样本的物种（例如，"human"，"mouse"）。
- **`additional_info`**: （可选）关于实验或样本的额外上下文。
- **`provider`**: 使用的 API 提供商（"openrouter"、"openai"、"anthropic"）。
- **`validator_involvement`**: 验证严格程度（"v1" 为中等，"v0" 为高）。
- **`reasoning`**: （可选）控制兼容模型的推理深度（"high"、"medium"、"low"）。Python 还接受字典格式：`{"effort": "high"}`。省略则使用标准模式。详见 [推理深度参数](setting-up-cassia.md#推理深度参数)。

更多参数详情请参阅 [批量处理参数详情](batch-processing.md#参数详情)。

### 返回值

函数返回一个元组：`(result, conversation_history, _)`

- **`result`**: 包含以下内容的字典：
  - `main_cell_type`: 主要细胞类型预测
  - `sub_cell_types`: 可能的亚型列表
- **`conversation_history`**: 与模型的完整对话，用于透明度

_注意:_ 使用 OpenRouter 时，请指定完整的模型 ID。
