---
title: Symphony Compare (Optional)
---

`symphonyCompare` 智能体充当虚拟专家小组来解决模棱两可的细胞类型注释。它协调多个 AI 模型（智能体“交响乐团”）来比较潜在的细胞类型，在讨论轮次中辩论它们的发现，并根据标记基因证据达成共识。

### 用法

此智能体在运行默认 CASSIA 流程后特别有用，如果您不确定特定聚类的身份。例如，区分浆细胞的不同亚型。

```python
# 这里的标记复制自 CASSIA 之前的结果。
marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

# 运行 Symphony Compare 分析
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
    marker_set = marker,
    species = "human",
    model_preset = "premium",  # 选项: "premium", "budget"
    output_basename = "plasma_cell_comparison",
    enable_discussion = True
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

### 参数详情

- **`tissue`**: 正在分析的组织类型（例如，“大肠”）。
- **`celltypes`**: 要比较的 2-4 种细胞类型的列表。
- **`marker_set`**: 逗号分隔的标记基因字符串。
- **`species`**: 样本的物种（默认：“人类”）。
- **`model_preset`**: 要使用的模型配置。
    - `"premium"` (默认): 高性能组合 (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4)。
    - `"budget"`: 具成本效益的模型 (DeepSeek V3.2, Grok 4 Fast, Kimi K2, Gemini 2.5 Flash)。
- **`output_basename`**: 输出文件的基本名称。
- **`enable_discussion`**: 是否启用模型之间的多轮辩论（默认：`True`）。
- **`max_discussion_rounds`**: 最大讨论轮数（默认：2）。
- **`consensus_threshold`**: 达成共识所需的模型比例（默认：0.8）。

### 输出文件
- `{output_basename}.csv`: 来自所有模型和轮次的详细比较结果、推理和分数。
- `{output_basename}_report.html`: 可视化辩论和共识过程的交互式 HTML 报告。
