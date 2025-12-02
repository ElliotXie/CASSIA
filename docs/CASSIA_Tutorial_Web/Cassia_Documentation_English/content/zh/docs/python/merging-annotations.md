---
title: 合并注释（可选）
---

合并注释将详细的细胞类型注释分组到多个粒度级别的更广泛类别中。这对于需要细胞类型层次分组的可视化和下游分析非常有用。

### 概述

CASSIA 提供两个合并函数：
- `merge_annotations()`: 在单个详细级别合并
- `merge_annotations_all()`: 同时在所有三个级别合并

### 详细级别

| 级别 | 列名 | 描述 | 示例 |
|-------|-------------|-------------|---------
| 广泛 | Merged_Grouping_1 | 一般类别 | "T 细胞"，"髓系细胞" |
| 详细 | Merged_Grouping_2 | 中等特异性 | "CD4 T 细胞"，"巨噬细胞" |
| 非常详细 | Merged_Grouping_3 | 规范化的特定名称 | "CD4+ 初始 T 细胞" |

### 单级别合并

使用 `merge_annotations()` 在特定详细级别进行合并：

```python
from CASSIA import merge_annotations

# 在广泛级别合并
result_df = merge_annotations(
    csv_path="batch_results_full.csv",
    output_path="merged_broad.csv",
    provider="openrouter",
    model="google/gemini-2.5-flash",
    detail_level="broad",  # 选项："broad"、"detailed"、"very_detailed"
    batch_size=10
)

# 检查结果
print(result_df[['Cluster ID', 'main_cell_type', 'Merged_Grouping_1']])
```

### 多级别合并（一次全部）

使用 `merge_annotations_all()` 通过并行处理同时在所有三个级别进行合并：

```python
from CASSIA import merge_annotations_all

# 同时在所有级别合并
result_df = merge_annotations_all(
    csv_path="batch_results_full.csv",
    output_path="merged_all.csv",
    provider="openrouter",
    model="google/gemini-2.5-flash",
    batch_size=10
)

# 结果包含所有三个分组列
print(result_df[['main_cell_type', 'Merged_Grouping_1', 'Merged_Grouping_2', 'Merged_Grouping_3']])
```

### 参数详情

- **`csv_path`**: CASSIA 批处理结果 CSV 文件的路径（`runCASSIA_batch` 的输出）。
- **`output_path`**: 带有合并分组的输出 CSV 文件路径。
- **`provider`**: API 提供商（"openrouter"、"openai"、"anthropic"）。
- **`model`**: 用于分组决策的 LLM 模型。
- **`detail_level`**: 分组级别（"broad"、"detailed"、"very_detailed"）- 仅适用于 `merge_annotations()`。
- **`batch_size`**: 每次 LLM 调用处理的细胞类型数量（默认：10）。

### 输出示例

对于浆细胞注释：

| 详细级别 | 分组 |
|--------------|----------|
| 广泛 | B 细胞 / 浆细胞 |
| 详细 | 浆细胞 |
| 非常详细 | 浆细胞 |

对于 CD8 阳性，alpha-beta T 细胞：

| 详细级别 | 分组 |
|--------------|----------|
| 广泛 | T 细胞 |
| 详细 | CD8 T 细胞 |
| 非常详细 | CD8+ alpha-beta T 细胞 |

### 使用场景

1. **可视化**：在不同粒度级别创建汇总图。
2. **跨数据集比较**：使用规范化的细胞类型名称比较结果。
3. **下游分析**：按广泛类别对细胞进行分组以进行统计分析。

### 与流水线集成

通过 `do_merge_annotations` 参数，合并功能自动可在 `runCASSIA_pipeline()` 中使用：

```python
CASSIA.runCASSIA_pipeline(
    ...,
    do_merge_annotations=True,  # 启用自动合并
    merge_model="google/gemini-2.5-flash",
    merge_provider="openrouter"
)
```

### 注意事项

- 合并过程使用基于 LLM 的分组，因此运行之间的结果可能略有不同。
- 建议使用较低的温度（例如 0.3）以获得更一致的分组。
- `merge_annotations_all()` 中的并行处理默认使用 3 个工作线程的 ThreadPoolExecutor。
