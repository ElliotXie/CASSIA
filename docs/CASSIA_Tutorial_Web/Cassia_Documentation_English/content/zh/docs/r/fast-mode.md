---
title: 快速模式
---

## 概述

CASSIA的快速模式提供了一个简化的、单行的解决方案，用于运行完整的分析流程。此模式在单个函数调用中结合了注释、评分和注释增强（用于修正低质量注释），使用优化的默认参数。

## 快速开始

```R
runCASSIA_pipeline(
    output_file_name = "my_analysis",
    tissue = "brain",
    species = "human",
    marker = marker_data,
    max_workers = 4
)
```

## 输入

| 输入 | 类型 | 描述 |
|------|------|------|
| `marker` | 数据框或路径 | 包含cluster和gene列的标记基因数据 |
| `tissue` | 字符串 | 样本的组织类型（如 "brain"、"lung"） |
| `species` | 字符串 | 样本的物种（如 "human"、"mouse"） |

## 参数

### 必需参数

| 参数 | 类型 | 描述 |
|------|------|------|
| `output_file_name` | 字符串 | 输出文件夹和文件的基本名称 |
| `tissue` | 字符串 | 样本的组织类型 |
| `species` | 字符串 | 样本的物种 |
| `marker` | 数据框/路径 | 标记基因数据 |

### 可选参数

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `max_workers` | 6 | 并行处理的工作进程数 |
| `overall_provider` | "openrouter" | 所有流程阶段的主提供商（"openai"、"anthropic"、"openrouter"） |
| `score_threshold` | 75 | 低于此分数（0-100）的注释将触发注释增强 |
| `additional_info` | NULL | 可选的实验上下文（如 "药物X处理"） |
| `validator_involvement` | "v1" | 验证严格程度（"v1" = 中等，"v0" = 高） |
| `merge_annotations` | TRUE | 如果为TRUE，将详细细胞类型合并为更广泛的类别 |
| `annotation_model` | "anthropic/claude-sonnet-4.5" | 用于初始细胞类型注释的模型 |
| `annotation_provider` | "openrouter" | 注释模型的提供商 |
| `score_model` | "openai/gpt-5.1" | 用于质量评分的模型 |
| `score_provider` | "openrouter" | 评分模型的提供商 |
| `annotationboost_model` | "anthropic/claude-sonnet-4.5" | 用于优化低置信度注释的模型 |
| `annotationboost_provider` | "openrouter" | 注释增强模型的提供商 |
| `merge_model` | "google/gemini-2.5-flash" | 用于合并步骤的模型 |
| `merge_provider` | "openrouter" | 合并模型的提供商 |
| `overall_reasoning` | NULL | 所有阶段的推理深度（"low"、"medium"、"high"） |
| `annotation_reasoning` | NULL | 仅覆盖注释阶段的推理级别 |
| `score_reasoning` | NULL | 仅覆盖评分阶段的推理级别 |
| `annotationboost_reasoning` | NULL | 仅覆盖注释增强阶段的推理级别 |
| `merge_reasoning` | NULL | 仅覆盖合并阶段的推理级别 |

## 输出

流程生成一个带时间戳的主文件夹（`CASSIA_Pipeline_{tissue}_{species}_{timestamp}`），包含三个组织有序的子文件夹：

### 文件夹结构

```
CASSIA_Pipeline_brain_human_20240115_143022/
├── 01_annotation_report/
│   └── {name}_report.html          # 交互式HTML报告
├── 02_annotation_boost/
│   └── {cluster_name}/             # 每个低分cluster一个文件夹
│       └── {name}_{cluster}_boosted_report.html
└── 03_csv_files/
    ├── {name}_summary.csv          # 初始注释结果
    ├── {name}_conversations.json   # 完整对话历史
    ├── {name}_scored.csv           # 带质量评分的结果
    ├── {name}_merged.csv           # 合并后的注释（如果启用）
    └── {name}_FINAL_RESULTS.csv    # 合并的最终结果
```

### 输出文件

| 文件夹 | 文件 | 描述 |
|--------|------|------|
| `01_annotation_report` | `{name}_report.html` | 包含所有注释的交互式HTML报告 |
| `02_annotation_boost` | 每个cluster的文件夹 | 低于阈值分数的cluster的增强分析 |
| `03_csv_files` | `{name}_FINAL_RESULTS.csv` | **主要输出** - 包含评分和合并注释的综合结果 |
| `03_csv_files` | `{name}_summary.csv` | 初始细胞类型注释 |
| `03_csv_files` | `{name}_scored.csv` | 带质量评分的注释 |
| `03_csv_files` | `{name}_merged.csv` | 更广泛的类别分组（如果 `do_merge_annotations = TRUE`） |
| `03_csv_files` | `{name}_conversations.json` | 完整的LLM对话历史，用于可重复性 |

### 将CASSIA结果添加到Seurat对象

```R
seurat_corrected <- add_cassia_to_seurat(
  seurat_obj = seurat_corrected, # 您想要添加CASSIA结果的Seurat对象
  cassia_results_path = "/FastAnalysisResults_scored.csv", # 评分结果保存的位置，指定路径
  cluster_col = "celltype", # Seurat对象中包含细胞类型的列
  cassia_cluster_col="True Cell Type" # 评分结果中包含真实细胞类型的列
)

# 这将向Seurat对象添加六个新列：一般细胞类型、所有三个亚细胞类型、最可能的细胞类型、第二可能的细胞类型、第三可能的细胞类型和混合细胞类型，以及每种细胞类型的质量评分。
```

### 性能提示

- 为获得最佳性能，根据系统的CPU核心调整 `max_workers`
- 使用 `additional_info` 提供相关的实验上下文
- 监控 `score_threshold` 以平衡严格性和处理量

---

接下来我们详细介绍每个功能...
