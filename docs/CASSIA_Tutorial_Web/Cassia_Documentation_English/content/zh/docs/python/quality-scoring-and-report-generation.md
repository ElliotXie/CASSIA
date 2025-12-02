---
title: 质量评分
---

质量评分有助于评估细胞类型注释的可靠性。CASSIA 通过 `runCASSIA_score_batch` 函数提供自动评分功能，该函数分析每个注释背后的推理和证据。

### 运行质量评分

```python
# 运行质量评分
CASSIA.runCASSIA_score_batch(
    input_file = output_name + "_full.csv",
    output_file = output_name + "_scored.csv",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter"
)
```

### 参数详情

- **`input_file`**: 完整注释结果 CSV 文件的路径（由 `runCASSIA_batch` 生成）。
- **`output_file`**: 输出评分结果 CSV 文件名。
- **`max_workers`**: 并行评分线程数。
- **`model`**: 用于质量评分的 LLM。建议使用 `claude-4.5-sonnet` 或 `gpt-4o` 等高能力模型以获得准确评分。
- **`provider`**: 模型的 API 提供商（例如，"openrouter"）。

### 解读评分

- **90-100**: 高置信度，证据充分。
- **76-89**: 良好置信度，证据充足。
- **<75**: 低置信度。这些聚类是使用注释增强智能体或比较智能体进行进一步分析的候选者。

### 报告生成

从您的分析生成详细报告。此步骤通常在质量评分之后进行。评分报告包括 CASSIA 的所有输出，包括结构化输出、对话历史记录和质量评分。

```python
# 生成质量报告
CASSIA.runCASSIA_generate_score_report(
    csv_path = output_name + "_scored.csv",
    index_name = output_name + "_report.html"
)
```

### 参数详情

- **`csv_path`**: 评分结果 CSV 文件的路径（例如，"my_annotation_scored.csv"）。
- **`index_name`**: 生成的报告索引文件的名称（例如，"my_report.html"）。

_从评分结果生成单独的报告和索引页面。_

### 替代报告生成函数

CASSIA 还提供独立的报告生成函数，可以从任何批处理结果创建 HTML 报告：

#### 从 CSV 文件生成

```python
from CASSIA.reports.generate_batch_report import generate_batch_html_report

# 从 CSV 文件生成 HTML 报告
output_path = generate_batch_html_report(
    full_csv_path="batch_results_full.csv",
    output_path="my_report.html",
    report_title="CASSIA 分析报告"
)
```

#### 从数据直接生成

```python
from CASSIA.reports.generate_batch_report import generate_batch_html_report_from_data

# 直接从数据生成 HTML 报告
output_path = generate_batch_html_report_from_data(
    rows=results_list,  # 字典列表
    output_path="my_report.html",
    report_title="CASSIA 分析报告"
)
```

这些函数创建包含以下功能的独立 HTML 报告：
- 交互式搜索功能
- 细胞类型筛选下拉菜单
- 详细信息的模态弹窗
- 嵌入式 CSS 和 JavaScript

