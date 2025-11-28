---
title: 注释增强扩展（可选）
---

该智能体是注释增强智能体的扩展版本。它使用自定义分析任务对特定簇执行额外分析。此功能通过允许使用大型语言模型对簇特定的生物学问题进行集中调查，扩展了基本注释功能。

## 用法
```r
runCASSIA_annottaionboost_additional_task(
    full_result_path,
    marker,
    output_name,
    cluster_name,
    major_cluster_info,
    num_iterations = 5,
    model = "anthropic/claude-3.5-sonnet",
    additional_task
)
```

## 参数
* `full_result_path`：字符串。先前分析的完整结果CSV文件的路径。应包含带有"_full.csv"后缀的输出名称。

* `marker`：包含分析中使用的未处理标记数据的对象。

* `output_name`：字符串。此函数生成的输出文件的名称。

* `cluster_name`：字符串。要分析的特定簇的名称（例如，"cd8-positive, alpha-beta t cell"）。

* `major_cluster_info`：字符串。关于正在研究的组织或系统的一般信息（例如，"Human Large Intestine"人类大肠）。

* `num_iterations`：整数。运行分析的迭代次数。默认为5。

* `model`：字符串。指定要使用的大型语言模型。目前仅支持"anthropic/claude-3.5-sonnet"。

* `additional_task`：字符串。对簇执行的自定义分析任务（例如，"infer the state of this T cell cluster"推断这个T细胞簇的状态）。

## 详情
此函数通过使用大型语言模型执行专业分析任务来增强簇注释。它特别适用于调查有关单个簇的特定生物学问题。该函数将标记数据与自定义分析目标集成，以提供对簇特征的更深入见解。

## 注意
这种分析方法的性能尚未进行广泛的基准测试。结果应谨慎解释，并通过其他方法进行验证。

## 示例
```r
runCASSIA_annottaionboost_additional_task(
    full_result_path = "output_name_full.csv",
    marker = markers_unprocessed,
    output_name = "T_cell_state",
    cluster_name = "cd8-positive, alpha-beta t cell",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-3.5-sonnet",
    additional_task = "infer the state of this T cell cluster"
)
```
