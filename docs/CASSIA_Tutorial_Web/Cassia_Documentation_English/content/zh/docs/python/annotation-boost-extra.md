---
title: 注释增强 Plus（可选）
---

这可以用来研究与聚类相关的特定问题，例如推断聚类的状态。这里我们以 cd8 阳性 alpha-beta T 细胞为例。请注意，该智能体的性能尚未经过基准测试，因此请谨慎对待结果。

### 用法

```python
# 目前仅支持 openrouter 作为提供商。

CASSIA.runCASSIA_annottaionboost_additional_task(
    full_result_path = output_name + "_full.csv",
    marker = unprocessed_markers,
    output_name = "T_cell_state",
    cluster_name = "cd8-positive, alpha-beta t cell",  # 高线粒体含量的聚类
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-3.5-sonnet",
    additional_task = "infer the state of this T cell cluster"
)
```

