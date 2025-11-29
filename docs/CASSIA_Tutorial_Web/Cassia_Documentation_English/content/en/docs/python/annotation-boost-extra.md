---
title: Annotation Boost Plus (Optional)
---

This can be used to study a given problem related to a cluster, such as infer the state of a cluster. Here we use the cd8-positive, alpha-beta t cell as an example. Note that the performance of this agent has not been benchmarked, so please be cautious with the results.

### Usage

```python
#only openrouter is supported as provider now.

CASSIA.runCASSIA_annottaionboost_additional_task(
    full_result_path = output_name + "_full.csv",
    marker = unprocessed_markers,
    output_name = "T_cell_state",
    cluster_name = "cd8-positive, alpha-beta t cell",  # Cluster with high mitochondrial content
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-3.5-sonnet",
    additional_task = "infer the state of this T cell cluster"
)
```

