---
title: 扩展分析
---

本教程涵盖 CASSIA 的高级功能，包括不确定性量化、注释增强和专用智能体。

## 可选步骤：不确定性量化

这对于研究注释的不确定性并可能提高准确性很有用。

注意：此步骤可能会很昂贵，因为将执行多次迭代。

```python
output_name="intestine_detailed"

# 运行多次迭代
iteration_results = CASSIA.runCASSIA_batch_n_times(
    n=2,
    marker=unprocessed_markers,
    output_name=output_name + "_Uncertainty",
    model="anthropic/claude-sonnet-4.5",
    provider="openrouter",
    tissue="large intestine",
    species="human",
    max_workers=6,
    batch_max_workers=3  # API 速率限制的保守设置
)


# 计算相似性评分
similarity_scores = CASSIA.runCASSIA_similarity_score_batch(
    marker=unprocessed_markers,
    file_pattern=output_name + "_Uncertainty_*_full.csv",
    output_name="intestine_uncertainty",
    max_workers=6,
    model="openai/gpt-5.1",
    provider="openrouter",
    main_weight=0.5,
    sub_weight=0.5
)
```

## 可选步骤：对选定聚类进行注释增强

单核细胞聚类有时被注释为免疫细胞和神经元/胶质细胞的混合群体。

这里我们使用注释增强智能体来更详细地测试这些假设。

```python
# 对高线粒体含量聚类运行增强验证
CASSIA.runCASSIA_annotationboost(
    full_result_path = output_name + "_full.csv",
    marker = unprocessed_markers,
    output_name = "monocyte_annotationboost",
    cluster_name = "monocyte",
    major_cluster_info = "Human Large Intestine",
    num_iterations = 5,
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter"
)
```

## 可选步骤：检索增强生成 (RAG)

这在您需要处理非常具体和详细的注释时特别有用。它可以显着提高注释的粒度和准确性。它自动提取标记信息并生成报告作为默认 CASSIA 流程的附加信息。

```bash
pip install cassia-rag
```

```python
from cassia_rag import run_complete_analysis
import os

os.environ["ANTHROPIC_API_KEY"] = "your-anthropic-key"
os.environ["OPENAI_API_KEY"] = "your-openai-key"

run_complete_analysis(
        tissue_type="Liver", # 您正在分析的组织
        target_species="Tiger", # 您正在分析的物种
        reference_species="Human", # 人类或小鼠，如果是其他物种，则使用人类代替小鼠
        model_choice='claude', # claude 或 gpt，强烈推荐 claude
        compare=True,  # 如果您想与参考物种进行比较，例如胎儿与成人，则设置为 True
        db_path="~/Canonical_Marker (1).csv", # 数据库路径
        max_workers=8
)
```

## 可选步骤：使用多个 LLM 比较亚型 (Symphony Compare)

此智能体在运行默认 CASSIA 流程后特别有用，如果您不确定特定聚类的身份。您可以使用此智能体获得更自信的亚型注释。这里我们以浆细胞聚类为例。区分它更像是普通浆细胞还是其他细胞类型。

```python
# 这里的标记复制自 CASSIA 之前的结果。
marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

# 运行 Symphony Compare
results = CASSIA.symphonyCompare(
    tissue = "large intestine",
    celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
    marker_set = marker,
    species = "human",
    model_preset = "premium",
    output_basename = "plasma_cell_subtype",
    enable_discussion = True
)

print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
```

## 可选步骤：亚群聚类

此智能体可以用来研究亚群聚类群体，例如 T 细胞群体或成纤维细胞聚类。我们建议首先应用默认的 CASSIA，然后在目标聚类上，应用 Seurat 流程进行亚群聚类并获取标记结果以在此处使用。这里我们展示了 cd8 阳性 alpha-beta T 细胞聚类的结果作为示例。此聚类是与其他细胞类型混合的 cd8 群体。

```python
CASSIA.runCASSIA_subclusters(marker = subcluster_results,
    major_cluster_info = "cd8 t cell",
    output_name = "subclustering_results",
    model = "anthropic/claude-sonnet-4.5",
    provider = "openrouter")
```

建议运行 CS 评分以获得更可信的答案。

```python
CASSIA.runCASSIA_n_subcluster(
    n=5, 
    marker=subcluster_results,
    major_cluster_info="cd8 t cell", 
    base_output_name="subclustering_results_n",
    model="anthropic/claude-sonnet-4.5",
    temperature=0,
    provider="openrouter",
    max_workers=5,
    n_genes=50
)

# 计算相似性评分
CASSIA.runCASSIA_similarity_score_batch(
    marker = subcluster_results,
    file_pattern = "subclustering_results_n_*.csv",
    output_name = "subclustering_uncertainty",
    max_workers = 6,
    model = "openai/gpt-5.1",
    provider = "openrouter",
    main_weight = 0.5,
    sub_weight = 0.5
)
```

