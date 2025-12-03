---
title: RAG（可选）
---

这在您需要处理非常具体和详细的注释时特别有用。它可以显着提高注释的粒度和准确性。它自动提取标记信息并生成报告作为默认 CASSIA 流程的附加信息。

### 安装

```bash
pip install cassia-rag
```

### 用法

```python
from cassia_rag import run_complete_analysis
import os

# 如果尚未设置，请设置 API 密钥。
os.environ["ANTHROPIC_API_KEY"] = "your-anthropic-key"
os.environ["OPENAI_API_KEY"] = "your-openai-key"

# 运行包装函数以触发多智能体管道。
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

### 输出文件

所有输出（中间和最终）都保存在 "TissueType_Species" 文件夹中的 txt 文件中。在我们的示例中，它是 Liver_Tiger 文件夹。
最终输出在 `summary_clean.txt` 文件中。该文件的内容稍后可用作 CASSIA 流程中的附加信息。

文件夹中还有一些其他文件，这些是中间输出。
以教程输入为例，文件为：
1. `liver_tiger_marker_analysis.txt` # 来自数据库的标记分析和解释
2. `final_ontology.txt` # 与组织类型和目标物种相关的本体
3. `cell_type_patterns_claude.txt` # 细胞类型模式分析
4. `summary.txt` # 原始摘要文件
5. `additional_considerations.txt` # 如果我们有与参考物种不同的物种，则需要考虑的其他事项。

