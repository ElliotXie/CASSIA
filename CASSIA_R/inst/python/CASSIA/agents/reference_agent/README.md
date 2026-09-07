# Reference Agent — 专家参考文档检索 Agent

CASSIA 的可选 agent 之一。在 `use_reference=True` 时，从一个专家手写/文献整理的参考知识库里挑出最相关的几份文档，把内容注入到 annotation prompt 里给主 agent 参考。

**设计哲学：开不开启由用户决定。** 用户清楚自己在做什么（比如正在对一群 macrophage 做 subtype 细分），所以不需要 LLM 再替他判断"要不要用 reference"。一旦启用，就直接进入选文件 + 抽内容的流程，行为可预期、每 cluster 少一次 API 调用。

---

## 当前设计

### Subclustering 的 agentic 工作流

Subclustering 场景现在优先走一个更 agentic 的流程。Reference agent
不是逐个 cluster 做 marker matching，而是先读取 lineage overview（例如
`references_brain/myeloid/macrophage/_overview.md`，相当于这个 lineage 的
`SKILL.md`/router），再同时查看本次所有 subcluster marker sets，形成全局
subtype landscape 假设，然后决定还需要读哪些更细的 reference 文档。

```
   All subcluster marker sets + parent cluster context
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  Read lineage overview/router                       │
│  e.g. myeloid/macrophage/_overview.md               │
└────────────┬────────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────┐
│  Plan reference reads                               │
│  LLM decides which detailed docs are needed          │
│  e.g. tam_pan_cancer.md, resident_like.md            │
└────────────┬────────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────┐
│  Read selected docs only                            │
│  Avoid stuffing the whole reference library          │
└────────────┬────────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────┐
│  Synthesize case-specific reference brief            │
│  Includes objective paper facts, cluster guidance,   │
│  and cross-cluster conflict resolution               │
└────────────┬────────────────────────────────────────┘
             │
             ▼
  <expert_reference><reference_brief>...</reference_brief></expert_reference>
```

这个流程更接近 annotation boost 的设计：LLM 使用工具式步骤读 overview、
选文档、读细节、再写 brief。最后注入主 subclustering prompt 的不是原始文档
拼接，而是一段针对当前数据集的专家参考简报。

### Single-cluster 兼容工作流

单 cluster annotation 仍保留原来的单步 reference 选择接口：

```
   Markers (top 20) + tissue/species + 可选 cell_type_hint
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  Reference Selection (complexity_scorer.py)         │
│                                                     │
│  LLM 一次调用，看 markers + context + _router.md     │
│  同时输出:                                           │
│   - preliminary_cell_type                           │
│   - cell_type_range (候选列表)                      │
│   - selected_references (0-3 个文件路径)            │
│   - reasoning                                       │
│                                                     │
│  若 selected_references 为空 → 视为"库里没合适的"     │
│    should_use_reference=False，不注入              │
└────────────┬────────────────────────────────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────┐
│  Section Extraction (section_extractor.py)          │
│                                                     │
│  解析选中的 markdown:                                │
│   - 按 cell type guess + marker 命中度筛选段落       │
│   - 按 depth ('detailed' / 'summary') 决定粒度       │
│   - 截断到 max_content_length (默认 8000)           │
└────────────┬────────────────────────────────────────┘
             │
             ▼
  <expert_reference> block 注入到主 annotation prompt
```

### 关键对外接口

```python
from CASSIA.agents.reference_agent import ReferenceAgent, format_reference_for_prompt

agent = ReferenceAgent(provider="openrouter", model="google/gemini-3.8-flash")

result = agent.get_reference_for_markers(
    markers=["SPP1", "MMP9", "VEGFA", "C1QA", ...],
    tissue="tumor",
    species="human",
    cell_type_hint="macrophage",   # 可选。subcluster 场景下直接告诉 agent 父级谱系
    depth="detailed",              # "detailed" 或 "summary"
    max_content_length=8000,
)

# result 字段:
# - should_use_reference: bool
# - content: str  (format 好的 markdown 片段)
# - references_used: List[str]
# - preliminary_cell_type, cell_type_range, selected_references, reasoning
```

与主管线集成：
```python
runCASSIA(
    marker_list=...,
    use_reference=True,
    reference_cell_type_hint="macrophage",   # 可选
    reference_provider=..., reference_model=...,
)
```
单 cluster (`runCASSIA`) 和 batch (`runCASSIA_batch`) 都支持。

### 为什么不再做 complexity gating

之前版本先让 LLM 判断"当前 markers 够不够复杂、需不需要 reference"，再决定要不要选文件。实践中这一步没太大价值：
- 用户开 `use_reference` 时通常就是知道场景需要 reference（例如全是 macrophage 的 subcluster）
- 多一次 LLM 调用，batch 跑几十个 cluster 时成本翻倍
- LLM 的 complexity 判断本身就不稳定，用户也难以预期哪些 cluster 会真的注入

现在统一到一次调用：LLM 看到 router，要么选出相关文件，要么返回空列表。"空列表"天然替代了"不需要 reference"。

---

## 目录结构

```
reference_agent/
├── reference_agent.py        # 主 orchestrator（ReferenceAgent 类）
├── complexity_scorer.py      # 单次 LLM 调用：推断 cell type + 选 reference
├── reference_selector.py     # 独立工具：按类别/marker 命中查找（list_available_references 用）
├── section_extractor.py      # Markdown section 解析 / marker 匹配
├── utils.py                  # 索引加载、路径、内容格式化
├── __init__.py
│
├── prompts/
│   └── extract_markers_from_paper.md   # 用于从文献抽取 marker 的提示词
│
├── references_brain/         # 知识库（专家手写 + 文献整理）
│   ├── _router.md            # 目录/路由表，LLM 看这个决定选哪个文件
│   ├── b_cell/_overview.md
│   ├── myeloid/_overview.md
│   └── t_cell/
│       ├── _overview.md
│       ├── cd4/_overview.md
│       └── cd8/_overview.md
│
└── macrophage_test/          # 实验子目录：为 macrophage 扩充知识库
    ├── papers/
    │   ├── candidate_papers.md   # 候选文献清单（注释版）
    │   ├── methods.md            # batch downloader 的输入格式
    │   └── downloads/            # 已下载的 PDF + supplements
    ├── reference_drafts/     # 手写 reference 草稿待填
    ├── test_data/            # 测试数据集待填
    └── tools/                # 复制过来的下载/后处理脚本
        ├── download_paper.py
        ├── batch_download_from_md.py
        └── postprocess_markdown.py
```

### Reference 文件格式

每个 `.md` 文件头部带 YAML frontmatter，正文是自由 markdown：

```markdown
---
id: myeloid_overview
category: myeloid
cell_types:
  - Macrophage
  - Monocyte
trigger_markers:
  - CD14
  - CD68
  - CSF1R
exclusion_markers:
  - CD3D
  - CD19
---

# Myeloid Cell Annotation Guide
## Overview
...
## Key Canonical Markers
...
## Subtype Differentiation
...
```

- `trigger_markers`：命中任一即可被 router 候选
- `exclusion_markers`：命中则排除（避免 T/B 细胞误选到髓系）
- 正文用 `## Subsection` 分块，`section_extractor.py` 会按 cell type guess / marker 命中度筛选相关段落，而不是全文塞进 prompt

---

## 当前进度

### 已完成
- 两步 ReAct 工作流跑通
- 与主 annotation 管线 (`runCASSIA`) 集成
- 测试：[Test/12_batch_with_reference/](../../../../Test/12_batch_with_reference/)
- 目录/路由/解析代码全部就位
- Macrophage pilot 的第一版 subtype 文档已加入：
  - `references_brain/myeloid/macrophage/_overview.md`
  - `references_brain/myeloid/macrophage/tam_pan_cancer.md`
  - `references_brain/myeloid/macrophage/resident_like.md`
  - `references_brain/myeloid/macrophage/inflammatory_interferon.md`
- Subclustering 入口支持 `use_reference=True`，可用
  `reference_cell_type_hint="macrophage"` 将 subtype reference 注入 prompt。

### 尚不完整（瓶颈在这里）
- `references_brain/` 仍然很稀疏。Macrophage 有第一版 subtype docs，
  但 fibroblast / broader immune-cell subtype references 还没有系统扩充。
- 需要继续扩充到更多亚型文档（如更细的 `tam_spp1.md`,
  `tam_trem2.md`, `fibroblast/inflammatory.md`, `fibroblast/myofibroblast.md`），
  才能覆盖更多真实 subclustering 场景。

### Macrophage 试点子目录 (`macrophage_test/`)
作为第一个系统化扩充的 pilot：
- 候选文献 15 篇，已成功下载 10 篇 PDF + 40 个 supplementary 文件
  - 核心命中：Cheng 2021 Cell pan-cancer 髓系 atlas、Mulder 2021 Immunity MoMac-VERSE、Ma 2024 Nat Commun TAM×ICB、Wang 2023 Cell 胎儿巨噬细胞 atlas 等
  - 缺的 3 篇在 Science Immunology / Blood 付费墙，无 bioRxiv 预印本
- 下一步要做的是**从 supplementary excel 里提取 marker 表** → 写成 `references_brain/myeloid/macrophage/tam_*.md` 一套文档

---

## 扩充知识库的推荐流程

1. 在 `*_test/papers/methods.md` 里列出候选文献（格式见该文件示例，`- **Field**: value`，冒号必须在 `**` 外）
2. 运行 `tools/batch_download_from_md.py methods.md -o downloads/` 自动下载 PDF + supplements
3. 从 supplementary excel / 正文 marker table 提取 subtype → marker 映射
4. 写成单个 `.md` 文件放进 `references_brain/<category>/<subfolder>/<subtype>.md`，注意带齐 YAML frontmatter
5. 把该文件追加到 `_router.md` 目录表
6. 跑一次 `Test/12_batch_with_reference/` 验证 agent 能正确选中新文档
7. 对 subclustering 跑 `Test/27_reference_subtype/`，再用
   `benchmark_macrophage_reference.py` 做 live Kimi benchmark（需要 `OPENROUTER_API_KEY`）。

---

## 已知待改进

- `references_brain/` 内容稀疏——是目前最大的瓶颈
- Router 当前是手写纯 markdown 目录，没做自动生成；加了新文档要手动更新 `_router.md`
- Section 抽取是 keyword-based，没有向量检索——对于正文很长的参考文档可能召回不全
- Complexity threshold (默认 40) 是拍脑袋定的，缺乏系统 benchmark
