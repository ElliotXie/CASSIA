---
title: 简介
---

欢迎来到 CASSIA Python 教程！本指南将带您了解如何在 Python 中使用 CASSIA 包进行细胞类型注释。


## 什么是 CASSIA？
CASSIA (Collaborative Agent System for Single-cell Interpretable Annotation) 旨在利用多智能体大语言模型 (LLM) 的能力来增强细胞类型注释。我们将在本教程中介绍每个智能体。

![CASSIA Agent Workflow](/images/agent.webp)

*图 1：CASSIA 的多智能体工作流程。(a) 用户提供输入，包括物种、组织、标记基因和模型选择。(b) 核心智能体通过注释、验证、格式化、评分和报告处理数据。(c) CASSIA 输出全面的注释结果，包括细胞类型、亚型、混合群体和质量评分。(d) 可选的专用智能体为复杂分析提供额外功能。*

## CASSIA 的表现如何？

下面的热图显示，CASSIA 在各种物种和组织中都优于大多数其他无参考方法。

![CASSIA Performance Comparison](/images/cassia-comparison.webp)

*图 2：CASSIA 与其他细胞类型注释方法在不同数据集上的性能比较。数值越高（红色）表示性能越好。*

让我们开始为您的单细胞分析工作流程设置 CASSIA！

