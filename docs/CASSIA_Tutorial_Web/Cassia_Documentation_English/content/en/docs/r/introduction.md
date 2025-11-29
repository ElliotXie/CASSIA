---
title: Introduction
---

Welcome to the CASSIA R tutorial! This guide will walk you through using the CASSIA package in R for cell type annotation. 


## What is CASSIA?
CASSIA (Collaborative Agent System for Single-cell Interpretable Annotation) is designed to enhance cell type annotation by leveraging the power of multi-agent large language models (LLMs). We will introduce each of the agent in this tutorial.

![CASSIA Agent Workflow](/images/agent.png)

*Figure 1: CASSIA's multi-agent workflow. (a) User provides input including species, tissue, marker genes, and model selection. (b) Core agents process the data through annotation, validation, formatting, scoring, and reporting. (c) CASSIA outputs comprehensive annotation results including cell types, subtypes, mixed populations, and quality scores. (d) Optional specialized agents provide additional capabilities for complex analyses.*

## How does CASSIA perform?

The heatmap below shows that CASSIA has outperformed most of the other reference-free methods in various species and tissues.

![CASSIA Performance Comparison](/images/cassia-comparison.jpg)

*Figure 2: Performance comparison between CASSIA and other cell type annotation methods across different datasets. Higher values (red) indicate better performance.*

Let's get started with setting up CASSIA for your single-cell analysis workflow!
