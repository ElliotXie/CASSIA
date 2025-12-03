---
title: Introduction to Optional Agents
---


The basic CASSIA workflow is sufficient for most situations. However, to handle some special cases, we now introduce several advanced agents.

  

## 🤔 Uncertainty Quantification Agent (UQ Agent)

Large Language Models (LLMs), while powerful, can occasionally produce variable outputs. The Uncertainty Quantification (UQ) Agent mitigates this by employing a consensus-based strategy. For example, if a model has an 80% probability of correctness on a single attempt, independently querying it five times and aggregating the results via majority voting can theoretically increase accuracy to 94.2%. The UQ Agent automates this iterative process to enhance the reliability and robustness of CASSIA's cell type annotations.

  

## 🚀 Annotation Boost Agent

This agent has appeared in the previous workflow, but you can also choose to apply it to clusters of your choice. The Annotation Boost Agent can read CASSIA's annotation reports, then continuously generate hypotheses and retrieve gene expression information from clusters to verify these hypotheses, ultimately optimizing the annotation results. In our tests, this agent performed exceptionally well and is one of CASSIA's core innovations.

  

## ⚖️ SymphonyCompare Agent

We know that the basic CASSIA workflow outputs three possible final annotations, ranked from most to least likely. Sometimes, we encounter ambiguous annotations. In such cases, the SymphonyCompare Agent can help identify the most probable annotation. This agent is actually a combination of several agents, each scoring the possible annotation results. The final annotation is chosen based on these scores.

  

## 🔬 Subclustering Agent

A single round of clustering is often not enough. Sometimes, we are particularly interested in certain clusters that need to be extracted for further analysis. The Subclustering Agent can compare multiple clusters simultaneously, providing more accurate and detailed annotations.

  

## 📚 RAG Agent (Retrieval-Augmented Generation Agent)

This agent is one of CASSIA's core innovations. Based on the basic information provided by the user, it automatically retrieves marker databases, cell ontology trees, and uses a PCA-based approach on cell types to generate a series of background information to assist the basic CASSIA workflow. When very detailed or novel annotations are needed, this agent significantly improves CASSIA's performance. Due to library dependencies, the RAG Agent is currently only available in Python. We will release an R version as soon as possible.
