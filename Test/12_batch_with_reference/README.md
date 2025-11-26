# Test 12: Batch Annotation with Reference Retrieval

## Purpose
Tests the `runCASSIA_batch_with_reference()` function with `use_reference=True`,
which enables the two-step ReAct reference retrieval workflow.

## What it Tests
- Per-cluster reference retrieval using LLM-based complexity assessment
- Two-step workflow: assess complexity -> select references from router
- Reference content injection into annotation prompts

## Expected Output
- All 6 clusters annotated successfully
- "Reference Used" column in output CSV
- Reference retrieval statistics printed during execution

## Run
```bash
python test_batch_with_reference.py
```
