# CASSIA model card — cross-model comparison (FULL)
Judge (fixed referee): composer-2.5 · single run per model (LLM nondeterminism applies).

| metric | composer-2.5 | opus4.8-xhigh | gpt5.5 |
|---|---|---|---|
| **Overall correct** | **62%** | **70%** | **62%** |
| Overall +partial | 72% | 82% | 82% |
| annot. cost (API-equiv) | $0.98 | $6.30 | $4.67 |
| |  |  |  |
| _by solve tier_ ||||
| basic | 88% | 96% | 88% |
| boost | 42% | 47% | 32% |
| frontier | 17% | 33% | 50% |
| |  |  |  |
| _by failure mode_ ||||
| canonical | 100% | 100% | 100% |
| subtype_resolution | 67% | 75% | 75% |
| near_neighbor_state | 69% | 77% | 54% |
| buried_marker | 50% | 67% | 50% |
| spatial_hard | 17% | 17% | 17% |
| rna_protein_discordant | 0% | 33% | 33% |
| |  |  |  |
| _by species_ ||||
| human | 67% | 69% | 67% |
| mouse | 50% | 71% | 50% |
