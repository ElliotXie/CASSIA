---
title: Uncertainty Quantification (Optional)
---


Uncertainty quantification in CASSIA helps assess annotation reliability through multiple analysis iterations and similarity scoring. This process is crucial for:
- Identifying robust cell type assignments
- Detecting mixed or ambiguous clusters
- Quantifying annotation confidence
- Understanding prediction variability

### Multiple Iteration Analysis

#### Basic Usage
```R

# Run multiple analyses
runCASSIA_batch_n_times(
    # Core parameters
    n = 5, #number of iteratioins
    marker = marker_data,
    output_name = "my_annottaion_repeat",
    
    # Model settings
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter",
    
    # Context information
    tissue = "brain",
    species = "human",
    additional_info = NULL,


    # Processing control
    max_workers = 4,        # Total parallel workers
    batch_max_workers = 2   # Workers per batch
)
```

> **⚠️ Cost Warning**: Running multiple iterations with LLM models can incur significant costs. Each iteration makes separate API calls, so the total cost will be approximately n times the cost of a single run. Consider starting with a smaller number of iterations for testing purposes.

#### Parameter Details

1. **Iteration Control**:
   - `n`: Number of analysis iterations
   - Recommended: 5 iterations for standard analysis
   - Consider more iterations for critical applications

2. **Resource Management**:
   - `max_workers`: Overall parallel processing limit
   - `batch_max_workers`: Workers per iteration
   - max_workers * batch_max_workers to match your number of cores.


### Similarity Score Calculation

#### Running Similarity Analysis
```R

# Calculate similarity scores
runCASSIA_similarity_score_batch(
    # Input parameters
    marker = marker_data,
    file_pattern = "my_annottaion_repeat_*_full.csv",
    output_name = "similarity_results",
    
    
    # Processing parameters
    max_workers = 4,
    model = "google/gemini-2.5-flash-preview",
    provider = "openrouter",
    
    # Scoring weights
    main_weight = 0.5, # Weight for main cell type
    sub_weight = 0.5  # Weight for subtypes
)
```

#### Scoring Parameters

1. **Weight Configuration**:
   - `main_weight`: Importance of main cell type match (0-1)
   - `sub_weight`: Importance of subtype match (0-1)
   - Weights should sum to 1.0

2. **File Management**:
   - `file_pattern`: Pattern to match iteration results
   - Uses * to match iteration numbers
   - Example:  if you have "my_annottaion_repeat_1_full.csv", "my_annottaion_repeat_2_full.csv", and "my_annottaion_repeat_3_full.csv", use "my_annottaion_repeat__full.csv" to match the pattern.

#### Output Interpretation

1. **Similarity Scores**:
   - Range: 0 (completely different) to 1 (identical)
   - Interpretation guidelines:
     - 0.9: High consistency
     - 0.75-0.9: Moderate consistency
     - <0.75: Low consistency

#### Troubleshooting

1. **Performance Issues**:
   - Reduce worker counts
   - Process in smaller batches

2. **Low Similarity Scores**:
   - Review marker gene quality
   - Use Annotation Boost function
   - Review cluster heterogeneity
   - Consider biological variability
   - Increase iteration count
   - Try subclustering
