#!/usr/bin/env python
"""
Test script for the merge_annotations function from merging_annotation.py
"""

import os
import sys
import time
import pandas as pd
from dotenv import load_dotenv
from merging_annotation import merge_annotations, merge_annotations_all

# Set pandas display options
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.width', 1000)     # Wider display
pd.set_option('display.max_columns', None)  # Show all columns

# Load environment variables from .env file if it exists
load_dotenv()

# Create a sample CSV with test data
def create_test_csv(filepath):
    """Create a sample CSV file for testing with the new column naming scheme"""
    data = {
        "True Cell Type": [1, 2, 3, 4, 5, 6, 7, 8],  # This is actually the cluster ID column
        "Predicted Main Cell Type": [
            "macrophage", 
            "CD4 T cell", 
            "B cell", 
            "dendritic cell", 
            "CD8 T cell", 
            "NK cell", 
            "epithelial cell", 
            "fibroblast"
        ],
        "Predicted Sub Cell Types": [
            "inflammatory macrophage, resident macrophage", 
            "naive CD4 T cell, memory CD4 T cell", 
            "memory B cell, plasma cell", 
            "plasmacytoid dendritic cell, conventional dendritic cell", 
            "cytotoxic CD8 T cell, exhausted CD8 T cell", 
            "CD56bright NK cell, CD56dim NK cell", 
            "type II pneumocyte, type I pneumocyte", 
            "activated fibroblast, quiescent fibroblast"
        ]
    }
    
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Created test CSV file at: {filepath}")
    return df

def test_sequential_vs_parallel():
    """Test sequential vs parallel processing of all three detail levels"""
    print("\n=== Testing Sequential vs Parallel Processing ===")
    
    # Create test data
    test_csv = "test_clusters.csv"
    create_test_csv(test_csv)
    
    # Additional context to help with annotations
    additional_context = """
    Cell type reference information:
    - Macrophages and dendritic cells are types of myeloid cells
    - CD4 and CD8 T cells belong to T lymphocyte lineage
    - B cells are part of the B lymphocyte lineage
    - NK cells are natural killer cells and considered part of innate lymphoid cells
    - Epithelial cells form the tissue lining organs
    - Fibroblasts are connective tissue cells
    """
    
    # Choose provider based on available API keys
    provider = "openrouter"
    model = "deepseek/deepseek-chat-v3-0324"
    
    # Check if OpenRouter API key is available
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("OpenRouter API key not found. Please set OPENROUTER_API_KEY environment variable.")
        return False
    
    try:
        # Measure time for sequential processing
        print("\nRunning sequential processing...")
        start_time_sequential = time.time()
        
        # Process each level sequentially
        merge_annotations(csv_path=test_csv, output_path=None, provider=provider, model=model, 
                          additional_context=additional_context, detail_level="broad")
        merge_annotations(csv_path=test_csv, output_path=None, provider=provider, model=model, 
                          additional_context=additional_context, detail_level="detailed")
        merge_annotations(csv_path=test_csv, output_path=None, provider=provider, model=model, 
                          additional_context=additional_context, detail_level="very_detailed")
        
        end_time_sequential = time.time()
        sequential_time = end_time_sequential - start_time_sequential
        print(f"Sequential processing completed in {sequential_time:.2f} seconds")
        
        # Measure time for parallel processing
        print("\nRunning parallel processing...")
        start_time_parallel = time.time()
        
        # Process all levels in parallel
        result_df = merge_annotations_all(csv_path=test_csv, output_path="all_groupings.csv", 
                                         provider=provider, model=model, 
                                         additional_context=additional_context)
        
        end_time_parallel = time.time()
        parallel_time = end_time_parallel - start_time_parallel
        print(f"Parallel processing completed in {parallel_time:.2f} seconds")
        
        # Calculate speedup
        speedup = sequential_time / parallel_time if parallel_time > 0 else float('inf')
        print(f"Parallel speedup: {speedup:.2f}x")
        
        # Display results
        print("\nResults from parallel processing:")
        print(result_df[["True Cell Type", "Predicted Main Cell Type", 
                         "Merged_Grouping_1", "Merged_Grouping_2", "Merged_Grouping_3"]].head())
        
        return True
        
    except Exception as e:
        print(f"Error in performance test: {str(e)}")
        return False
    finally:
        # Clean up test file
        if os.path.exists(test_csv):
            try:
                os.remove(test_csv)
                print(f"Removed test file: {test_csv}")
            except:
                pass
        if os.path.exists("all_groupings.csv"):
            try:
                os.remove("all_groupings.csv")
                print(f"Removed test file: all_groupings.csv")
            except:
                pass

def test_parallel_all_levels():
    """Test parallel processing of all three detail levels."""
    print("\n=== Testing Parallel Processing of All Detail Levels ===")
    
    # Create test data
    test_csv = "test_clusters.csv"
    create_test_csv(test_csv)
    
    # Additional context to help with annotations
    additional_context = """
    Cell type reference information:
    - Macrophages and dendritic cells are types of myeloid cells
    - CD4 and CD8 T cells belong to T lymphocyte lineage
    - B cells are part of the B lymphocyte lineage
    - NK cells are natural killer cells and considered part of innate lymphoid cells
    - Epithelial cells form the tissue lining organs
    - Fibroblasts are connective tissue cells
    """
    
    # Choose provider based on available API keys
    provider = "openrouter"
    model = "deepseek/deepseek-chat-v3-0324"
    
    # Check if OpenRouter API key is available
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("OpenRouter API key not found. Please set OPENROUTER_API_KEY environment variable.")
        return False
    
    try:
        print("Processing all three detail levels in parallel...")
        
        # Process all levels in parallel and save to CSV
        all_results_df = merge_annotations_all(
            csv_path=test_csv,
            output_path="all_groupings.csv",
            provider=provider,
            model=model,
            additional_context=additional_context,
            batch_size=8
        )
        
        # Display the results
        print("\nResults from parallel processing (all three detail levels):")
        
        # Create a formatted comparison table
        comparison_df = pd.DataFrame({
            'Cluster': all_results_df['True Cell Type'],
            'Cell Type': all_results_df['Predicted Main Cell Type'],
            'Subtype': all_results_df['Predicted Sub Cell Types'].apply(
                lambda x: x.split(',')[0].strip() if isinstance(x, str) else x
            ),
            'Broad': all_results_df['Merged_Grouping_1'],
            'Detailed': all_results_df['Merged_Grouping_2'],
            'Very Detailed': all_results_df['Merged_Grouping_3']
        })
        
        # Save comparison to CSV
        comparison_df.to_csv("all_levels_comparison.csv", index=False)
        print("Saved comparison to all_levels_comparison.csv")
        
        # Print comparison
        print("\nCell groupings at all three detail levels:")
        for _, row in comparison_df.iterrows():
            print(f"Cluster {row['Cluster']} ({row['Cell Type']} / {row['Subtype']}):")
            print(f"  - Broad:        {row['Broad']}")
            print(f"  - Detailed:     {row['Detailed']}")
            print(f"  - Very Detailed: {row['Very Detailed']}")
            print()
        
        return True
        
    except Exception as e:
        print(f"Error processing all levels in parallel: {str(e)}")
        return False

def main():
    """Run test for merge_annotations"""
    print("Testing merge_annotations function...")
    
    # Check if any API keys are available
    has_api_key = False
    for env_var in ["OPENAI_API_KEY", "ANTHROPIC_API_KEY", "OPENROUTER_API_KEY"]:
        if os.environ.get(env_var):
            has_api_key = True
            break
    
    if not has_api_key:
        print("No API keys found. Please set at least one of: OPENAI_API_KEY, ANTHROPIC_API_KEY, OPENROUTER_API_KEY")
        return 1
    
    # Run sequential vs parallel test if --perf flag is provided
    if len(sys.argv) > 1 and sys.argv[1] == "--perf":
        test_sequential_vs_parallel()
        return 0
    
    # Run all levels in parallel if --all flag is provided
    if len(sys.argv) > 1 and sys.argv[1] == "--all":
        test_parallel_all_levels()
        return 0
        
    # Create test data
    test_csv = "test_clusters.csv"
    create_test_csv(test_csv)
    
    # Additional context to help with annotations
    additional_context = """
    Cell type reference information:
    - Macrophages and dendritic cells are types of myeloid cells
    - CD4 and CD8 T cells belong to T lymphocyte lineage
    - B cells are part of the B lymphocyte lineage
    - NK cells are natural killer cells and considered part of innate lymphoid cells
    - Epithelial cells form the tissue lining organs
    - Fibroblasts are connective tissue cells
    """
    
    # Choose provider based on available API keys
    provider = "openrouter"
    model = "deepseek/deepseek-chat-v3-0324"
    
    # Check if OpenRouter API key is available
    if not os.environ.get("OPENROUTER_API_KEY"):
        print("OpenRouter API key not found. Please set OPENROUTER_API_KEY environment variable.")
        return 1
    
    # Process with broad detail level (default)
    print("\n=== Testing with BROAD detail level ===")
    try:
        broad_df = merge_annotations(
            csv_path=test_csv,
            output_path="broad_groupings.csv",
            provider=provider,
            model=model,
            additional_context=additional_context,
            batch_size=8,  # Process all rows in one batch for this test
            detail_level="broad"
        )
        
        # Display results
        print("\nBroad groupings:")
        broad_df.to_csv("broad_groupings_readable.csv", index=False)
        print(f"Saved broad grouping results to broad_groupings_readable.csv")
        
        # Print row by row for reliable console output
        for idx, row in broad_df.iterrows():
            # Extract first subtype for display
            subtype = row['Predicted Sub Cell Types'].split(',')[0].strip() if isinstance(row['Predicted Sub Cell Types'], str) else row['Predicted Sub Cell Types']
            print(f"Cluster {row['True Cell Type']}: {row['Predicted Main Cell Type']} / {subtype} → {row['Merged_Grouping_1']}")
    
    except Exception as e:
        print(f"Error with broad groupings: {str(e)}")
    
    # Process with detailed level
    print("\n=== Testing with DETAILED detail level ===")
    try:
        detailed_df = merge_annotations(
            csv_path=test_csv,
            output_path="detailed_groupings.csv",
            provider=provider,
            model=model,
            additional_context=additional_context,
            batch_size=8,  # Process all rows in one batch for this test
            detail_level="detailed"
        )
        
        # Display results
        print("\nDetailed groupings:")
        detailed_df.to_csv("detailed_groupings_readable.csv", index=False)
        print(f"Saved detailed grouping results to detailed_groupings_readable.csv")
        
        # Print row by row for reliable console output
        for idx, row in detailed_df.iterrows():
            # Extract first subtype for display
            subtype = row['Predicted Sub Cell Types'].split(',')[0].strip() if isinstance(row['Predicted Sub Cell Types'], str) else row['Predicted Sub Cell Types']
            print(f"Cluster {row['True Cell Type']}: {row['Predicted Main Cell Type']} / {subtype} → {row['Merged_Grouping_2']}")
    
    except Exception as e:
        print(f"Error with detailed groupings: {str(e)}")
        
    # Process with very detailed level
    print("\n=== Testing with VERY DETAILED level ===")
    try:
        very_detailed_df = merge_annotations(
            csv_path=test_csv,
            output_path="very_detailed_groupings.csv",
            provider=provider,
            model=model,
            additional_context=additional_context,
            batch_size=8,  # Process all rows in one batch for this test
            detail_level="very_detailed"
        )
        
        # Display results
        print("\nVery detailed groupings:")
        very_detailed_df.to_csv("very_detailed_groupings_readable.csv", index=False)
        print(f"Saved very detailed grouping results to very_detailed_groupings_readable.csv")
        
        # Print row by row for reliable console output
        for idx, row in very_detailed_df.iterrows():
            # Extract first subtype for display
            subtype = row['Predicted Sub Cell Types'].split(',')[0].strip() if isinstance(row['Predicted Sub Cell Types'], str) else row['Predicted Sub Cell Types']
            print(f"Cluster {row['True Cell Type']}: {row['Predicted Main Cell Type']} / {subtype} → {row['Merged_Grouping_3']}")
    
        # Print comparison of all three levels
        print("\n=== Comparison of all detail levels ===")
        comparison_df = pd.DataFrame({
            'Cluster': broad_df['True Cell Type'],
            'Cell Type': broad_df['Predicted Main Cell Type'],
            'Subtype': broad_df['Predicted Sub Cell Types'].apply(lambda x: x.split(',')[0].strip() if isinstance(x, str) else x),
            'Broad Grouping': broad_df['Merged_Grouping_1'],
            'Detailed Grouping': detailed_df['Merged_Grouping_2'],
            'Very Detailed Grouping': very_detailed_df['Merged_Grouping_3']
        })
        
        comparison_df.to_csv("comparison_all_groupings.csv", index=False)
        print(f"Saved full comparison to comparison_all_groupings.csv")
        
        # Print comparison
        print("\nComparison of all grouping levels:")
        for _, row in comparison_df.iterrows():
            print(f"Cluster {row['Cluster']} ({row['Cell Type']} / {row['Subtype']}):")
            print(f"  - Broad:        {row['Broad Grouping']}")
            print(f"  - Detailed:     {row['Detailed Grouping']}")
            print(f"  - Very Detailed: {row['Very Detailed Grouping']}")
            print()
            
    except Exception as e:
        print(f"Error with very detailed groupings: {str(e)}")
        
    # Update cleanup files list
    cleanup_files = [
        test_csv,
        "broad_groupings.csv",
        "broad_groupings_readable.csv",
        "detailed_groupings.csv", 
        "detailed_groupings_readable.csv",
        "very_detailed_groupings.csv",
        "very_detailed_groupings_readable.csv",
        "comparison_groupings.csv"
    ]
    
    # Keep the full comparison file for reference, remove others
    for file in cleanup_files:
        if os.path.exists(file):
            try:
                os.remove(file)
                print(f"Removed test file: {file}")
            except:
                pass
    
    print(f"\nFull comparison file kept for reference: comparison_all_groupings.csv")
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 
