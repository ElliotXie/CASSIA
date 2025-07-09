#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example: Using Total Summary Report Style with Custom Models

This script demonstrates how to use the new total_summary report style 
with custom model providers like DeepSeek for CASSIA annotation boost.

The total_summary style creates gene-focused reports that are more readable
and organize findings by biological function rather than iteration sequence.
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from CASSIA.annotation_boost import runCASSIA_annotationboost
import pandas as pd

def example_total_summary_with_deepseek():
    """
    Example showing how to use total_summary report style with DeepSeek custom model.
    """
    print("=== Example: Total Summary Report Style with DeepSeek ===")
    
    # Set up DeepSeek API key (replace with your actual key)
    api_key = "sk-afb39114f1334ba486505d9425937d16"  # Replace with your key
    os.environ["CUSTERMIZED_API_KEY"] = api_key
    
    # Example paths (replace with your actual data paths)
    full_result_path = "data/example_results.csv"  # Your results CSV
    marker_data_path = "data/processed.csv"        # Your marker data
    
    # Load marker data
    try:
        marker_data = pd.read_csv(marker_data_path)
        print(f"✅ Loaded marker data: {marker_data.shape}")
    except FileNotFoundError:
        print(f"❌ Marker data file not found: {marker_data_path}")
        print("Please update the path to your actual marker data file")
        return
    
    # Run annotation boost with different configurations
    configurations = [
        {
            "name": "DeepSeek + Depth-First + Total Summary",
            "provider": "https://api.deepseek.com",
            "search_strategy": "depth",
            "report_style": "total_summary",
            "description": "Focused analysis with gene-focused summary"
        },
        {
            "name": "DeepSeek + Breadth-First + Total Summary", 
            "provider": "https://api.deepseek.com",
            "search_strategy": "breadth",
            "report_style": "total_summary",
            "description": "Comprehensive analysis with gene-focused summary"
        },
        {
            "name": "OpenRouter + Depth-First + Total Summary",
            "provider": "openrouter",
            "search_strategy": "depth", 
            "report_style": "total_summary",
            "description": "Comparison with standard provider"
        }
    ]
    
    for config in configurations:
        print(f"\n{'-'*60}")
        print(f"Running: {config['name']}")
        print(f"Description: {config['description']}")
        print(f"{'-'*60}")
        
        try:
            result = runCASSIA_annotationboost(
                full_result_path=full_result_path,
                marker=marker_data,
                cluster_name="monocyte",  # Change to your cluster of interest
                major_cluster_info="Human Large Intestine",
                output_name=f"example_{config['search_strategy']}_{config['report_style']}",
                num_iterations=3,  # Reduced for example
                model="deepseek-chat" if "deepseek" in config['provider'] else "google/gemini-2.5-flash-preview",
                provider=config['provider'],
                temperature=0,
                conversation_history_mode="final",
                search_strategy=config['search_strategy'],
                report_style=config['report_style']
            )
            
            if isinstance(result, dict) and result.get('status') == 'success':
                print(f"✅ {config['name']} completed successfully!")
                print(f"   Summary report: {result.get('summary_report_path', 'N/A')}")
                print(f"   Raw conversation: {result.get('raw_text_path', 'N/A')}")
            else:
                print(f"⚠️ {config['name']} completed with issues")
                
        except Exception as e:
            print(f"❌ Error with {config['name']}: {str(e)}")
    
    print(f"\n{'='*60}")
    print("EXAMPLE COMPLETE")
    print("="*60)
    print("Key differences you should see in the reports:")
    print("• Total summary reports organize by gene function, not iteration")
    print("• Gene groups show what was tested and what was learned")
    print("• More readable format for understanding biological conclusions")
    print("• Custom providers (DeepSeek) work seamlessly with gene-focused summaries")

def show_command_line_examples():
    """
    Show command line examples for using total_summary with custom models.
    """
    print("\n=== Command Line Examples ===")
    print("Use these commands to test total_summary with custom models:")
    print()
    
    examples = [
        "# Basic total_summary with DeepSeek:",
        "python CASSIA_python_tutorial.py --step boost --provider https://api.deepseek.com --api_key YOUR_API_KEY --report_style total_summary",
        "",
        "# Depth-first search + total_summary with custom model:",
        "python CASSIA_python_tutorial.py --step boost --search_strategy depth --report_style total_summary --provider https://api.deepseek.com --api_key YOUR_API_KEY",
        "",
        "# Test annotation boost with additional task + total_summary:",
        "python CASSIA_python_tutorial.py --step boost_task --report_style total_summary --provider https://api.deepseek.com --api_key YOUR_API_KEY",
        "",
        "# Compare strategies with gene-focused summaries:",
        "python CASSIA_python_tutorial.py --step compare_strategies --report_style total_summary --provider https://api.deepseek.com --api_key YOUR_API_KEY",
        "",
        "# Test all providers with total_summary:",
        "python CASSIA_python_tutorial.py --step test_boost --report_style total_summary",
        "",
        "# New dedicated test for total_summary:",
        "python CASSIA_python_tutorial.py --step test_total_summary --cluster monocyte",
    ]
    
    for example in examples:
        print(example)

if __name__ == "__main__":
    print("CASSIA Total Summary Report Style Example")
    print("=" * 50)
    
    # Show command line examples
    show_command_line_examples()
    
    # Ask user if they want to run the full example
    print(f"\n{'='*50}")
    response = input("Do you want to run the full example? (y/N): ").strip().lower()
    
    if response in ['y', 'yes']:
        example_total_summary_with_deepseek()
    else:
        print("Example skipped. Use the command line examples above to test the functionality.") 