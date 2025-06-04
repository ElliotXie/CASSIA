#!/usr/bin/env python3
"""
Simple example demonstrating the enhanced cell type comparison functionality.

This example shows how to:
1. Compare multiple cell types using structured output
2. Generate beautiful HTML reports
3. Extract scores and reasoning from AI models

Requirements:
- Set OPENROUTER_API_KEY environment variable
- Install CASSIA package

Usage:
    python example_enhanced_comparison.py
"""

import os
import CASSIA

def run_comparison_example():
    """Run a simple comparison example."""
    
    print("🧬 Enhanced Cell Type Comparison Example")
    print("=" * 50)
    
    # Check for API key
    if not os.environ.get('OPENROUTER_API_KEY'):
        print("❌ OPENROUTER_API_KEY not set!")
        print("\n💡 To run this example, set your API key:")
        print("   os.environ['OPENROUTER_API_KEY'] = 'your_key_here'")
        print("   # or: export OPENROUTER_API_KEY='your_key_here'")
        return
    
    # Example parameters
    print("📋 Analysis Parameters:")
    tissue = "peripheral blood"
    species = "human"
    celltypes = ["T cell", "B cell", "NK cell", "Monocyte"]
    marker_set = "CD3, CD4, CD8, CD19, CD20, CD56, CD16, CD14"
    
    print(f"  🧬 Tissue: {tissue}")
    print(f"  🐭 Species: {species}")
    print(f"  🎯 Cell Types: {celltypes}")
    print(f"  🧪 Markers: {marker_set}")
    
    print(f"\n🚀 Running enhanced comparison...")
    
    try:
        # Run the enhanced comparison
        results = CASSIA.compareCelltypes(
            tissue=tissue,
            celltypes=celltypes,
            marker_set=marker_set,
            species=species,
            model_list=["anthropic/claude-3.5-haiku"],  # Fast model for demo
            output_file="example_results.csv",
            generate_html_report=True  # This is the new feature!
        )
        
        print("✅ Comparison completed!")
        
        # Display structured results
        print(f"\n📊 Structured Results:")
        print("-" * 30)
        
        for result in results['results']:
            model = result['model']
            extracted_scores = result['extracted_scores']
            
            print(f"\n🤖 Model: {model}")
            for celltype, data in extracted_scores.items():
                score = data['score']
                reasoning = data['reasoning'][:100] + "..." if len(data['reasoning']) > 100 else data['reasoning']
                
                print(f"  📈 {celltype}: {score}/100")
                print(f"     💭 {reasoning}")
                print()
        
        # Show generated files
        csv_file = results['csv_file']
        html_file = csv_file.replace('.csv', '_report.html')
        
        print(f"📁 Generated Files:")
        print(f"  📊 CSV Report: {csv_file}")
        print(f"  🌐 HTML Report: {html_file}")
        print(f"\n💡 Open the HTML file in your browser for an interactive report!")
        
    except Exception as e:
        print(f"❌ Error: {e}")


def show_extraction_example():
    """Show how the extraction works with mock data."""
    print(f"\n🔍 Extraction Example (Mock Data)")
    print("=" * 40)
    
    # Import the extraction function
    from CASSIA.cell_type_comparison import extract_celltype_scores
    
    # Mock response in the structured format
    mock_response = """
<celltype>T cell</celltype>
<reasoning>
CD3 is the definitive T cell marker, present on all T cells. CD4 and CD8 
identify T helper and cytotoxic T cell subsets respectively. This combination 
provides strong evidence for T cell identification.
</reasoning>
<score>95</score>

<celltype>B cell</celltype>
<reasoning>
CD19 and CD20 are classic B cell markers. However, the presence of T cell 
markers suggests this might be a mixed population or contamination.
</reasoning>
<score>70</score>

<celltype>NK cell</celltype>
<reasoning>
CD56 is present which is the primary NK cell marker. CD16 is also found on 
NK cells. This provides good evidence for NK cell presence.
</reasoning>
<score>85</score>
"""
    
    celltypes = ["T cell", "B cell", "NK cell"]
    extracted = extract_celltype_scores(mock_response, celltypes)
    
    print("✅ Extracted from structured response:")
    for celltype, data in extracted.items():
        print(f"  🎯 {celltype}: Score {data['score']}")
        reasoning_short = data['reasoning'][:60] + "..." if len(data['reasoning']) > 60 else data['reasoning']
        print(f"     📝 {reasoning_short}")
        print()


def main():
    """Main function."""
    print("🎉 CASSIA Enhanced Cell Type Comparison")
    print("This example demonstrates the new structured output and HTML visualization features.")
    print()
    
    # Show extraction example (works without API key)
    show_extraction_example()
    
    # Run live comparison (requires API key)
    run_comparison_example()
    
    print("\n✨ Example completed!")
    print("\n📚 Key Features Demonstrated:")
    print("  ✅ Structured prompt with <celltype>, <reasoning>, <score> tags")
    print("  ✅ Automatic extraction of scores and reasoning")
    print("  ✅ Beautiful HTML report generation")
    print("  ✅ Enhanced CSV output with structured data")
    print("  ✅ Multiple model support and comparison")


if __name__ == "__main__":
    main() 