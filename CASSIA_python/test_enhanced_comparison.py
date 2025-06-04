#!/usr/bin/env python3
"""
Test file for the enhanced cell type comparison functionality.

This script demonstrates:
1. Structured output with <celltype>, <reasoning>, and <score> tags
2. Score and reasoning extraction
3. HTML report generation
4. CSV output with structured data
5. Multiple test scenarios

Usage:
    python test_enhanced_comparison.py

Requirements:
    - OPENROUTER_API_KEY environment variable (for live API tests)
    - CASSIA package installed
"""

import os
import sys
import traceback
from pathlib import Path


def test_extraction_functionality():
    """Test the score extraction function with various mock responses."""
    print("\n" + "="*60)
    print("🔍 TESTING SCORE EXTRACTION FUNCTIONALITY")
    print("="*60)
    
    try:
        from CASSIA.cell_type_comparison import extract_celltype_scores
        
        # Test Case 1: Perfect structured response
        print("\n📋 Test Case 1: Perfect Structured Response")
        print("-" * 40)
        
        mock_response_perfect = """
<celltype>T cell</celltype>
<reasoning>
The marker set contains CD3, CD4, and CD8 which are classic T cell markers. 
CD3 is found on all T cells, while CD4 and CD8 are specific to helper and 
cytotoxic T cells respectively. This marker combination strongly indicates 
T cell population.
</reasoning>
<score>95</score>

<celltype>B cell</celltype>
<reasoning>
CD19 and CD20 are B cell specific markers present in the set. However, the 
presence of T cell markers (CD3, CD4, CD8) suggests this might be a mixed 
population or the markers are from different cell types in the sample.
</reasoning>
<score>70</score>

<celltype>NK cell</celltype>
<reasoning>
NK cells typically express CD16 and CD56, but CD56 is not in this marker set. 
CD16 is present but it's also found on some monocytes. The lack of NK-specific 
markers makes this less likely to be an NK cell population.
</reasoning>
<score>30</score>
"""
        
        celltypes = ["T cell", "B cell", "NK cell"]
        extracted = extract_celltype_scores(mock_response_perfect, celltypes)
        
        print("✅ Extracted results:")
        for celltype, data in extracted.items():
            score = data['score']
            reasoning = data['reasoning'][:100] + "..." if len(data['reasoning']) > 100 else data['reasoning']
            print(f"  📊 {celltype}: Score {score}")
            print(f"     💭 Reasoning: {reasoning}")
            print()
        
        # Test Case 2: Malformed response (fallback testing)
        print("\n📋 Test Case 2: Malformed Response (Fallback)")
        print("-" * 40)
        
        mock_response_malformed = """
T cell seems to match well based on CD3, CD4, CD8 markers. Score: 85
B cell has some markers like CD19 but mixed with T cell markers. Score: 60
NK cell doesn't match well, only CD16 present. Score: 25
"""
        
        extracted_fallback = extract_celltype_scores(mock_response_malformed, celltypes)
        print("✅ Fallback extraction results:")
        for celltype, data in extracted_fallback.items():
            print(f"  📊 {celltype}: Score {data['score']}, Reasoning: {data['reasoning'][:50]}...")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in extraction testing: {e}")
        traceback.print_exc()
        return False


def test_html_generation():
    """Test HTML report generation with mock data."""
    print("\n" + "="*60)
    print("🎨 TESTING HTML REPORT GENERATION")
    print("="*60)
    
    try:
        from CASSIA.cell_type_comparison import generate_comparison_html_report
        
        # Create mock results data
        mock_results = [
            {
                'model': 'anthropic/claude-3.5-sonnet',
                'tissue': 'blood',
                'species': 'human',
                'cell_types': 'T cell, B cell, NK cell',
                'response': 'Mock response from Claude showing structured analysis...',
                'extracted_scores': {
                    'T cell': {'score': '95', 'reasoning': 'Strong T cell markers CD3, CD4, CD8 present. High confidence match.'},
                    'B cell': {'score': '65', 'reasoning': 'B cell markers CD19, CD20 present but mixed with T cell markers.'},
                    'NK cell': {'score': '25', 'reasoning': 'Limited NK markers, only CD16 present without CD56.'}
                },
                'status': 'success'
            },
            {
                'model': 'openai/gpt-4o',
                'tissue': 'blood', 
                'species': 'human',
                'cell_types': 'T cell, B cell, NK cell',
                'response': 'Mock response from GPT-4 with different analysis...',
                'extracted_scores': {
                    'T cell': {'score': '90', 'reasoning': 'CD3+ cells with CD4/CD8 subsets clearly indicate T cell lineage.'},
                    'B cell': {'score': '70', 'reasoning': 'CD19+CD20+ markers suggest B cell presence despite T cell markers.'},
                    'NK cell': {'score': '35', 'reasoning': 'CD16 expression could indicate NK cells but lacks definitive markers.'}
                },
                'status': 'success'
            },
            {
                'model': 'google/gemini-pro',
                'tissue': 'blood',
                'species': 'human', 
                'cell_types': 'T cell, B cell, NK cell',
                'response': 'Mock response from Gemini with third perspective...',
                'extracted_scores': {
                    'T cell': {'score': '88', 'reasoning': 'T cell receptor complex (CD3) and co-receptors (CD4, CD8) strongly support T cell identity.'},
                    'B cell': {'score': '60', 'reasoning': 'B cell lineage markers present but contamination with T cell markers reduces confidence.'},
                    'NK cell': {'score': '40', 'reasoning': 'Some evidence for NK cells via CD16 but insufficient for definitive classification.'}
                },
                'status': 'success'
            }
        ]
        
        # Generate HTML report
        output_file = "test_comparison_report.html"
        html_content = generate_comparison_html_report(mock_results, output_file)
        
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            print(f"✅ HTML report generated successfully!")
            print(f"📁 File: {output_file}")
            print(f"📏 Size: {file_size:,} bytes")
            print(f"🌐 Open in browser to view the interactive report")
            
            # Show a snippet of the HTML
            print(f"\n📄 HTML Preview (first 200 chars):")
            print(html_content[:200] + "...")
            
            return output_file
        else:
            print("❌ HTML file was not created")
            return None
            
    except Exception as e:
        print(f"❌ Error in HTML generation: {e}")
        traceback.print_exc()
        return None


def test_live_api_comparison():
    """Test the full comparison functionality with live API calls."""
    print("\n" + "="*60)
    print("🚀 TESTING LIVE API COMPARISON")
    print("="*60)
    
    # Check for API key
    api_key = os.environ.get('OPENROUTER_API_KEY')
    if not api_key:
        print("⚠️  OPENROUTER_API_KEY not found in environment variables")
        print("💡 To test live API functionality, set your API key:")
        print("   export OPENROUTER_API_KEY='your_key_here'")
        print("   # or in Python: os.environ['OPENROUTER_API_KEY'] = 'your_key'")
        return False
    
    try:
        from CASSIA.cell_type_comparison import compareCelltypes
        
        print("🔑 API Key found - proceeding with live test")
        print("\n📋 Test Parameters:")
        
        # Test parameters
        tissue = "peripheral blood"
        species = "human"
        celltypes = ["T cell", "B cell", "Monocyte"]
        marker_set = "CD3, CD4, CD8A, CD19, CD20, CD14, CD16, FCGR3A"
        
        print(f"  🧬 Tissue: {tissue}")
        print(f"  🐭 Species: {species}")
        print(f"  🎯 Cell Types: {celltypes}")
        print(f"  🧪 Markers: {marker_set}")
        
        # Use a faster model for testing
        test_models = ["anthropic/claude-3.5-haiku"]
        
        print(f"\n🤖 Using models: {test_models}")
        print("🔄 Running comparison...")
        
        # Run the comparison
        results = compareCelltypes(
            tissue=tissue,
            celltypes=celltypes,
            marker_set=marker_set,
            species=species,
            model_list=test_models,
            output_file="live_test_comparison.csv",
            generate_html_report=True
        )
        
        print("\n✅ Live API test completed successfully!")
        
        # Display results
        print("\n📊 RESULTS SUMMARY:")
        print("-" * 30)
        
        for result in results['results']:
            model = result['model']
            status = result['status']
            extracted_scores = result.get('extracted_scores', {})
            
            print(f"\n🤖 Model: {model} ({status})")
            
            if extracted_scores:
                for celltype, data in extracted_scores.items():
                    score = data.get('score', 'N/A')
                    reasoning = data.get('reasoning', 'N/A')
                    reasoning_preview = reasoning[:80] + "..." if len(reasoning) > 80 else reasoning
                    
                    print(f"  📈 {celltype}: Score {score}")
                    print(f"     💭 {reasoning_preview}")
            else:
                print("  ⚠️  No structured scores extracted")
        
        # File outputs
        csv_file = results.get('csv_file')
        html_file = csv_file.replace('.csv', '_report.html') if csv_file else None
        
        print(f"\n📁 Generated Files:")
        if csv_file and os.path.exists(csv_file):
            print(f"  📊 CSV: {csv_file} ({os.path.getsize(csv_file):,} bytes)")
        if html_file and os.path.exists(html_file):
            print(f"  🌐 HTML: {html_file} ({os.path.getsize(html_file):,} bytes)")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in live API test: {e}")
        traceback.print_exc()
        return False


def test_edge_cases():
    """Test edge cases and error handling."""
    print("\n" + "="*60)
    print("🧪 TESTING EDGE CASES")
    print("="*60)
    
    try:
        from CASSIA.cell_type_comparison import compareCelltypes, extract_celltype_scores
        
        # Test Case 1: Invalid number of cell types
        print("\n📋 Test Case 1: Invalid Cell Type Count")
        print("-" * 40)
        
        try:
            compareCelltypes("tissue", ["only_one"], "markers", generate_html_report=False)
            print("❌ Should have raised ValueError for insufficient cell types")
        except ValueError as e:
            print(f"✅ Correctly caught error: {e}")
        
        # Test Case 2: Empty extraction
        print("\n📋 Test Case 2: Empty Response Extraction")
        print("-" * 40)
        
        empty_response = "No structured data here"
        result = extract_celltype_scores(empty_response, ["T cell", "B cell"])
        print(f"✅ Empty extraction result: {result}")
        
        # Test Case 3: Partial extraction
        print("\n📋 Test Case 3: Partial Response Extraction") 
        print("-" * 40)
        
        partial_response = """
<celltype>T cell</celltype>
<reasoning>Good T cell markers</reasoning>
<score>85</score>

<celltype>B cell</celltype>
<reasoning>Some B cell evidence</reasoning>
<!-- Missing score tag -->
"""
        
        partial_result = extract_celltype_scores(partial_response, ["T cell", "B cell", "NK cell"])
        print("✅ Partial extraction results:")
        for ct, data in partial_result.items():
            print(f"  {ct}: Score='{data['score']}', Reasoning length={len(data['reasoning'])}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in edge case testing: {e}")
        traceback.print_exc()
        return False


def cleanup_test_files():
    """Clean up generated test files."""
    test_files = [
        "test_comparison_report.html",
        "live_test_comparison.csv", 
        "live_test_comparison_report.html"
    ]
    
    print("\n🧹 Cleaning up test files...")
    for file in test_files:
        if os.path.exists(file):
            os.remove(file)
            print(f"  🗑️  Removed {file}")


def main():
    """Run all tests."""
    print("🧪 ENHANCED CELL TYPE COMPARISON - COMPREHENSIVE TEST SUITE")
    print("=" * 70)
    print("This script tests the enhanced cell type comparison functionality")
    print("including structured output, extraction, and HTML report generation.")
    
    test_results = []
    
    # Run all tests
    tests = [
        ("Score Extraction", test_extraction_functionality),
        ("HTML Generation", test_html_generation),
        ("Edge Cases", test_edge_cases),
        ("Live API", test_live_api_comparison)
    ]
    
    for test_name, test_func in tests:
        try:
            print(f"\n🔄 Running {test_name} test...")
            result = test_func()
            test_results.append((test_name, result))
            if result:
                print(f"✅ {test_name} test PASSED")
            else:
                print(f"⚠️  {test_name} test completed with warnings")
        except Exception as e:
            print(f"❌ {test_name} test FAILED: {e}")
            test_results.append((test_name, False))
    
    # Summary
    print("\n" + "="*70)
    print("📋 TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for _, result in test_results if result)
    total = len(test_results)
    
    for test_name, result in test_results:
        status = "✅ PASS" if result else "❌ FAIL/WARN"
        print(f"  {status} {test_name}")
    
    print(f"\n🎯 Overall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests completed successfully!")
    else:
        print("⚠️  Some tests had issues - check the output above")
    
    # Show usage example
    print(f"\n💡 USAGE EXAMPLE:")
    print("=" * 30)
    print("""
import CASSIA

# Enhanced comparison with structured output and HTML report
results = CASSIA.compareCelltypes(
    tissue="blood",
    celltypes=["T cell", "B cell", "NK cell"],
    marker_set="CD3, CD4, CD8, CD19, CD20, CD16, CD56",
    species="human",
    generate_html_report=True  # Creates beautiful HTML visualization
)

# Access structured results
for result in results['results']:
    model = result['model']
    scores = result['extracted_scores']
    for celltype, data in scores.items():
        print(f"{model} - {celltype}: {data['score']}")
""")
    
    cleanup_test_files()
    print("\n✨ Test suite completed!")


if __name__ == "__main__":
    main() 